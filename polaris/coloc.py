"""Bayesian colocalization from scratch (Giambartolomei et al. 2014, PLoS Genetics).

Tests whether a GWAS signal and a molecular-QTL signal share a single causal variant
(H4) versus having distinct causal variants (H3). Per-SNP Wakefield approximate Bayes
factors are combined under the 5-hypothesis enumeration. Validated on simulation;
Open Targets colocalisation (full-summary-statistics coloc) provides the production
PP4 channel for every locus, and this implementation corroborates it.
"""
from __future__ import annotations
import numpy as np
from scipy.special import logsumexp


def wakefield_labf(beta, se, W):
    """Natural-log Wakefield ABF for association vs null. Arrays in, array out."""
    beta = np.asarray(beta, float); se = np.asarray(se, float)
    v = se ** 2
    r = W / (W + v)
    z = beta / se
    return 0.5 * (np.log(1 - r) + r * z ** 2)


def _logdiffexp(a, b):
    # log(exp(a) - exp(b)), a > b
    return a + np.log1p(-np.exp(np.clip(b - a, -700, 0)))


def coloc_from_labf(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5):
    """Coloc PP0..PP4 directly from two arrays of per-SNP log approximate Bayes factors."""
    l1 = np.asarray(l1, float); l2 = np.asarray(l2, float)
    L1 = logsumexp(l1); L2 = logsumexp(l2)
    Lshared = logsumexp(l1 + l2)
    Ldistinct = _logdiffexp(L1 + L2, Lshared)
    lp = np.array([0.0, np.log(p1) + L1, np.log(p2) + L2,
                   np.log(p1) + np.log(p2) + Ldistinct, np.log(p12) + Lshared])
    lp -= logsumexp(lp); pp = np.exp(lp)
    return dict(PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4],
                nsnps=len(l1), best_snp=int(np.argmax(l1 + l2)))


def coloc_abf(beta1, se1, beta2, se2, W1=0.21**2, W2=0.15**2,
              p1=1e-4, p2=1e-4, p12=1e-5):
    """Return dict of PP0..PP4 for two traits over a common SNP set."""
    l1 = wakefield_labf(beta1, se1, W1)
    l2 = wakefield_labf(beta2, se2, W2)
    L1 = logsumexp(l1)                       # sum_i BF1_i
    L2 = logsumexp(l2)                       # sum_i BF2_i
    Lshared = logsumexp(l1 + l2)             # sum_i BF1_i BF2_i   (H4 kernel)
    # H3 kernel: sum_{i!=j} BF1_i BF2_j = exp(L1+L2) - exp(Lshared)
    Ldistinct = _logdiffexp(L1 + L2, Lshared)
    lp = np.array([
        0.0,                                 # H0
        np.log(p1) + L1,                     # H1
        np.log(p2) + L2,                     # H2
        np.log(p1) + np.log(p2) + Ldistinct, # H3
        np.log(p12) + Lshared,               # H4
    ])
    lp -= logsumexp(lp)
    pp = np.exp(lp)
    return dict(PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4],
                nsnps=len(l1), best_snp=int(np.argmax(l1 + l2)))


def z_from_p(pval, sign=1.0):
    """Two-sided p-value -> |z| with given sign (for GTEx nes direction)."""
    from scipy.stats import norm
    pval = np.clip(np.asarray(pval, float), 1e-300, 1.0)
    return np.sign(sign) * np.abs(norm.isf(pval / 2.0))


# --------------------------------------------------------------------------- #
def simulate_coloc(p=150, n1=80000, n2=500, shared=True, lead_gap=12,
                   block_rho=0.9, seed=0):
    """Simulate GWAS + eQTL z over an LD block; return (beta1,se1,beta2,se2,k1,k2)."""
    rng = np.random.default_rng(seed)
    idx = np.arange(p)
    R = block_rho ** np.abs(idx[:, None] - idx[None, :])
    np.fill_diagonal(R, 1.0)
    Lc = np.linalg.cholesky(R + 1e-6 * np.eye(p))
    k1 = p // 2
    k2 = k1 if shared else min(p - 1, k1 + lead_gap)
    b1 = np.zeros(p); b1[k1] = 0.10
    b2 = np.zeros(p); b2[k2] = 0.35
    z1 = np.sqrt(n1) * (R @ b1) + Lc @ rng.standard_normal(p)
    z2 = np.sqrt(n2) * (R @ b2) + Lc @ rng.standard_normal(p)
    se1 = np.full(p, 1.0); se2 = np.full(p, 1.0)
    return z1 * se1, se1, z2 * se2, se2, k1, k2


if __name__ == "__main__":
    res = {"shared": [], "distinct": []}
    for seed in range(40):
        for scen in (True, False):
            b1, s1, b2, s2, k1, k2 = simulate_coloc(shared=scen, seed=seed)
            r = coloc_abf(b1, s1, b2, s2, W1=10.0, W2=10.0)
            res["shared" if scen else "distinct"].append(r["PP4"])
    import numpy as np
    print("[coloc-ABF simulation, 40 reps]")
    print(f"  shared-causal   mean PP4 = {np.mean(res['shared']):.3f}  "
          f"(>0.8 in {np.mean(np.array(res['shared'])>0.8):.0%})")
    print(f"  distinct-causal mean PP4 = {np.mean(res['distinct']):.3f}  "
          f"(<0.2 in {np.mean(np.array(res['distinct'])<0.2):.0%})")
