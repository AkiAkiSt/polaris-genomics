"""SuSiE-RSS fine-mapping, implemented from scratch (Wang et al. 2020 JRSS-B;
Zou et al. 2022 PLoS Genet). Sum of Single Effects regression on z-scores + LD.

Used as an independent, transparent fine-mapper: validated on simulation (recovers
causal variants with calibrated credible sets) and applied to real loci to corroborate
the Open Targets credible sets that feed the POLARIS posterior.
"""
from __future__ import annotations
import numpy as np


def _ser(r, R_unused, V, pi, est_prior=True, n_inner=3):
    """Single-effect regression in z-space. r=residual z (p,), V=prior var.
    Returns alpha (p,), mu (p,), post_var (p,), V, lbf (p,)."""
    p = r.shape[0]
    shat2 = 1.0
    for _ in range(n_inner if est_prior else 1):
        post_var = (V * shat2) / (V + shat2)
        post_mean = post_var * r / shat2
        lbf = 0.5 * np.log(shat2 / (shat2 + V)) + 0.5 * (r ** 2) * (V / (shat2 * (shat2 + V)))
        lbf = np.clip(lbf, -700, 700)
        w = np.log(pi) + lbf
        w -= w.max()
        alpha = np.exp(w); alpha /= alpha.sum()
        if est_prior:
            V_new = float(np.sum(alpha * (post_mean ** 2 + post_var)))
            if not np.isfinite(V_new) or V_new <= 1e-9:
                break
            V = V_new
    post_var = (V * shat2) / (V + shat2)
    post_mean = post_var * r / shat2
    lbf = 0.5 * np.log(shat2 / (shat2 + V)) + 0.5 * (r ** 2) * (V / (shat2 * (shat2 + V)))
    lbf = np.clip(lbf, -700, 700)
    w = np.log(pi) + lbf; w -= w.max()
    alpha = np.exp(w); alpha /= alpha.sum()
    return alpha, post_mean, np.full(p, post_var), V, lbf


def susie_rss(z, R, L=10, max_iter=200, tol=1e-5, V0=0.2, est_prior=True,
              coverage=0.95, min_purity=0.5, pi=None):
    """Fine-map from z-scores (p,) and LD matrix R (p,p, diag=1).
    Returns dict with pip, alpha (L,p), mu (L,p), credible sets, prior vars."""
    z = np.asarray(z, float)
    R = np.asarray(R, float)
    p = z.shape[0]
    if pi is None:
        pi = np.full(p, 1.0 / p)
    alpha = np.zeros((L, p))
    mu = np.zeros((L, p))
    pv = np.zeros((L, p))
    V = np.full(L, float(V0))
    bbar = np.zeros(p)
    hist = []
    for it in range(max_iter):
        a_prev = alpha.copy()
        for l in range(L):
            b_others = bbar - alpha[l] * mu[l]
            r = z - R @ b_others
            alpha[l], mu[l], pv[l], V[l], _ = _ser(r, R, max(V[l], 1e-6), pi, est_prior)
            bbar = b_others + alpha[l] * mu[l]
        delta = np.abs(alpha - a_prev).max()
        hist.append(delta)
        if delta < tol:
            break
    pip = 1.0 - np.prod(1.0 - alpha, axis=0)
    pip = np.clip(pip, 0, 1)

    # credible sets per effect (cumulative coverage, purity-filtered)
    sets = []
    for l in range(L):
        order = np.argsort(-alpha[l])
        cum = np.cumsum(alpha[l][order])
        k = int(np.searchsorted(cum, coverage)) + 1
        idx = order[:k]
        if len(idx) == 0:
            continue
        sub = np.abs(R[np.ix_(idx, idx)])
        purity = sub[np.triu_indices(len(idx), 1)].min() if len(idx) > 1 else 1.0
        if purity >= min_purity and alpha[l].max() > 1.0 / p:
            sets.append(dict(effect=l, idx=idx.tolist(), coverage=float(cum[k - 1]),
                             purity=float(purity), top=int(idx[0]), top_alpha=float(alpha[l][idx[0]])))
    # dedupe sets that point to the same top SNP
    seen, ded = set(), []
    for s in sorted(sets, key=lambda s: -s["top_alpha"]):
        if s["top"] in seen:
            continue
        seen.add(s["top"]); ded.append(s)
    return dict(pip=pip, alpha=alpha, mu=mu, post_var=pv, V=V, n_iter=it + 1,
                converged=(hist[-1] < tol if hist else False), hist=hist, sets=ded)


# --------------------------------------------------------------------------- #
def simulate_locus(p=200, n=50000, n_causal=2, block_rho=0.92, seed=0):
    """Simulate a realistic LD block + z-scores with known causal variants."""
    rng = np.random.default_rng(seed)
    # AR(1)-style banded LD with a couple of sub-blocks
    idx = np.arange(p)
    R = block_rho ** np.abs(idx[:, None] - idx[None, :])
    # add a second correlation ridge to make it non-trivial
    R = 0.85 * R + 0.15 * (block_rho ** np.abs(((idx[:, None] % 50) - (idx[None, :] % 50))))
    np.fill_diagonal(R, 1.0)
    # nearest PD
    w, Vt = np.linalg.eigh(R); w = np.clip(w, 1e-4, None); R = (Vt * w) @ Vt.T
    d = np.sqrt(np.diag(R)); R = R / np.outer(d, d)
    causal = rng.choice(p, n_causal, replace=False)
    b = np.zeros(p)
    b[causal] = rng.choice([-1, 1], n_causal) * rng.uniform(0.06, 0.12, n_causal)
    mean_z = np.sqrt(n) * (R @ b)
    Lc = np.linalg.cholesky(R + 1e-6 * np.eye(p))
    z = mean_z + Lc @ rng.standard_normal(p)
    return z, R, sorted(causal.tolist()), b


if __name__ == "__main__":
    # validation on simulation
    aucs, hits, fdrs = [], [], []
    for seed in range(30):
        z, R, causal, b = simulate_locus(seed=seed)
        out = susie_rss(z, R, L=10)
        pip = out["pip"]
        top_sets = out["sets"]
        # recall: each causal captured by some credible set
        captured = sum(any(c in s["idx"] for s in top_sets) for c in causal)
        hits.append(captured / len(causal))
        # false discovery among sets: sets not containing any causal
        if top_sets:
            fdrs.append(np.mean([0 if any(c in s["idx"] for c in causal) else 1 for s in top_sets]))
        # PIP of causals
        aucs.append(float(np.mean(pip[causal])))
    print(f"[SuSiE-RSS simulation, 30 reps]")
    print(f"  mean PIP at true causals : {np.mean(aucs):.3f}")
    print(f"  causal recall (in a CS)  : {np.mean(hits):.3f}")
    print(f"  credible-set FDR         : {np.mean(fdrs):.3f}")
    z, R, causal, b = simulate_locus(seed=1)
    out = susie_rss(z, R, L=10)
    print(f"  example: causals={causal}, #sets={len(out['sets'])}, "
          f"top SNPs={[s['top'] for s in out['sets']]}, conv={out['converged']} iters={out['n_iter']}")
