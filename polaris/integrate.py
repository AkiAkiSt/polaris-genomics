"""POLARIS evidence-integration engine (Step 5) - the novel core.

Two guardrails from the proposal, made mathematical:

(1) Functional != Pathogenic. The disease-causal-via-gene posterior FACTORIZES:
        P(causal_g | E) = P(F=1 | E_molecular) * P(D_g=1 | E_genetic,linkage)
    so a variant can be highly functional yet not disease-linked, and vice versa.

(2) Triangulate; distrust any single oracle. Each evidence channel is a noisy annotator
    of a latent truth. A semi-supervised Gaussian latent-class model (continuous
    Dawid-Skene; Dawid & Skene 1979) infers, by EM, each channel's class-conditional
    distributions -> its discriminability d' -> its reliability weight. The per-sample
    log-posterior decomposes additively into auditable per-channel contributions.

Calibration (Platt / isotonic; Brier, ECE) turns scores into trustworthy probabilities.
"""
from __future__ import annotations
import numpy as np
from scipy.special import expit, logsumexp

LOG2PI = np.log(2 * np.pi)


def _gauss_ll(x, mu, var):
    return -0.5 * (LOG2PI + np.log(var) + (x - mu) ** 2 / var)


def latent_class_em(X, y_anchor=None, max_iter=300, tol=1e-6, var_floor=1e-3,
                    prior=0.5, ridge=1e-9):
    """Semi-supervised Gaussian naive-Bayes latent-class model.

    X         : (n, k) standardized features, NaN allowed (missing-at-random, omitted).
    y_anchor  : (n,) in {0,1,nan}; nan = unlabeled. Anchors pin class identity.
    Returns dict: p1 (responsibilities), logit, per-channel mu/var/dprime/weight, pi, contrib.
    """
    X = np.asarray(X, float)
    n, k = X.shape
    obs = ~np.isnan(X)
    Xz = np.where(obs, X, 0.0)
    if y_anchor is None:
        y_anchor = np.full(n, np.nan)
    y_anchor = np.asarray(y_anchor, float)
    anc = ~np.isnan(y_anchor)

    # init responsibilities
    with np.errstate(invalid="ignore", divide="ignore"):
        seed = np.nanmean(np.where(obs, X, np.nan), axis=1)
    seed = np.nan_to_num(seed, nan=0.0)
    seed = (seed - seed.mean()) / (seed.std() + 1e-9)
    r1 = expit(seed)
    r1[anc] = y_anchor[anc]
    r1 = np.clip(np.nan_to_num(r1, nan=prior), 1e-6, 1 - 1e-6)
    pi = prior
    mu1 = mu0 = var1 = var0 = None
    prev = None
    for it in range(max_iter):
        # M-step (weighted moments per class, per feature, over observed entries)
        w1 = r1[:, None] * obs
        w0 = (1 - r1)[:, None] * obs
        s1 = w1.sum(0) + ridge
        s0 = w0.sum(0) + ridge
        mu1 = (w1 * Xz).sum(0) / s1
        mu0 = (w0 * Xz).sum(0) / s0
        var1 = np.maximum((w1 * (Xz - mu1) ** 2).sum(0) / s1, var_floor)
        var0 = np.maximum((w0 * (Xz - mu0) ** 2).sum(0) / s0, var_floor)
        pi = np.clip(r1.mean(), 1e-4, 1 - 1e-4)
        # E-step
        ll1 = np.where(obs, _gauss_ll(Xz, mu1, var1), 0.0).sum(1) + np.log(pi)
        ll0 = np.where(obs, _gauss_ll(Xz, mu0, var0), 0.0).sum(1) + np.log(1 - pi)
        logit = ll1 - ll0
        r1 = expit(logit)
        r1[anc] = y_anchor[anc]
        if prev is not None and np.abs(logit - prev).max() < tol:
            break
        prev = logit
    dprime = (mu1 - mu0) / np.sqrt(0.5 * (var1 + var0))
    weight = dprime ** 2
    weight = weight / (weight.sum() + 1e-12)
    # per-channel additive contribution to the logit (auditable)
    contrib = np.where(obs, _gauss_ll(Xz, mu1, var1) - _gauss_ll(Xz, mu0, var0), 0.0)
    return dict(p1=r1, logit=logit, pi=pi, mu1=mu1, mu0=mu0, var1=var1, var0=var0,
                dprime=dprime, weight=weight, contrib=contrib, n_iter=it + 1)


# ----------------------------- calibration -------------------------------- #
def platt(scores, y):
    from sklearn.linear_model import LogisticRegression
    s = np.asarray(scores, float).reshape(-1, 1)
    lr = LogisticRegression().fit(s, y)
    return lambda x: lr.predict_proba(np.asarray(x, float).reshape(-1, 1))[:, 1]


def isotonic(scores, y):
    from sklearn.isotonic import IsotonicRegression
    ir = IsotonicRegression(out_of_bounds="clip").fit(scores, y)
    return lambda x: ir.predict(x)


def brier(p, y):
    return float(np.mean((np.asarray(p) - np.asarray(y)) ** 2))


def ece(p, y, bins=10):
    p = np.asarray(p, float); y = np.asarray(y, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < bins - 1 else p <= edges[i + 1])
        if m.sum():
            e += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(e)


def reliability_curve(p, y, bins=10):
    p = np.asarray(p, float); y = np.asarray(y, float)
    edges = np.linspace(0, 1, bins + 1)
    xs, ys, ns = [], [], []
    for i in range(bins):
        m = (p >= edges[i]) & (p <= edges[i + 1] if i == bins - 1 else p < edges[i + 1])
        if m.sum():
            xs.append(p[m].mean()); ys.append(y[m].mean()); ns.append(int(m.sum()))
    return np.array(xs), np.array(ys), np.array(ns)


# ----------------------------- self-test ---------------------------------- #
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n = 800
    y = rng.binomial(1, 0.4, n)
    # 3 informative channels of differing reliability + 2 pure-noise channels
    reliab = [2.2, 1.1, 0.6, 0.0, 0.0]
    X = np.column_stack([rng.normal(d * y, 1.0) for d in reliab])
    X[rng.random(X.shape) < 0.1] = np.nan  # 10% missing
    # reveal 15% of labels as anchors
    ya = np.full(n, np.nan); idx = rng.choice(n, int(0.15 * n), replace=False); ya[idx] = y[idx]
    out = latent_class_em((X - np.nanmean(X, 0)) / np.nanstd(X, 0), ya)
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(y, out["p1"])
    print("[POLARIS latent-class EM self-test]")
    print("  inferred d' per channel:", np.round(out["dprime"], 2), " (true:", reliab, ")")
    print("  inferred weights       :", np.round(out["weight"], 2))
    print(f"  AUROC(latent truth)    : {auc:.3f}   iters={out['n_iter']}  pi={out['pi']:.2f}")
    cal = isotonic(out["logit"], y)
    p = cal(out["logit"])
    print(f"  Brier={brier(p,y):.3f}  ECE={ece(p,y):.3f}")
