"""Advanced geometry & random-matrix tools for the POLARIS evidence space.

diffusion_map  - Coifman & Lafon (2006). Embeds variants by a Markov diffusion whose
                 transition operator is the heat semigroup e^{-tL} of the graph Laplacian
                 L (a discretization of the heat/diffusion PDE du/dt = -L u). Eigen-
                 decomposition of L yields nonlinear coordinates of the evidence manifold.
sinkhorn       - Cuturi (2013). Entropic-regularized optimal transport; used to attribute
                 a variant's causal probability mass across tissues (Wasserstein geometry).
mp_*           - Marchenko-Pastur (1967) random-matrix null for the oracle correlation
                 matrix: eigenvalues above the MP upper edge are genuine shared signal,
                 quantifying how INDEPENDENT the evidence channels are (triangulation).
"""
from __future__ import annotations
import numpy as np
from scipy.spatial.distance import pdist, squareform


# --------------------------- diffusion maps -------------------------------- #
def diffusion_map(X, n_comp=3, eps=None, alpha=1.0, t=1.0):
    X = np.asarray(X, float)
    D2 = squareform(pdist(X, "sqeuclidean"))
    if eps is None:
        # local scale: median squared distance to the k-th nearest neighbour
        k = min(7, D2.shape[0] - 1)
        knn = np.sort(D2, axis=1)[:, 1:k + 1].mean(axis=1)
        eps = np.median(knn) if knn.size else 1.0
        eps = max(eps, 1e-9)
    K = np.exp(-D2 / eps)
    d = K.sum(1)
    Ka = K / np.outer(d ** alpha, d ** alpha)        # anisotropic (alpha=1 -> Laplace-Beltrami)
    dd = Ka.sum(1)
    s = np.sqrt(dd)
    Ms = Ka / np.outer(s, s)
    Ms = (Ms + Ms.T) / 2.0
    w, V = np.linalg.eigh(Ms)
    idx = np.argsort(-w); w = np.clip(w[idx], 0, None); V = V[:, idx]
    psi = V / s[:, None]
    coords = psi[:, 1:n_comp + 1] * (w[1:n_comp + 1] ** t)
    return dict(coords=coords, evals=w, psi=psi, eps=float(eps),
                spectral_gap=float(w[1] - w[2]) if len(w) > 2 else np.nan)


def diffusion_distance(dm, i, j, t=1.0):
    w = dm["evals"][1:]; psi = dm["psi"][:, 1:]
    return float(np.sqrt(np.sum((w ** (2 * t)) * (psi[i] - psi[j]) ** 2)))


# --------------------------- optimal transport ----------------------------- #
def sinkhorn(a, b, C, reg=0.05, n_iter=1000, tol=1e-9):
    a = np.asarray(a, float); b = np.asarray(b, float); C = np.asarray(C, float)
    a = a / a.sum(); b = b / b.sum()
    K = np.exp(-C / reg)
    u = np.ones_like(a)
    for _ in range(n_iter):
        up = u
        v = b / (K.T @ u + 1e-300)
        u = a / (K @ v + 1e-300)
        if np.max(np.abs(u - up)) < tol:
            break
    P = u[:, None] * K * v[None, :]
    return dict(plan=P, distance=float(np.sum(P * C)), marg_row=P.sum(1), marg_col=P.sum(0))


def wasserstein_barycenter(dists, C, reg=0.05, weights=None, n_iter=400):
    """Entropic Wasserstein barycenter of histograms (columns of `dists`)."""
    dists = np.asarray(dists, float)
    n, m = dists.shape
    weights = np.ones(m) / m if weights is None else np.asarray(weights) / np.sum(weights)
    K = np.exp(-C / reg)
    u = np.ones((n, m)); bary = np.ones(n) / n
    for _ in range(n_iter):
        v = dists / (K.T @ u + 1e-300)
        Ku = K @ v
        bary = np.exp((np.log(Ku + 1e-300) @ weights))
        u = bary[:, None] / (Ku + 1e-300)
    return bary / bary.sum()


# --------------------------- Marchenko-Pastur ------------------------------ #
def mp_edges(gamma, sigma=1.0):
    return sigma ** 2 * (1 - np.sqrt(gamma)) ** 2, sigma ** 2 * (1 + np.sqrt(gamma)) ** 2


def mp_pdf(x, gamma, sigma=1.0):
    lm, lp = mp_edges(gamma, sigma)
    x = np.asarray(x, float); out = np.zeros_like(x)
    m = (x > lm) & (x < lp)
    out[m] = np.sqrt((lp - x[m]) * (x[m] - lm)) / (2 * np.pi * gamma * sigma ** 2 * x[m])
    return out


def oracle_spectrum(Xz):
    """Eigenspectrum of the channel correlation matrix vs MP null.
    Xz: (n, k) standardized channels. Returns evals + count above MP edge."""
    Xz = np.asarray(Xz, float)
    Xz = Xz[~np.isnan(Xz).any(1)]
    n, k = Xz.shape
    Xz = (Xz - Xz.mean(0)) / (Xz.std(0) + 1e-9)
    Cmat = np.corrcoef(Xz, rowvar=False)
    ev = np.sort(np.linalg.eigvalsh(Cmat))[::-1]
    gamma = k / n
    lm, lp = mp_edges(gamma, 1.0)
    return dict(evals=ev, gamma=float(gamma), mp_lo=float(lm), mp_hi=float(lp),
                n_signal=int((ev > lp).sum()), n=n, k=k,
                effective_independence=float(k / max(ev[0], 1.0)))


if __name__ == "__main__":
    from scipy.stats import spearmanr
    rng = np.random.default_rng(0)
    # diffusion map on an OPEN manifold (Archimedean spiral) embedded in noisy 3D
    th = np.sort(rng.uniform(0.5, 4 * np.pi, 300))
    X = np.column_stack([th * np.cos(th), th * np.sin(th), 0.4 * rng.standard_normal(300)])
    dm = diffusion_map(X, n_comp=2)
    rho = spearmanr(th, dm["coords"][:, 0]).correlation
    print(f"[diffusion map] |Spearman(coord1, arc-length)| = {abs(rho):.3f} (eps={dm['eps']:.2f})")
    # sinkhorn vs exact 1D OT (sorted matching)
    a = np.ones(20) / 20; b = np.ones(20) / 20
    xs = np.sort(rng.random(20)); ys = np.sort(rng.random(20))
    Cm = (xs[:, None] - ys[None, :]) ** 2
    sk = sinkhorn(a, b, Cm, reg=0.005)
    exact = np.mean((xs - ys) ** 2)
    print(f"[sinkhorn OT] W2^2 entropic={sk['distance']:.4f}  exact-sorted={exact:.4f}")
    # MP: 3 correlated signal channels + 3 noise among k=6, n=400
    n = 400
    f = rng.standard_normal((n, 2))
    sig = np.column_stack([f[:, 0], f[:, 0] + 0.5 * rng.standard_normal(n), f[:, 1]])
    noise = rng.standard_normal((n, 3))
    Xz = np.column_stack([sig, noise])
    sp = oracle_spectrum(Xz)
    print(f"[Marchenko-Pastur] gamma={sp['gamma']:.3f} MP_hi={sp['mp_hi']:.2f} "
          f"#signal eigenvalues above edge = {sp['n_signal']}  (top evals "
          f"{np.round(sp['evals'][:4],2)})")
