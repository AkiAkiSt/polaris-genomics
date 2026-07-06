"""Transparent, information-theoretic TF-binding-disruption engine (POLARIS C1).

For every candidate variant we ask a first-principles biophysical question: does the
single-base substitution create or destroy a transcription-factor binding site?

For each JASPAR position weight matrix we compute, over every placement covering the
variant and on both strands, the log-likelihood ratio LLR(seq) = sum_i log2 p(s_i|i)/b0(s_i).
By the Berg-von Hippel approximation LLR is proportional to binding free energy, so the
allelic difference dLLR = LLR(alt) - LLR(ref) is a calibrated DDG-of-binding proxy.
A motif match is assigned an exact p-value under an order-0 Markov background via
convolution of the per-position score distributions (the TFM-Pvalue construction).

No deep network, no opaque weights: every score is auditable and yields a mechanism
annotation such as "abolishes a FOXA2 motif core (dLLR = -7.3 bits, p_ref = 4e-5)".
"""
from __future__ import annotations
import numpy as np
import pathlib, json

BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}
COMP = {0: 3, 1: 2, 2: 1, 3: 0}  # A<->T, C<->G
HG38_BG = np.array([0.295, 0.205, 0.205, 0.295])  # A C G T


def seq_to_idx(seq):
    return np.array([B2I.get(b.upper(), -1) for b in seq], dtype=np.int64)


def load_jaspar(path):
    """Parse JASPAR-format PFM file -> {matrix_id: (name, counts 4xw)}."""
    motifs = {}
    mid = name = None
    rows = []
    for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if mid and rows:
                motifs[mid] = (name, np.array(rows, float))
            parts = line[1:].split()
            mid, name = parts[0], (parts[1] if len(parts) > 1 else parts[0])
            rows = []
        elif "[" in line:
            nums = line[line.index("[") + 1:line.index("]")].split()
            rows.append([float(x) for x in nums])
    if mid and rows:
        motifs[mid] = (name, np.array(rows, float))
    return motifs


def build_pwm(counts, bg=HG38_BG, pseudo=0.8):
    """counts 4xw (ACGT) -> dict(llr, llr_rc, ic, total_ic, w, maxscore, surv, surv_lo, step)."""
    w = counts.shape[1]
    N = counts.sum(0)
    p = (counts + pseudo * bg[:, None]) / (N + pseudo)  # 4xw probabilities
    llr = np.log2(p / bg[:, None])                       # 4xw bits
    ic = np.sum(p * np.log2(p / bg[:, None]), axis=0)    # per-position information (bits)
    # reverse-complement log-odds
    llr_rc = llr[[3, 2, 1, 0], :][:, ::-1]
    surv, lo, step = _null_survival(llr, bg)
    return dict(llr=llr, llr_rc=llr_rc, ic=ic, total_ic=float(ic.sum()), w=w,
                maxscore=float(llr.max(0).sum()), surv=surv, surv_lo=lo, step=step)


def _null_survival(llr, bg, step=0.05):
    """Exact score survival function under order-0 background by convolution."""
    w = llr.shape[1]
    pmf = np.array([1.0]); offset = 0.0
    for i in range(w):
        v = llr[:, i]
        vmin = v.min()
        k = np.round((v - vmin) / step).astype(int)
        local = np.zeros(k.max() + 1)
        for b in range(4):
            local[k[b]] += bg[b]
        pmf = np.convolve(pmf, local)
        offset += vmin
    surv = np.cumsum(pmf[::-1])[::-1]  # P(score >= grid point)
    return surv, offset, step


def match_pvalue(pwm, score):
    idx = int(np.round((score - pwm["surv_lo"]) / pwm["step"]))
    idx = min(max(idx, 0), len(pwm["surv"]) - 1)
    return float(pwm["surv"][idx])


def _scan_strand(llr, seq_idx, c, w):
    """Best LLR over all placements of width w covering index c. Returns (score, start)."""
    W = seq_idx.shape[0]
    starts = np.arange(max(0, c - w + 1), min(c, W - w) + 1)
    if len(starts) == 0:
        return -np.inf, -1, None
    cols = np.arange(w)
    idxmat = seq_idx[starts[:, None] + cols[None, :]]      # (n_s, w)
    valid = (idxmat >= 0).all(axis=1)
    scores = np.full(len(starts), -np.inf)
    if valid.any():
        gathered = llr[idxmat[valid], cols[None, :]]        # (n_valid, w)
        scores[valid] = gathered.sum(axis=1)
    j = int(np.argmax(scores))
    return float(scores[j]), int(starts[j]), idxmat[j]


def scan_variant(seq, var_local, ref, alt, pwms, ref_pthresh=1e-3, top_k=12):
    """Scan all motifs at one variant.
    seq: window string; var_local: 0-based index of variant in seq; ref/alt: alleles.
    Returns aggregate dict + per-motif records (significant ones)."""
    seq_idx = seq_to_idx(seq)
    ref_ok = (0 <= var_local < len(seq)) and seq[var_local].upper() == ref.upper()
    seq_ref = seq_idx.copy()
    seq_alt = seq_idx.copy()
    if 0 <= var_local < len(seq):
        seq_alt[var_local] = B2I.get(alt.upper(), -1)
        seq_ref[var_local] = B2I.get(ref.upper(), seq_ref[var_local])
    recs = []
    for mid, pwm in pwms.items():
        w = pwm["w"]
        best = None
        for llr, strand in ((pwm["llr"], "+"), (pwm["llr_rc"], "-")):
            sr, st, _ = _scan_strand(llr, seq_ref, var_local, w)
            sa, _, _ = _scan_strand(llr, seq_alt, var_local, w)
            # placement-matched alt score at the ref-best start
            if st >= 0:
                cols = np.arange(w)
                aidx = seq_alt[st:st + w]
                sa_match = float(llr[aidx, cols].sum()) if (aidx >= 0).all() else -np.inf
            else:
                sa_match = -np.inf
            cand = (sr, sa_match, sa, st, strand)
            if best is None or sr > best[0]:
                best = cand
        sr, sa_match, sa_best, st, strand = best
        if not np.isfinite(sr):
            continue
        p_ref = match_pvalue(pwm, sr)
        p_alt = match_pvalue(pwm, sa_best)
        d_disrupt = sa_match - sr            # effect on the reference site (<=0 = loss)
        # information at the hit position
        pos_in_motif = (var_local - st) if strand == "+" else (st + w - 1 - var_local)
        ic_hit = float(pwm["ic"][pos_in_motif]) if 0 <= pos_in_motif < w else 0.0
        recs.append(dict(mid=mid, name=pwm["name"], strand=strand, w=w,
                         llr_ref=sr, llr_alt_site=sa_match, llr_alt_best=sa_best,
                         p_ref=p_ref, p_alt=p_alt, d_disrupt=d_disrupt,
                         d_create=sa_best - sr, ic_hit=ic_hit,
                         frac_max=sr / pwm["maxscore"] if pwm["maxscore"] > 0 else 0.0))
    # aggregate
    sig = [r for r in recs if r["p_ref"] <= ref_pthresh]
    disr = sorted(sig, key=lambda r: r["d_disrupt"])           # most negative first
    crea = sorted([r for r in recs if r["p_alt"] <= ref_pthresh and r["d_create"] > 0],
                  key=lambda r: -r["d_create"])
    top_disrupt = disr[0] if disr else None
    top_create = crea[0] if crea else None
    # headline disruption magnitude: significance-weighted |dLLR|
    mag = 0.0
    if top_disrupt:
        mag = max(mag, -top_disrupt["d_disrupt"] * min(1.0, -np.log10(top_disrupt["p_ref"]) / 4))
    if top_create:
        mag = max(mag, top_create["d_create"] * min(1.0, -np.log10(top_create["p_alt"]) / 4))
    out = dict(
        ref_ok=ref_ok, n_sig=len(sig),
        max_abs_dLLR=float(max([abs(r["d_disrupt"]) for r in sig], default=0.0)),
        motif_disruption=float(mag),
        top_disrupt_tf=(top_disrupt["name"] if top_disrupt else None),
        top_disrupt_dLLR=(top_disrupt["d_disrupt"] if top_disrupt else None),
        top_disrupt_pref=(top_disrupt["p_ref"] if top_disrupt else None),
        top_disrupt_ic=(top_disrupt["ic_hit"] if top_disrupt else None),
        top_create_tf=(top_create["name"] if top_create else None),
        top_create_dLLR=(top_create["d_create"] if top_create else None),
    )
    recs_top = sorted(sig, key=lambda r: -abs(r["d_disrupt"]))[:top_k]
    return out, recs_top


def load_pwms(jaspar_path, cache_npz=None, bg=HG38_BG):
    motifs = load_jaspar(jaspar_path)
    pwms = {}
    for mid, (name, counts) in motifs.items():
        pwm = build_pwm(counts, bg)
        pwm["name"] = name
        pwms[mid] = pwm
    return pwms
