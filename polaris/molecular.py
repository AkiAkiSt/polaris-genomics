"""Step 3 orchestrator: assemble a transparent, multimodal molecular-effect profile
for each candidate variant. Every channel is auditable (no black-box oracle):

  motif       - information-theoretic TF-binding disruption (polaris.motifs)
  conservation- phyloP100way, phastCons100way, phastCons470way (Zoonomia mammals)
  regulatory  - ENCODE SCREEN cCRE class (enhancer / promoter / CTCF) + Ensembl reg build
  splicing    - transparent donor/acceptor consensus delta (cryptic-splice proxy)
  frequency   - gnomAD allele frequency -> rarity / selection proxy
"""
from __future__ import annotations
import numpy as np
from . import apis, utils, motifs as M

# Compact splice consensus log-odds (donor 9mer MAG|GTRAGT, acceptor core YAG|G)
# frequencies from canonical U2 splice signals; used only as a transparent proxy.
_DONOR = np.array([  # positions -3..+6 ; rows A C G T
    [0.34, 0.36, 0.18, 0.12], [0.60, 0.13, 0.14, 0.13], [0.09, 0.03, 0.80, 0.08],
    [0.00, 0.00, 1.00, 0.00], [0.00, 0.00, 0.00, 1.00], [0.59, 0.03, 0.35, 0.03],
    [0.71, 0.08, 0.12, 0.09], [0.06, 0.05, 0.84, 0.05], [0.15, 0.16, 0.23, 0.46]]).T
_ACC = np.array([  # YAG|G core
    [0.10, 0.30, 0.10, 0.50], [0.05, 0.05, 0.05, 0.85], [1.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 1.00, 0.00], [0.30, 0.20, 0.35, 0.15]]).T
_BG = M.HG38_BG


def _splice_pwm(freq):
    p = (freq + 0.01) / (freq + 0.01).sum(0, keepdims=True)
    return np.log2(p / _BG[:, None])


_DLO = _splice_pwm(_DONOR)
_ALO = _splice_pwm(_ACC)


def _splice_delta(seq_idx, c):
    out = 0.0
    for lo in (_DLO, _ALO):
        w = lo.shape[1]
        for s in range(max(0, c - w + 1), min(c, len(seq_idx) - w) + 1):
            sub = seq_idx[s:s + w]
            if (sub < 0).any():
                continue
            ref = float(lo[sub, np.arange(w)].sum())
            alt_idx = sub.copy()
            # ref vs alt handled by caller swapping; here compute match strength
            out = max(out, ref)
    return out


def conservation(chrom, pos):
    uc = utils.ucsc_chrom(chrom)
    z0 = pos - 1
    res = {}
    try:
        pp = apis.ucsc_phylop(uc, z0 - 8, z0 + 9)
        vals = {d["start"]: d["value"] for d in (pp if isinstance(pp, list) else [])}
        res["phylop"] = vals.get(z0, np.nan)
        res["phylop_w17"] = float(np.nanmean(list(vals.values()))) if vals else np.nan
        res["phylop_max"] = float(np.nanmax(list(vals.values()))) if vals else np.nan
    except Exception:
        res["phylop"] = res["phylop_w17"] = res["phylop_max"] = np.nan
    for track, key in (("phastCons100way", "phastcons"), ("phastCons470way", "phastcons470")):
        try:
            pc = apis.ucsc_track(track, uc, z0, z0 + 1)
            res[key] = pc[0]["value"] if isinstance(pc, list) and pc else np.nan
        except Exception:
            res[key] = np.nan
    return res


def regulatory(chrom, pos):
    uc = utils.ucsc_chrom(chrom); z0 = pos - 1
    out = dict(ccre_class=None, in_ccre=0, ccre_enhancer=0, ccre_promoter=0,
               ccre_ctcf=0, ccre_zscore=np.nan, ens_reg=None)
    try:
        feats = apis.ucsc_ccre(uc, z0 - 250, z0 + 251)
        feats = feats if isinstance(feats, list) else []
        hit = None
        for f in feats:
            if f.get("chromStart", 1e18) <= z0 < f.get("chromEnd", -1):
                hit = f; break
        if hit is None and feats:  # nearest within window
            hit = min(feats, key=lambda f: min(abs(z0 - f.get("chromStart", z0)),
                                               abs(z0 - f.get("chromEnd", z0))))
            out["in_ccre"] = 0
        elif hit is not None:
            out["in_ccre"] = 1
        if hit:
            cc = (hit.get("ccre") or hit.get("encodeLabel") or "").lower()
            out["ccre_class"] = hit.get("ccre")
            out["ccre_enhancer"] = int("els" in cc or "enh" in cc)
            out["ccre_promoter"] = int("pls" in cc or "prom" in cc)
            out["ccre_ctcf"] = int("ctcf" in cc)
            out["ccre_zscore"] = hit.get("zScore", np.nan)
    except Exception:
        pass
    try:
        reg = apis.ens_regulatory(chrom, pos - 1, pos + 1)
        if isinstance(reg, list) and reg:
            out["ens_reg"] = reg[0].get("description") or reg[0].get("feature_type")
    except Exception:
        pass
    return out


def frequency(vid):
    out = dict(gnomad_af=np.nan, rarity=np.nan, ultra_rare=0)
    try:
        g = apis.gnomad_variant(utils.to_gnomad(vid))
        if g and g.get("genome"):
            af = g["genome"].get("af")
            out["gnomad_af"] = af
            if af is not None and af > 0:
                out["rarity"] = float(min(-np.log10(af), 8))
                out["ultra_rare"] = int(af < 1e-4)
            elif af == 0 or af is None:
                out["rarity"] = 8.0; out["ultra_rare"] = 1
    except Exception:
        pass
    return out


def annotate_variant(chrom, pos, ref, alt, vid, pwms, flank=25):
    uc = utils.ucsc_chrom(chrom)
    feat = dict(chrom=chrom, pos=pos, ref=ref, alt=alt, vid=vid)
    # sequence window (0-based half-open), variant at local index = flank
    try:
        seq = apis.ucsc_sequence(uc, pos - 1 - flank, pos - 1 + flank + 1)
    except Exception:
        seq = None
    if seq and len(seq) == 2 * flank + 1:
        mo, recs = M.scan_variant(seq, flank, ref, alt, pwms)
        feat.update({f"motif_{k}": v for k, v in mo.items()})
        feat["_motif_records"] = recs
        si = M.seq_to_idx(seq)
        sref = _splice_delta(si, flank)
        si_alt = si.copy(); si_alt[flank] = M.B2I.get(alt.upper(), -1)
        feat["splice_ref"] = sref
        feat["splice_delta"] = abs(_splice_delta(si_alt, flank) - sref)
    else:
        feat["motif_motif_disruption"] = np.nan
        feat["_motif_records"] = []
    feat.update(conservation(chrom, pos))
    feat.update(regulatory(chrom, pos))
    feat.update(frequency(vid))
    return feat
