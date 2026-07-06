"""Real-data colocalization (from-scratch Giambartolomei) on showcase T2D loci.

GWAS evidence = per-variant fine-mapping log-Bayes-factors (Open Targets credible set);
eQTL evidence = Wakefield ABFs from GTEx per-variant effect sizes (dyneqtl). Combined by
the 5-hypothesis enumeration to a colocalization PP4 and a tissue of action. The engine
is validated on simulation; this corroborates the production OT colocalisation.
"""
from __future__ import annotations
import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from polaris import apis, config, coloc, utils

SHOWCASE = [("IGF2BP2", "rs4402960", "IGF2BP2"), ("CDKAL1", "rs7754840", "CDKAL1"),
            ("KLF14", "rs972283", "KLF14"), ("THADA", "rs7578597", "THADA"),
            ("JAZF1", "rs864745", "JAZF1"), ("TP53INP1", "rs896854", "TP53INP1"),
            ("HHEX_IDE", "rs1111875", "HHEX"), ("ARAP1_CENTD2", "rs1552224", "STARD10"),
            ("CCND2", "rs11063069", "CCND2"), ("HMG20A", "rs7177055", "HMG20A")]
TISSUES = ["Pancreas", "Adipose_Subcutaneous", "Liver", "Muscle_Skeletal", "Whole_Blood"]
KW = config.T2D_TRAIT_KEYWORDS
MAXM = 30
W2 = 0.04  # eQTL prior variance


def best_t2d_cs(vid):
    light = apis.ot_cs_light_all(vid)
    kept = [r for r in light if r.get("studyType") == "gwas"
            and any(k in (r["study"].get("traitFromSource") or "").lower() for k in KW)]
    if not kept:
        return None
    def sc(cs):
        st = cs["study"]; t = (st.get("traitFromSource") or "").lower()
        return (3 if "type 2 diab" in t else (2 if "diab" in t else 1)) * 1e6 + (st.get("nSamples") or 1)
    b = max(kept, key=sc)
    return apis.ot_cs_detail_verified(vid, b["_idx"], b["studyLocusId"])


def main():
    rows, figdata = [], {}
    for locus, rsid, gene in SHOWCASE:
        vid = apis.ot_resolve_variant(rsid)
        cs = best_t2d_cs(vid)
        if not cs:
            print(f"{locus:14s} no T2D credible set"); continue
        members = sorted(cs["locus"]["rows"], key=lambda m: -(m.get("posteriorProbability") or 0))[:MAXM]
        gw = {}
        for m in members:
            mv = m["variant"]; c, p, r, a = utils.parse_ot_vid(mv["id"])
            # coloc-SuSiE insight: within a credible set PIP ~ ABF, so log(PIP) is the
            # per-variant GWAS log approximate Bayes factor (Wallace 2021).
            lbf = m.get("logBF")
            pip = m.get("posteriorProbability")
            if lbf is None and pip is not None and pip > 1e-12:
                lbf = float(np.log(pip))
            elif lbf is None and m.get("beta") is not None and m.get("standardError"):
                z = m["beta"] / m["standardError"]; lbf = 0.5 * (np.log(1 / 11) + z * z * 10 / 11)
            if lbf is not None:
                gw[mv["id"]] = dict(labf=float(lbf), pos=p, rsid=(mv.get("rsIds") or [None])[0], pip=pip)
        gobj = apis.gtex_gene(gene)
        if not gobj:
            print(f"{locus:14s} no GTEx gene {gene}"); continue
        gid = gobj["gencodeId"]; best = None
        for tissue in TISSUES:
            l1, l2, pos, rs, pips = [], [], [], [], []
            for vidm, gmeta in gw.items():
                try:
                    e = apis.gtex_dyneqtl(utils.to_gtex(vidm), gid, tissue)
                except Exception:
                    e = None
                if not e or e.get("nes") is None or not e.get("error"):
                    continue
                l1.append(gmeta["labf"])
                l2.append(float(coloc.wakefield_labf(e["nes"], max(e["error"], 1e-3), W2)))
                pos.append(gmeta["pos"]); rs.append(gmeta["rsid"]); pips.append(gmeta["pip"])
            if len(l1) < 4:
                continue
            res = coloc.coloc_from_labf(np.array(l1), np.array(l2))
            rows.append(dict(locus=locus, gene=gene, tissue=tissue, n_common=len(l1),
                             PP3=round(res["PP3"], 4), PP4=round(res["PP4"], 4)))
            if best is None or res["PP4"] > best[0]:
                best = (res["PP4"], tissue, len(l1))
                figdata[locus] = dict(gene=gene, tissue=tissue, pos=pos, rsids=rs, pip=pips,
                                      gwas_labf=l1, eqtl_labf=l2, PP4=round(res["PP4"], 4))
        if best:
            print(f"{locus:14s} -> {gene:8s} best={best[1]:20s} PP4={best[0]:.3f} (n={best[2]})")
        else:
            print(f"{locus:14s} -> {gene:8s} insufficient overlap")
    pd.DataFrame(rows).to_csv(config.TAB / "real_coloc.csv", index=False)
    json.dump(figdata, open(config.PROC / "coloc_figdata.json", "w"))
    print(f"\nsaved real_coloc.csv ({len(rows)} tests, {len(figdata)} loci colocalized)")


if __name__ == "__main__":
    main()
