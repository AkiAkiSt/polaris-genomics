"""G1 Generalization: run the ENTIRE POLARIS pipeline on a second disease (CAD)
with its own textbook controls (SORT1; distal PHACTR1->EDN1). Reuses every engine and
reports the same metrics, demonstrating transfer rather than overfitting to T2D."""
from __future__ import annotations
import sys, json, time, pathlib, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from scipy.special import expit
from polaris import apis, config, utils, motifs as M, integrate as IN
from polaris.molecular import annotate_variant

KW = config.CAD_TRAIT_KEYWORDS


def best_cs(vid):
    light = apis.ot_cs_light_all(vid)
    kept = [r for r in light if r.get("studyType") == "gwas"
            and any(k in (r["study"].get("traitFromSource") or "").lower() for k in KW)]
    if not kept:
        return None, 0
    def sc(cs):
        st = cs["study"]; t = (st.get("traitFromSource") or "").lower()
        pr = 3 if "coronary" in t or "myocard" in t else 1
        return pr * 1e6 + (st.get("nSamples") or 1)
    b = max(kept, key=sc)
    return apis.ot_cs_detail_verified(vid, b["_idx"], b["studyLocusId"]), len(kept)


def main():
    pwms = M.load_pwms(str(config.EXT / "JASPAR2024_CORE_vertebrates_nr.txt"))
    rows, l2g_rows, gene_rows = [], [], []
    t0 = time.time()
    for rsid, locus, eff, control, coding, note in config.CAD_LOCI:
        vid = apis.ot_resolve_variant(rsid)
        if not vid:
            print(f"{locus}: unresolved"); continue
        cs, ncs = best_cs(vid)
        chrom, pos, ref, alt = utils.parse_ot_vid(vid)
        members = (cs["locus"]["rows"] if cs else [{"variant": {"id": vid, "rsIds": [rsid]},
                   "posteriorProbability": None, "is95CredibleSet": None}])
        l2g = sorted((cs["l2GPredictions"]["rows"] if cs else []), key=lambda x: -x["score"])
        best_l2g = (l2g[0]["target"]["approvedSymbol"], l2g[0]["score"]) if l2g else (None, None)
        for p in l2g:
            l2g_rows.append(dict(locus=locus, gene=p["target"]["approvedSymbol"], l2g_score=p["score"]))
        # genes for nearest/distance
        try:
            gj = apis.ens_overlap_genes(str(chrom), pos - 600000, pos + 600000)
        except Exception:
            gj = []
        gl = [dict(locus=locus, gene=g.get("external_name") or g["id"],
                   tss=(g["start"] if g.get("strand", 1) == 1 else g["end"]),
                   dist=abs(pos - (g["start"] if g.get("strand", 1) == 1 else g["end"])))
              for g in gj if g.get("biotype") == "protein_coding"]
        gl.sort(key=lambda d: d["dist"])
        for rank, d in enumerate(gl):
            d["is_nearest"] = int(rank == 0); gene_rows.append(d)
        for mi, m in enumerate(members):
            mv = m["variant"]; mc, mp, mr, ma = utils.parse_ot_vid(mv["id"])
            try:
                f = annotate_variant(str(mc), mp, mr, ma, mv["id"], pwms)
            except Exception:
                f = {}
            f.pop("_motif_records", None)
            f.update(dict(locus=locus, rsid=(mv.get("rsIds") or [None])[0], vid=mv["id"],
                          pip=m.get("posteriorProbability"), is95=m.get("is95CredibleSet"),
                          known_effector=eff, is_control=control, coding_lead=coding,
                          best_l2g_gene=best_l2g[0], best_l2g_score=best_l2g[1], is_lead=(mv["id"] == vid)))
            rows.append(f)
        print(f"  {locus:12s} cs={ncs} members={len(members)} L2G={best_l2g[0]} ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows); df["pip"] = df["pip"].fillna(0.0)
    df.to_csv(config.PROC / "cad_integrated.csv", index=False)

    # ---- metrics (mirror T2D) ----
    Me = {}
    POS, NEG = df.pip > 0.5, df.pip < 0.05
    for c, key in [("phylop_max", "conservation"), ("rarity", "rarity")]:
        a, b = df.loc[POS, c].dropna(), df.loc[NEG, c].dropna()
        if len(a) > 2 and len(b) > 2:
            u, p = mannwhitneyu(a, b)
            Me[f"{key}_auc"] = round(u / (len(a) * len(b)), 3); Me[f"{key}_p"] = float(p)
    fcols = ["motif_motif_disruption", "phylop_max", "phastcons470", "ccre_enhancer", "splice_delta"]
    for c in fcols:
        df[c] = df.get(c, np.nan)
    Mf = (df[fcols].values - np.nanmean(df[fcols].values, 0)) / (np.nanstd(df[fcols].values, 0) + 1e-9)
    y = np.full(len(df), np.nan); y[POS.values] = 1; y[NEG.values] = 0
    em = IN.latent_class_em(Mf, y)
    Me["em_reliability"] = {fcols[i]: round(float(em["weight"][i]), 3) for i in range(len(fcols))}
    # gene top-1 via L2G on controls
    l2gdf = pd.DataFrame(l2g_rows); gdf = pd.DataFrame(gene_rows)
    known = dict(df.groupby("locus").known_effector.first()); ctl = dict(df.groupby("locus").is_control.first())
    ctrl_loci = [L for L in df.locus.unique() if ctl.get(L)]
    hit_l2g = hit_near = n = 0
    for L in ctrl_loci:
        sub = l2gdf[l2gdf.locus == L]
        if len(sub):
            pick = sub.sort_values("l2g_score", ascending=False).iloc[0].gene
            hit_l2g += int(pick == known[L])
        ns = gdf[gdf.locus == L]
        if len(ns):
            hit_near += int(ns.sort_values("dist").iloc[0].gene == known[L])
        n += 1
    Me["gene_top1_L2G"] = round(hit_l2g / max(n, 1), 3)
    Me["gene_top1_nearest"] = round(hit_near / max(n, 1), 3)
    Me["n_control_loci"] = n; Me["n_loci"] = df.locus.nunique(); Me["n_variants"] = len(df)
    # SORT1 + PHACTR1 mechanism
    for rs, lab in [("rs12740374", "SORT1"), ("rs9349379", "PHACTR1")]:
        r = df[df.rsid == rs]
        if len(r):
            Me[f"{lab}_top_TF"] = str(r.iloc[0].get("motif_top_disrupt_tf"))
            Me[f"{lab}_dLLR"] = (round(float(r.iloc[0]["motif_top_disrupt_dLLR"]), 2)
                                 if pd.notna(r.iloc[0].get("motif_top_disrupt_dLLR")) else None)
    json.dump(Me, open(config.TAB / "cad_metrics.json", "w"), indent=2, default=float)
    print("\n=== CAD generalization metrics ==="); print(json.dumps(Me, indent=2, default=float))


if __name__ == "__main__":
    main()
