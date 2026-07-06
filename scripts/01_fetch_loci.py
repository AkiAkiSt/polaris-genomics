"""POLARIS Step 1-2: pull T2D GWAS loci + fine-mapped credible sets (Open Targets).

For each curated lead variant we resolve its Open Targets id, paginate its credible
sets, keep only T2D/glycemic GWAS sets, pick the best-powered fine-mapped set, and
emit analysis-ready tables of candidate causal variants (with PIPs), L2G predictions,
and colocalisations.
"""
from __future__ import annotations
import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pandas as pd
from polaris import apis, config, utils

KW = config.T2D_TRAIT_KEYWORDS
EXACT = ("type 2 diabetes", "type ii diabetes")


def trait_priority(trait):
    t = (trait or "").lower()
    if any(e in t for e in EXACT):
        return 3
    if "diabetes" in t:
        return 2
    return 1  # glycemic adjacent (glucose/insulin/hba1c)


def score_cs(cs):
    st = cs["study"]
    n = st.get("nSamples") or ((st.get("nCases") or 0) + (st.get("nControls") or 0)) or 1
    return (trait_priority(st.get("traitFromSource")) * 1000
            + utils.finemap_rank(cs.get("finemappingMethod")) * 100
            + utils.safe_log10(n))


def main():
    loci_rows, cand_rows, l2g_rows, coloc_rows = [], [], [], []
    for rsid, locus, eff, control, coding, note in config.CURATED_LOCI:
        sys.stdout.write(f"[locus] {locus:16s} {rsid:13s} ")
        sys.stdout.flush()
        try:
            vid = apis.ot_resolve_variant(rsid)
        except Exception as e:
            vid = None
            sys.stdout.write(f"resolve-error {str(e)[:40]}\n"); continue
        if not vid:
            sys.stdout.write("UNRESOLVED\n"); continue
        core = apis.ot_variant_core(vid)
        light = apis.ot_cs_light_all(vid)
        kept = [r for r in light if r.get("studyType") == "gwas"
                and any(k in (r["study"].get("traitFromSource") or "").lower() for k in KW)]
        chrom, pos, ref, alt = utils.parse_ot_vid(vid)

        chosen_light = max(kept, key=score_cs) if kept else None
        chosen = (apis.ot_cs_detail_verified(vid, chosen_light["_idx"], chosen_light["studyLocusId"])
                  if chosen_light else None)
        n_members = 0
        best_l2g_gene = best_l2g_score = None
        best_coloc_h4 = best_coloc_gene = None

        if chosen:
            members = chosen["locus"]["rows"]
            n_members = len(members)
            # L2G
            l2g = sorted(chosen["l2GPredictions"]["rows"], key=lambda x: -x["score"]) if chosen["l2GPredictions"]["rows"] else []
            if l2g:
                best_l2g_gene = l2g[0]["target"]["approvedSymbol"]; best_l2g_score = l2g[0]["score"]
                for p in l2g:
                    l2g_rows.append(dict(locus=locus, rsid=rsid, study_id=chosen["study"]["id"],
                                         gene=p["target"]["approvedSymbol"], l2g_score=p["score"]))
            # coloc
            for cr in chosen["colocalisation"]["rows"]:
                g = (cr.get("otherStudyLocus") or {}).get("study", {}) or {}
                coloc_rows.append(dict(locus=locus, rsid=rsid, study_id=chosen["study"]["id"],
                                       method=cr.get("colocalisationMethod"), qtl_type=cr.get("rightStudyType"),
                                       h3=cr.get("h3"), h4=cr.get("h4"), clpp=cr.get("clpp"),
                                       qtl_trait=g.get("traitFromSource"), qtl_project=g.get("projectId")))
                if cr.get("h4") is not None and (best_coloc_h4 is None or cr["h4"] > best_coloc_h4):
                    best_coloc_h4 = cr["h4"]; best_coloc_gene = g.get("traitFromSource")
            # candidate causal variants (credible-set members)
            for m in members:
                mv = m["variant"]
                mc, mp, mr, ma = utils.parse_ot_vid(mv["id"])
                cand_rows.append(dict(
                    locus=locus, lead_rsid=rsid, rsid=(mv.get("rsIds") or [None])[0], vid=mv["id"],
                    chrom=mc, pos=mp, ref=mr, alt=ma, pip=m.get("posteriorProbability"),
                    is95=m.get("is95CredibleSet"), logBF=m.get("logBF"), beta=m.get("beta"),
                    is_lead=(mv["id"] == vid), finemap=chosen.get("finemappingMethod"),
                    study_id=chosen["study"]["id"], study_trait=chosen["study"].get("traitFromSource"),
                    n_samples=chosen["study"].get("nSamples"),
                    known_effector=eff, is_control=control, coding_lead=coding,
                    best_l2g_gene=best_l2g_gene, best_l2g_score=best_l2g_score,
                    best_coloc_h4=best_coloc_h4, best_coloc_gene=best_coloc_gene))
        else:
            # fallback: lead variant as singleton candidate
            cand_rows.append(dict(
                locus=locus, lead_rsid=rsid, rsid=rsid, vid=vid, chrom=chrom, pos=pos, ref=ref, alt=alt,
                pip=None, is95=None, logBF=None, beta=None, is_lead=True, finemap="lead-only",
                study_id=None, study_trait=None, n_samples=None,
                known_effector=eff, is_control=control, coding_lead=coding,
                best_l2g_gene=None, best_l2g_score=None, best_coloc_h4=None, best_coloc_gene=None))

        loci_rows.append(dict(
            locus=locus, lead_rsid=rsid, vid_lead=vid, chrom=chrom, pos_lead=pos,
            most_severe=(core.get("mostSevereConsequence") or {}).get("label"),
            known_effector=eff, is_control=control, coding_lead=coding, note=note,
            n_t2d_credible_sets=len(kept),
            chosen_study=(chosen["study"]["id"] if chosen else None),
            chosen_trait=(chosen["study"].get("traitFromSource") if chosen else None),
            chosen_finemap=(chosen.get("finemappingMethod") if chosen else "lead-only"),
            n_members=n_members, best_l2g_gene=best_l2g_gene, best_l2g_score=best_l2g_score,
            best_coloc_h4=best_coloc_h4, best_coloc_gene=best_coloc_gene))
        sys.stdout.write(f"vid={vid} T2D-cs={len(kept)} members={n_members} "
                         f"L2G={best_l2g_gene}({best_l2g_score if best_l2g_score is None else round(best_l2g_score,2)}) "
                         f"coloc_h4={None if best_coloc_h4 is None else round(best_coloc_h4,2)}\n")

    loci = pd.DataFrame(loci_rows); cand = pd.DataFrame(cand_rows)
    l2g = pd.DataFrame(l2g_rows); coloc = pd.DataFrame(coloc_rows)
    loci.to_csv(config.PROC / "loci.csv", index=False)
    cand.to_csv(config.PROC / "candidates.csv", index=False)
    l2g.to_csv(config.PROC / "l2g.csv", index=False)
    coloc.to_csv(config.PROC / "coloc_ot.csv", index=False)

    print("\n=== SUMMARY ===")
    print(f"loci processed       : {len(loci)}")
    print(f"loci with T2D cred set: {(loci['n_members']>0).sum()}")
    print(f"candidate variants    : {len(cand)}")
    print(f"  in 95% credible set : {cand['is95'].fillna(False).sum()}")
    print(f"L2G predictions       : {len(l2g)}")
    print(f"colocalisations       : {len(coloc)}")
    if len(cand):
        top = cand.dropna(subset=['pip']).sort_values('pip', ascending=False).head(8)
        print("\nTop candidate variants by PIP:")
        for _, r in top.iterrows():
            print(f"  {r['locus']:14s} {str(r['rsid']):13s} PIP={r['pip']:.3f} "
                  f"L2G={r['best_l2g_gene']} known={r['known_effector']}")


if __name__ == "__main__":
    main()
