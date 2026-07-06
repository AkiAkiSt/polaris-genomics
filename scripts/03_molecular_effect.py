"""POLARIS Step 3: transparent multimodal molecular-effect feature matrix.

Runs the white-box molecular-effect orchestrator over every candidate causal variant
and writes data/processed/molecular_features.csv plus motif_records.json (for logos).
"""
from __future__ import annotations
import sys, json, time, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pandas as pd, numpy as np
from polaris import apis, config, motifs as M

def main(test=False):
    cand = pd.read_csv(config.PROC / "candidates.csv")
    if test:
        cand = cand[cand.is_lead].head(4)
    pwms = M.load_pwms(str(config.EXT / "JASPAR2024_CORE_vertebrates_nr.txt"))
    print(f"loaded {len(pwms)} PWMs; annotating {len(cand)} variants")
    rows, motif_json = [], {}
    t0 = time.time()
    for i, (_, r) in enumerate(cand.iterrows()):
        from polaris.molecular import annotate_variant
        try:
            f = annotate_variant(str(r["chrom"]), int(r["pos"]), str(r["ref"]),
                                 str(r["alt"]), r["vid"], pwms)
        except Exception as e:
            f = dict(chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"], vid=r["vid"],
                     _error=str(e)[:80], _motif_records=[])
        recs = f.pop("_motif_records", [])
        motif_json[r["vid"]] = recs
        f["locus"] = r["locus"]; f["rsid"] = r["rsid"]; f["pip"] = r["pip"]
        rows.append(f)
        if (i + 1) % 25 == 0 or i == len(cand) - 1:
            el = time.time() - t0
            sys.stdout.write(f"  [{i+1}/{len(cand)}] {el:.0f}s "
                             f"({el/(i+1):.2f}s/var) last={r['locus']}\n"); sys.stdout.flush()
    df = pd.DataFrame(rows)
    if not test:
        df.to_csv(config.PROC / "molecular_features.csv", index=False)
        json.dump(motif_json, open(config.PROC / "motif_records.json", "w"))
        print(f"\nsaved molecular_features.csv ({df.shape}) + motif_records.json")
        # quick sanity
        ok = df.get("motif_ref_ok")
        if ok is not None:
            print(f"ref-base concordance: {ok.fillna(False).mean():.1%}")
        print("mean motif_disruption:", round(float(np.nanmean(df['motif_motif_disruption'])), 2))
        print("variants in a cCRE:", int(df['in_ccre'].fillna(0).sum()))
        print("median phyloP100:", round(float(np.nanmedian(df['phylop'])), 3))
    else:
        cols = [c for c in df.columns if not c.startswith("_")]
        print(df[["locus", "rsid", "motif_top_disrupt_tf", "motif_top_disrupt_dLLR",
                  "phylop", "phastcons", "ccre_class", "gnomad_af"]].to_string())
    return df


if __name__ == "__main__":
    main(test=("--test" in sys.argv))
