"""Fetch protein-coding genes around each locus and compute variant->TSS distances.
Feeds the POLARIS target-gene linkage model (nearest-gene + distance channels)."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pandas as pd
from polaris import apis, config, utils

WIN = 600_000


def main():
    loci = pd.read_csv(config.PROC / "loci.csv")
    rows = []
    for _, r in loci.iterrows():
        chrom, pos = str(r["chrom"]), int(r["pos_lead"])
        try:
            genes = apis.ens_overlap_genes(chrom, pos - WIN, pos + WIN)
        except Exception as e:
            sys.stdout.write(f"{r['locus']}: gene fetch error {str(e)[:40]}\n"); continue
        gl = []
        for g in genes:
            if g.get("biotype") != "protein_coding":
                continue
            strand = g.get("strand", 1)
            tss = g["start"] if strand == 1 else g["end"]
            gl.append(dict(locus=r["locus"], lead_pos=pos,
                           gene=g.get("external_name") or g.get("id"), gene_id=g.get("id"),
                           start=g["start"], end=g["end"], strand=strand, tss=tss,
                           dist_to_tss=abs(pos - tss),
                           dist_to_body=0 if g["start"] <= pos <= g["end"] else min(abs(pos - g["start"]), abs(pos - g["end"]))))
        gl.sort(key=lambda d: d["dist_to_tss"])
        for rank, d in enumerate(gl):
            d["tss_rank"] = rank
            d["is_nearest"] = int(rank == 0)
            rows.append(d)
        sys.stdout.write(f"{r['locus']:16s} genes={len(gl)} nearest={gl[0]['gene'] if gl else None}\n")
    df = pd.DataFrame(rows)
    df.to_csv(config.PROC / "genes.csv", index=False)
    print(f"\nsaved genes.csv {df.shape}; loci={df.locus.nunique()}")


if __name__ == "__main__":
    main()
