"""Evidence-landscape heatmap: top POLARIS hypotheses x multimodal evidence channels,
annotated with nominated gene, known effector, and transparent mechanism."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from polaris import viz, config
C = viz.C

df = pd.read_csv(config.PROC / "integrated.csv")
df = df[df.coding_lead != True]                       # focus on non-coding
top = df.sort_values("POLARIS", ascending=False).drop_duplicates("locus").head(22).copy()

chans = [("pip", "PIP\n(fine-map)"), ("P_functional", "P(func)\n(molecular)"),
         ("phylop_max", "phyloP\n(cons.)"), ("phastcons470", "phastCons\n470"),
         ("motif_motif_disruption", "motif\ndisrupt"), ("ccre_enhancer", "enhancer\ncCRE"),
         ("best_l2g_score", "L2G\n(linkage)")]
Mtx = np.zeros((len(top), len(chans)))
for j, (c, _) in enumerate(chans):
    v = top[c].astype(float).values
    v = np.nan_to_num(v, nan=np.nanmedian(v))
    lo, hi = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    Mtx[:, j] = np.clip((v - lo) / (hi - lo + 1e-9), 0, 1)

fig, ax = plt.subplots(figsize=(13.5, 9.6))
im = ax.imshow(Mtx, aspect="auto", cmap=viz.cmap("sequential"), norm=Normalize(0, 1))
ax.set_xticks(range(len(chans))); ax.set_xticklabels([c[1] for c in chans], fontsize=9)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f"{r.locus}  {r.rsid}" for r in top.itertuples()], fontsize=8.5)
ax.xaxis.tick_top()
# group separators
for x in (0.5, 5.5):
    ax.axvline(x, color="white", lw=3)
ax.text(0, -1.3, "STATISTICAL", fontsize=8, color=C["navy"], fontweight="bold")
ax.text(1.4, -1.3, "MOLECULAR (transparent)", fontsize=8, color=C["teal"], fontweight="bold")
ax.text(6, -1.3, "LINKAGE", fontsize=8, color=C["violet"], fontweight="bold")

# right-side annotations: gene call (hit/miss) + mechanism TF
for i, r in enumerate(top.itertuples()):
    hit = (str(r.polaris_gene) == str(r.known_effector))
    mark = "OK" if hit else "x"
    col = C["grass"] if hit else C["crimson"]
    ax.text(len(chans) - 0.35, i, f"-> {r.polaris_gene} [{mark}]", fontsize=7.8, va="center", color=col)
    tf = r.motif_top_disrupt_tf
    if isinstance(tf, str):
        ax.text(len(chans) + 1.15, i, tf, fontsize=7.3, va="center", color=C["slate"], style="italic")
ax.text(len(chans) - 0.35, -1.0, "gene call", fontsize=7.5, color=C["ink"])
ax.text(len(chans) + 1.15, -1.0, "disrupted TF", fontsize=7.5, color=C["ink"])
ax.set_xlim(-0.5, len(chans) + 2.2)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.13)
cb.set_label("normalized evidence (per channel)", fontsize=9)
ax.set_title("Evidence landscape of top POLARIS non-coding causal hypotheses in type 2 diabetes",
             fontsize=13, fontweight="bold", pad=38)
viz.save(fig, "fig_evidence_heatmap")
print(f"heatmap saved: {len(top)} loci. gene-call hits: {(top.polaris_gene==top.known_effector).sum()}/{len(top)}")
