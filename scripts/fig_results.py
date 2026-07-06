"""Main results figure: reliability, key findings, random-matrix oracle analysis,
evidence manifold, calibration, gene nomination."""
from __future__ import annotations
import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from polaris import viz, geometry as GE, integrate as IN, config
C = viz.C

df = pd.read_csv(config.PROC / "integrated.csv")
M = json.load(open(config.TAB / "polaris_metrics.json"))
rel = pd.read_csv(config.TAB / "channel_reliability.csv")

fig = plt.figure(figsize=(15.5, 9.6))
gs = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.30)

# ---- A: EM-inferred channel reliability ----
axA = fig.add_subplot(gs[0, 0])
r = rel.sort_values("weight")
labels = [s.replace("motif_motif_disruption", "motif disruption").replace("phylop_max", "phyloP (cons.)")
          .replace("phastcons470", "phastCons470").replace("ccre_enhancer", "enhancer cCRE")
          .replace("splice_delta", "splice delta") for s in r.channel]
bars = axA.barh(labels, r.weight, color=[C["teal"] if w == r.weight.max() else C["slate"] for w in r.weight])
axA.set_xlabel("inferred reliability weight"); axA.set_title("A  Which oracle to trust (latent-truth EM)")
for b, d in zip(bars, r.dprime):
    axA.text(b.get_width() + 0.01, b.get_y() + b.get_height() / 2, f"d'={d:.2f}", va="center", fontsize=8.5)
axA.set_xlim(0, max(r.weight) * 1.25)

# ---- B: the two key findings (conservation enrichment + rarity inversion) ----
axB = fig.add_subplot(gs[0, 1])
pos, neg = df[df.pip > 0.5], df[df.pip < 0.05]
data = [pos["phylop_max"].dropna(), neg["phylop_max"].dropna(),
        pos["rarity"].dropna(), neg["rarity"].dropna()]
bp = axB.boxplot(data, positions=[1, 2, 4, 5], widths=0.7, patch_artist=True, showfliers=False)
for i, patch in enumerate(bp["boxes"]):
    patch.set_facecolor([C["teal"], C["mist"], C["amber"], C["mist"]][i]); patch.set_alpha(0.8)
for med in bp["medians"]:
    med.set_color(C["ink"])
axB.set_xticks([1.5, 4.5]); axB.set_xticklabels(["conservation\n(phyloP)", "rarity\n(-log10 AF)"])
axB.set_ylabel("value")
axB.set_title("B  Causal-variant signatures")
ce, ri = M["finding_conservation_enrichment"], M["finding_rarity_inversion"]
axB.text(1.5, axB.get_ylim()[1] * 0.93, f"AUC {ce['auc']}\np={ce['p']:.0e}", ha="center", fontsize=8, color=C["teal"])
axB.text(4.5, axB.get_ylim()[1] * 0.93, f"AUC {ri['auc']}\nINVERTED\np={ri['p']:.0e}", ha="center", fontsize=8, color=C["crimson"])
axB.text(0.5, -0.22, "causal common-variant signals are MORE common, not rarer (functional != pathogenic)",
         transform=axB.transAxes, ha="center", fontsize=7.5, style="italic", color=C["slate"])
axB.legend([bp["boxes"][0], bp["boxes"][1]], ["fine-mapped causal", "LD passenger"], fontsize=8, loc="upper right")

# ---- C: Marchenko-Pastur oracle spectrum ----
axC = fig.add_subplot(gs[0, 2])
ev = np.array(M["oracle_spectrum"]["evals"]); gamma = M["oracle_spectrum"]["gamma"]
mp_hi = M["oracle_spectrum"]["mp_hi"]
xs = np.linspace(0.01, max(ev) * 1.1, 400)
axC.hist(ev, bins=12, density=True, color=C["mist"], edgecolor=C["slate"], alpha=0.9, label="oracle eigenvalues")
axC.plot(xs, GE.mp_pdf(xs, gamma, 1.0), color=C["navy"], lw=2, label="Marchenko-Pastur null")
axC.axvline(mp_hi, color=C["crimson"], ls="--", lw=1.5, label=f"MP edge {mp_hi:.2f}")
for e in ev[ev > mp_hi]:
    axC.axvline(e, color=C["amber"], lw=2.2)
axC.set_xlabel("eigenvalue"); axC.set_ylabel("density")
axC.set_title("C  Oracles are largely independent")
axC.legend(fontsize=8)
axC.text(0.97, 0.62, f"{M['oracle_spectrum']['n_signal_dims']} of {M['oracle_spectrum']['n_channels']}\nchannels\nshared",
         transform=axC.transAxes, ha="right", fontsize=9, color=C["amber"])

# ---- D: PCA of the evidence space (interpretable axes) ----
axD = fig.add_subplot(gs[1, 0])
gcols = ["P_functional", "pip", "phylop_max", "motif_motif_disruption", "phastcons470", "best_l2g_score", "rarity"]
clean = df.dropna(subset=gcols).copy()
Xz = (clean[gcols].values - clean[gcols].mean().values) / (clean[gcols].std().values + 1e-9)
U, S, Vt = np.linalg.svd(Xz - Xz.mean(0), full_matrices=False)
pcs = (Xz - Xz.mean(0)) @ Vt[:2].T
var = (S ** 2 / (S ** 2).sum())[:2]
sc = axD.scatter(pcs[:, 0], pcs[:, 1], c=clean["POLARIS"], cmap=viz.cmap("sequential"),
                 s=22, alpha=0.85, edgecolor="none")
ctrlmask = (clean["is_control"] == True).values
axD.scatter(pcs[ctrlmask, 0], pcs[ctrlmask, 1], s=70, facecolor="none",
            edgecolor=C["crimson"], lw=1.4, label="control loci")
plt.colorbar(sc, ax=axD, label="POLARIS score", fraction=0.046)
axD.set_xlabel(f"PC1 ({var[0]*100:.0f}% var)"); axD.set_ylabel(f"PC2 ({var[1]*100:.0f}% var)")
axD.set_title("D  Evidence-space geometry (PCA)")
axD.legend(fontsize=8, loc="best"); axD.grid(False)

# ---- E: calibration reliability curve ----
axE = fig.add_subplot(gs[1, 1])
try:
    co = np.load(config.PROC / "calib_oof.npy")
    xs2, ys2, ns = IN.reliability_curve(co[:, 0], co[:, 1], bins=6)
    axE.plot([0, 1], [0, 1], ls="--", color=C["slate"], lw=1)
    axE.plot(xs2, ys2, "o-", color=C["teal"], ms=8, lw=2)
    for x, y, n in zip(xs2, ys2, ns):
        axE.text(x, y + 0.04, str(n), fontsize=7, ha="center", color=C["slate"])
    axE.set_xlim(0, 1); axE.set_ylim(0, 1)
except Exception:
    pass
axE.set_xlabel("predicted P(causal)"); axE.set_ylabel("observed frequency")
axE.set_title("E  Calibration (locus-aware CV)")
axE.text(0.05, 0.9, f"ECE = {M.get('calib_ece')}\nBrier = {M.get('calib_brier')}",
         transform=axE.transAxes, fontsize=9.5, color=C["teal"], va="top")

# ---- F: gene nomination top-1 ----
axF = fig.add_subplot(gs[1, 2])
methods = [("POLARIS\n(L2G linkage)", M["gene_top1_POLARIS(=L2G)"], C["teal"]),
           ("L2G+distance\n(naive)", M["gene_top1_naive_L2G+distance_combo"], C["slate"]),
           ("nearest\ngene", M["gene_top1_nearest"], C["slate"]),
           ("distance\nto TSS", M["gene_top1_distance"], C["slate"])]
names = [m[0] for m in methods]; vals = [m[1] for m in methods]; cols = [m[2] for m in methods]
axF.bar(names, vals, color=cols, alpha=0.9)
axF.set_ylabel("target-gene top-1 accuracy"); axF.set_ylim(0, 1)
axF.set_title("F  Gene nomination (13 control loci)")
for i, v in enumerate(vals):
    axF.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
axF.text(0.5, -0.30, "the 2 misses are the known-hard cases: distal FTO->IRX3 (POLARIS flags it)\nand imprinted DUSP8",
         transform=axF.transAxes, ha="center", fontsize=7.5, style="italic", color=C["slate"])

fig.suptitle("POLARIS — calibrated triangulation of non-coding variant evidence in type 2 diabetes",
             fontsize=15.5, fontweight="bold", y=0.985)
viz.save(fig, "fig_results_main")
print(f"saved fig_results_main. PCA PC1/PC2 var = {var[0]:.2f}/{var[1]:.2f} | n clean variants: {len(clean)}")
