"""Figure 1 - POLARIS conceptual architecture & factorized posterior."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from polaris import viz, config
C = viz.C


def box(ax, x, y, w, h, title, sub, fc, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=fc, ec="none", zorder=2))
    ax.text(x + w / 2, y + h * 0.66, title, ha="center", va="center", fontsize=10,
            fontweight="bold", color=tc, zorder=3)
    ax.text(x + w / 2, y + h * 0.28, sub, ha="center", va="center", fontsize=7.6, color=tc, zorder=3)


def arrow(ax, x1, y1, x2, y2, col=None):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16,
                                 lw=2, color=col or C["slate"], zorder=1))


fig, ax = plt.subplots(figsize=(14.5, 8.4))
ax.set_xlim(0, 14.5); ax.set_ylim(0, 8.4); ax.axis("off")

ax.text(7.25, 8.05, "POLARIS", ha="center", fontsize=20, fontweight="bold", color=C["navy"])
ax.text(7.25, 7.6, "Posterior-calibrated triangulation of non-coding variant evidence -> mechanism & target gene",
        ha="center", fontsize=10.5, color=C["slate"])

# pipeline flow
y = 5.7; h = 1.15; w = 2.5
box(ax, 0.2, y, w, h, "1  GWAS loci", "T2D credible sets\n(Open Targets + SuSiE)", C["navy"])
box(ax, 3.05, y, w, h, "2  Molecular effect", "TF-motif dLLR, conservation,\ncCRE  (transparent)", C["teal"])
box(ax, 5.9, y, w, h, "3  Gene linkage", "L2G, colocalization,\nnearest / distance", C["violet"])
box(ax, 8.75, y, w, h, "4  Integrate", "reliability-weighted\nlatent-truth EM", C["amber"], tc=C["ink"])
box(ax, 11.6, y, 2.7, h, "5  Calibrated output", "ranked variant->gene\n+ mechanism + CI", C["crimson"])
for x in (2.7, 5.55, 8.4, 11.25):
    arrow(ax, x, y + h / 2, x + 0.35, y + h / 2)

# triangulation concept (three converging arrows into posterior)
ax.text(7.25, 4.75, "Triangulation: credibility = convergence across independent channels",
        ha="center", fontsize=10, style="italic", color=C["ink"])

# factorized posterior box
fb = FancyBboxPatch((1.6, 2.55), 11.3, 1.5, boxstyle="round,pad=0.02,rounding_size=0.03",
                    fc=C["mist"], ec=C["navy"], lw=1.5)
ax.add_patch(fb)
ax.text(7.25, 3.65, "Factorized posterior  (functional  $\\neq$  pathogenic)", ha="center",
        fontsize=11.5, fontweight="bold", color=C["navy"])
ax.text(7.25, 3.05,
        r"$P(\mathrm{causal\ via}\ g\,|\,E)\;=\;P(F\,|\,E_{\mathrm{molecular}})\;\times\;\mathrm{PIP}\;\times\;P(D_g\,|\,E_{\mathrm{linkage}})$",
        ha="center", fontsize=13.5, color=C["ink"])
for xx, lab in ((4.35, "is it functional?"), (7.05, "is it causal?"), (9.95, "which gene?")):
    ax.text(xx, 2.68, lab, ha="center", fontsize=7.8, style="italic", color=C["slate"])

# guardrails
ax.add_patch(FancyBboxPatch((0.4, 0.5), 6.5, 1.6, boxstyle="round,pad=0.02,rounding_size=0.03",
                            fc="white", ec=C["teal"], lw=1.4))
ax.text(3.65, 1.85, "Guardrail 1 - functional $\\neq$ pathogenic", ha="center", fontsize=9.5,
        fontweight="bold", color=C["teal"])
ax.text(3.65, 1.15, "molecular function and disease causation are\nseparate factors; rarity is INVERTED for\ncommon-variant complex traits",
        ha="center", fontsize=8.2, color=C["ink"])

ax.add_patch(FancyBboxPatch((7.6, 0.5), 6.5, 1.6, boxstyle="round,pad=0.02,rounding_size=0.03",
                            fc="white", ec=C["amber"], lw=1.4))
ax.text(10.85, 1.85, "Guardrail 2 - distrust any single oracle", ha="center", fontsize=9.5,
        fontweight="bold", color=C["amber"])
ax.text(10.85, 1.15, "each channel is a noisy annotator; the EM learns\nits reliability (phyloP > motif > ...); channels are\nlargely independent (Marchenko-Pastur)",
        ha="center", fontsize=8.2, color=C["ink"])

viz.save(fig, "fig_pipeline")
print("pipeline schematic saved")
