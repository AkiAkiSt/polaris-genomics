"""Generate a 1200x630 social-preview (OG) card for the POLARIS site."""
import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from polaris import config

m = json.load(open(config.TAB / "polaris_metrics.json"))
gene = m["gene_top1_POLARIS(=L2G)"]; ece = m["calib_ece"]
fig = plt.figure(figsize=(12, 6.3), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 12); ax.set_ylim(0, 6.3); ax.axis("off")
# gradient-ish dark background
ax.add_patch(plt.Rectangle((0, 0), 12, 6.3, color="#0b1f3a", zorder=0))
rng = np.random.default_rng(3)
xs, ys = rng.uniform(0, 12, 90), rng.uniform(0, 6.3, 90)
ax.scatter(xs, ys, s=rng.uniform(1, 9, 90), c="#cfe3f2", alpha=0.35, zorder=1)
# triangulation motif to a bright star
px, py = 9.6, 4.5
for _ in range(3):
    sx, sy = rng.uniform(1, 8, 1)[0], rng.uniform(1, 5.5, 1)[0]
    ax.plot([sx, px], [sy, py], color="#7fd8c8", lw=1, alpha=0.35, zorder=1)
ax.scatter([px], [py], s=220, c="#7fd8c8", alpha=0.25, zorder=2)
ax.scatter([px], [py], s=40, c="#eafaf6", zorder=3)
ax.text(0.8, 4.35, "POLARIS", fontsize=76, fontweight="bold", color="white",
        family="DejaVu Serif", zorder=4)
ax.text(0.85, 3.35, "Calibrated triangulation of non-coding variant evidence",
        fontsize=23, color="#c3d6ea", zorder=4)
ax.text(0.85, 2.78, "into mechanism and target gene", fontsize=23, color="#c3d6ea", zorder=4)
chips = [(f"{gene*100:.0f}%", "gene top-1"), ("829", "variants"),
         (f"{ece:.3f}", "calibration"), ("2", "diseases")]
for i, (v, l) in enumerate(chips):
    x = 0.85 + i * 2.75
    ax.add_patch(FancyBboxPatch((x, 0.7), 2.4, 1.35, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="#12263f", ec="#7fd8c8", lw=1, zorder=4))
    ax.text(x + 0.25, 1.45, v, fontsize=30, fontweight="bold", color="white", family="DejaVu Serif", zorder=5)
    ax.text(x + 0.27, 0.98, l, fontsize=13, color="#9fbdd6", zorder=5)
ax.text(0.85, 0.28, "transparent  ·  honest  ·  reproducible from public APIs",
        fontsize=14, color="#7fd8c8", zorder=4)
fig.savefig(config.ROOT / "website" / "assets" / "og.png", dpi=100, facecolor="#0b1f3a")
print("saved og.png")
