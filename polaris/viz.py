"""Publication-grade plotting identity for POLARIS (matplotlib, Agg).
Cohesive palette, clean axes, TextPath sequence logos, consistent saving."""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
import matplotlib.font_manager as fm
import numpy as np
from . import config

C = config.VIS
plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.titlesize": 13, "axes.titleweight": "bold", "axes.labelsize": 11.5,
    "axes.edgecolor": C["ink"], "axes.linewidth": 1.0, "axes.labelcolor": C["ink"],
    "xtick.color": C["ink"], "ytick.color": C["ink"], "text.color": C["ink"],
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e6edf4", "grid.linewidth": 0.8,
    "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
})

SEQUENTIAL = ["#0b1f3a", "#163a5f", "#1f6f8b", "#1f9e89", "#7ecb9c", "#f0e6a8"]
DIVERGING = ["#b3123a", "#e08a9b", "#f3e6ea", "#cfe3ee", "#3b7fb0", "#0b3d63"]
BASE_COLORS = {"A": "#3a8d5b", "C": "#1f6f8b", "G": "#f0a202", "T": "#b3123a"}


def cmap(name="sequential"):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, SEQUENTIAL if name == "sequential" else DIVERGING)


def new(figsize=(7, 4.3), n=1, m=1, **kw):
    fig, ax = plt.subplots(n, m, figsize=figsize, **kw)
    return fig, ax


def save(fig, name, also_svg=True):
    p = config.FIG / f"{name}.png"
    fig.savefig(p, facecolor="white")
    if also_svg:
        fig.savefig(config.FIG / f"{name}.svg", facecolor="white")
    plt.close(fig)
    return p


def letter(ax, base, x, y, w, h):
    """Draw a single scaled sequence-logo letter."""
    tp = TextPath((0, 0), base, size=1, prop=fm.FontProperties(family="DejaVu Sans", weight="bold"))
    bb = tp.get_extents()
    sx = w / (bb.width or 1); sy = h / (bb.height or 1)
    t = Affine2D().translate(-bb.x0, -bb.y0).scale(sx, sy).translate(x, y)
    ax.add_patch(PathPatch(tp.transformed(t), fc=BASE_COLORS.get(base, "#333"), ec="none"))


def sequence_logo(ax, prob, ic=None, start=0, highlight=None, title=None):
    """prob: 4xw probabilities (ACGT). ic: per-position bits (else computed uniform)."""
    bases = "ACGT"
    w = prob.shape[1]
    if ic is None:
        bg = np.array([.295, .205, .205, .295])
        ic = np.sum(prob * np.log2((prob + 1e-9) / bg[:, None]), 0)
    for i in range(w):
        order = np.argsort(prob[:, i])
        y = 0.0
        for j in order:
            h = prob[j, i] * max(ic[i], 0)
            if h > 1e-3:
                letter(ax, bases[j], start + i + 0.05, y, 0.9, h)
            y += h
    if highlight is not None:
        ax.axvspan(start + highlight + 0.0, start + highlight + 1.0, color=C["amber"], alpha=0.18, zorder=0)
    ax.set_xlim(start, start + w)
    ax.set_ylim(0, 2.05)
    ax.set_ylabel("bits")
    ax.grid(False)
    if title:
        ax.set_title(title)


def despine_all(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
