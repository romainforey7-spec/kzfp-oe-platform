"""SUPERSEDED - not used by any delivered figure.

An earlier version of the AP-MS intensity figures. Superseded by build_fig6_fig.py; this file still sets type at 6.5 pt.
No builder in this repository imports it. Kept for the record; running it
would produce output below the 8 pt legibility floor.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

C_GFP, C_FL, C_DS = "#9E9E9E", "#2C6E9B", "#E08214"
C_SCAN_DOM = "#1B5E3F"
COND_LABEL = {"GFP": "GFP control", "FL": "full length", "dSCAN": "\u0394SCAN"}


def save(fig, stem, outdir="out"):
    fig.savefig(f"{outdir}/{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{outdir}/{stem}.png", dpi=300, bbox_inches="tight")
    return f"{outdir}/{stem}.svg"


def check_overlaps(fig):
    r = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    return [(a.get_text(), b.get_text())
            for i, (a, ba) in enumerate(tx) for b, bb in tx[i + 1:] if ba.overlaps(bb)]


def condition_bars(ax, genes, bait, condmean, obsmat, dsub, cond_cols, scan_set,
                   show_points=True, fontsize=6.5):
    """Horizontal grouped bars: GFP / full length / dSCAN mean normalised log2 intensity.

    Individual samples are overlaid as points so the reader can see the n.
    """
    conds = ["GFP", f"{bait}_FL", f"{bait}_dSCAN"]
    cols = [C_GFP, C_FL, C_DS]
    y = np.arange(len(genes))
    h = 0.26
    offs = [h, 0, -h]
    for cd, cl, off in zip(conds, cols, offs):
        vals = [condmean.loc[g, cd] for g in genes]
        ax.barh(y + off, vals, h, color=cl, zorder=3,
                edgecolor="white", linewidth=0.4)
        if show_points:
            scols = [c for c in dsub.column if cond_cols[c] == cd]
            for gi, g in enumerate(genes):
                pts = [obsmat.loc[g, c] for c in scols]
                obs = [p for p in pts if np.isfinite(p)]
                if obs:
                    ax.scatter(obs, [y[gi] + off] * len(obs), s=3.2, color="#222222",
                               zorder=5, linewidths=0)
    labs = [f"{g}*" if g in scan_set else g for g in genes]
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=fontsize)
    for i, g in enumerate(genes):
        if g in scan_set:
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.set_xlabel("normalised log$_2$ intensity", fontsize=6.5)
    ax.margins(y=0.03)
    return y
