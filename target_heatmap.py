"""SUPERSEDED - not used by any delivered figure.

An earlier version of the KZFP expression bar plus target-gene heatmap used in Fig 2C, 3C, 4C and 5C. Superseded by the panel functions inside the per-figure builders; this file still sets type at 3.2 and 5.4 pt.
No builder in this repository imports it. Kept for the record; running it
would produce output below the 8 pt legibility floor.
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram

# published red / white / purple divergent scale (pheatmap-like)
CMAP = LinearSegmentedColormap.from_list(
    "kzfp", ["#0000CD", "#8A7DD1", "#C9C3E3", "#E8E8E8", "#F0C0B4", "#E06B52", "#CC1F0A"])
BAR = "#4D4D4D"


def pick_cell_types(expr, gene, key, n=3):
    """Top-n and bottom-n cell types / tissues by nTPM, descending."""
    r = expr[expr["Gene name"] == gene][[key, "nTPM"]].copy()
    r = r.sort_values(["nTPM", key], ascending=[False, True]).reset_index(drop=True)
    hi = r.head(n)
    lo = r.sort_values(["nTPM", key], ascending=[True, True]).head(n)
    lo = lo.sort_values(["nTPM", key], ascending=[False, True])
    sel = pd.concat([hi, lo]).drop_duplicates(key)
    return sel[key].tolist(), dict(zip(sel[key], sel.nTPM))


def target_matrix(expr, targets, cols, key):
    m = expr[expr["Gene name"].isin(targets) & expr[key].isin(cols)]
    W = m.pivot_table(index="Gene name", columns=key, values="nTPM", aggfunc="mean")
    W = W.reindex(columns=cols)
    W = W.dropna(how="all")
    return np.log(W + 1)


def heatmap_panel(fig, gs, expr, kzfp, targets, key, bar_label, hm_label, box_label,
                  vmax=None, show_row_labels=True, row_fontsize=3.2):
    cols, barvals = pick_cell_types(expr, kzfp, key)
    L = target_matrix(expr, targets, cols, key)
    C = L.sub(L.mean(axis=1), axis=0)                    # row-centred log(nTPM+1)
    Z = linkage(C.values, method="complete", metric="euclidean")
    dn = dendrogram(Z, no_plot=True)
    order = dn["leaves"]
    C = C.iloc[order]; L = L.iloc[order]
    if vmax is None:
        vmax = float(np.ceil(np.nanmax(np.abs(C.values))))

    axbar = fig.add_subplot(gs[0, 1])
    axden = fig.add_subplot(gs[1, 0])
    axhm = fig.add_subplot(gs[1, 1])
    axbox = fig.add_subplot(gs[2, 1])

    axbar.bar(range(len(cols)), [barvals[c] for c in cols], width=0.70, color=BAR, edgecolor="none")
    axbar.set_xlim(-0.5, len(cols) - 0.5); axbar.set_xticks([])
    axbar.set_ylabel(bar_label, fontsize=5.4, linespacing=1.15)
    axbar.tick_params(labelsize=5.4)
    for sp in ("top", "right"):
        axbar.spines[sp].set_visible(False)

    with plt.rc_context({"lines.linewidth": 0.4}):
        dendrogram(Z, ax=axden, orientation="left", no_labels=True,
                   link_color_func=lambda *a: "black")
    axden.invert_yaxis(); axden.axis("off")

    axhm.imshow(C.values, aspect="auto", cmap=CMAP, vmin=-vmax, vmax=vmax, interpolation="nearest")
    axhm.set_xticks([]); axhm.set_yticks([])
    if show_row_labels:
        axhm.set_yticks(range(len(C)))
        axhm.set_yticklabels(C.index, fontsize=row_fontsize, style="italic")
        axhm.yaxis.tick_right(); axhm.tick_params(axis="y", length=0, pad=1)
    for sp in axhm.spines.values():
        sp.set_visible(False)

    cax = fig.add_axes([axhm.get_position().x0 - 0.155, axhm.get_position().y0 + 0.02,
                        0.012, axhm.get_position().height * 0.34])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax), cmap=CMAP),
                      cax=cax, ticks=[-vmax, 0, vmax])
    cb.ax.tick_params(labelsize=5.0, length=1.5, pad=1); cb.outline.set_linewidth(0.4)
    cax.set_ylabel(hm_label, fontsize=5.4, labelpad=1, linespacing=1.15)
    cax.yaxis.set_label_position("left")

    data = [L[c].dropna().values for c in cols]
    bp = axbox.boxplot(data, positions=range(len(cols)), widths=0.62, patch_artist=True,
                       showfliers=False, whis=0, showcaps=False,
                       medianprops=dict(color="black", lw=0.6),
                       boxprops=dict(facecolor="#D9D9D9", edgecolor="#6E6E6E", lw=0.4),
                       whiskerprops=dict(lw=0))
    for i, v in enumerate(data):
        axbox.plot([i], [np.mean(v)], marker="_", color="red", ms=4.0, mew=0.8)
        axbox.plot([i], [np.median(v)], marker="+", color="black", ms=2.6, mew=0.6)
    axbox.set_xlim(-0.5, len(cols) - 0.5)
    axbox.set_xticks(range(len(cols)))
    axbox.set_xticklabels(cols, rotation=45, ha="right", fontsize=5.4)
    axbox.set_ylabel(box_label, fontsize=5.4, linespacing=1.15)
    axbox.tick_params(labelsize=5.4)
    for sp in ("top", "right"):
        axbox.spines[sp].set_visible(False)
    return dict(kzfp=kzfp, cell_types=cols, kzfp_nTPM={c: round(float(barvals[c]), 2) for c in cols},
                n_targets_requested=len(targets), n_targets_plotted=int(len(C)),
                missing=sorted(set(targets) - set(C.index)),
                box_median={c: round(float(np.median(L[c].dropna())), 3) for c in cols},
                box_mean={c: round(float(np.mean(L[c].dropna())), 3) for c in cols})
