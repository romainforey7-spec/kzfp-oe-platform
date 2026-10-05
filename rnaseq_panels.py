"""Panel builders for the KZFP RNA-seq figures (volcano + alluvial), in the
original published visual style. Colours sampled from the published TIFFs."""
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from scipy.stats import fisher_exact

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "axes.titlesize": 6,
})

# --- published palettes (sampled from Figure 2.TIF) ---------------------------
VOLC = {("DOWN", True): "#00008B", ("DOWN", False): "#ADD8E6",
        ("UP",   True): "#8B0000", ("UP",   False): "#FFA07A",
        ("NS",   True): "#696969", ("NS",   False): "#D3D3D3"}
ALLUV = {("DOWN", "<1Kb"): "#4C4CAE", ("DOWN", "1 to 10Kb"): "#4C4CFF",
         ("DOWN", "10 to 50Kb"): "#C5E4ED",
         ("UP", "<1Kb"): "#AE4C4C", ("UP", "1 to 10Kb"): "#FF4C4C",
         ("UP", "10 to 50Kb"): "#FF9494"}
BINS = ["<1Kb", "1 to 10Kb", "10 to 50Kb"]
FS = {"star": 7.0, "leg": 6.0, "legs": 5.6, "node": 5.2}
NODE_GREY = "#E6E6E6"


def bubble_area(v):
    """ggplot-like area scale: legend keys 20/40/60 of the (score/10) value."""
    v = np.asarray(v, dtype=float)
    return (1.15 + 0.78 * np.sqrt(np.clip(v, 0, None))) ** 2


def stars(pv):
    return "ns" if pv > 0.05 else ("*" if pv > 0.01 else ("**" if pv > 0.001 else "***"))


def signed_log2(fc):
    fc = np.asarray(fc, dtype=float)
    return np.sign(fc) * np.log2(np.abs(fc))


def gene_table(merged, window):
    """One row per gene: DE status + nearest genuine peak (dist >= 0) + score."""
    g = merged.copy()
    g["dist"] = pd.to_numeric(g["dist"], errors="coerce")
    g.loc[g["dist"] < 0, "dist"] = np.nan          # bedtools "no feature" sentinel
    g["pscore"] = pd.to_numeric(g["pscore"], errors="coerce")
    g = g[g.ensembl.notna() & (g.ensembl != ".") & (g.ensembl != "ensembl")]
    g["foldChange"] = pd.to_numeric(g["foldChange"], errors="coerce")
    g["padj"] = pd.to_numeric(g["padj"], errors="coerce")
    g = g[g.foldChange.notna() & g.padj.notna()]
    g = g.sort_values("dist", na_position="last").drop_duplicates("ensembl")
    g["cat"] = np.where(g.Direction == "DOWN", "DOWN",
                        np.where(g.Direction == "UP", "UP", "NS"))
    g["haspeak"] = g["dist"].notna() & (g["dist"] <= window)
    g["l2fc"] = signed_log2(g.foldChange)
    g["nl2p"] = -np.log2(pd.to_numeric(g.padj, errors="coerce").clip(lower=1e-300))
    return g


def alluvial_table(merged, maxdist=50000):
    """One row per (DE gene x repeat overlapped by its nearest peak)."""
    g = merged.copy()
    g["dist"] = pd.to_numeric(g["dist"], errors="coerce")
    s = g[g.Direction.isin(["UP", "DOWN"]) & (g["dist"] >= 0) & (g["dist"] <= maxdist)].copy()
    s["bin"] = pd.cut(s["dist"], [-1, 1000, 10000, maxdist], labels=BINS).astype(str)
    s["assgn"] = s.rclass.fillna(".").replace({".": "Not RE", "nan": "Not RE"})
    return s


def group_assign(s, min_frac=0.045, keep=None):
    """Collapse rare repeat classes into 'Diverse REs', keeping 'Not RE' separate."""
    n = len(s)
    vc = s["assgn"].value_counts()
    keep = set(keep or [])
    keep |= {"Not RE"} | set(vc[vc / n >= min_frac].index)
    s = s.copy()
    s["assign2"] = np.where(s["assgn"].isin(keep), s["assgn"], "Diverse REs")
    return s


# ---------------------------------------------------------------- volcano panel
def volcano_panel(ax, axbar, g, kzfp, window_kb, axleg=None, ymax=None, claim="DOWN"):
    order = [("NS", False), ("NS", True), ("UP", False), ("DOWN", False),
             ("UP", True), ("DOWN", True)]
    for cat, hp in order:
        m = (g.cat == cat) & (g.haspeak == hp)
        if not m.any():
            continue
        sub = g[m]
        if hp:
            sz = bubble_area(sub.pscore.fillna(200) / 10.0)
            ax.scatter(sub.l2fc, sub.nl2p, s=sz, c=VOLC[(cat, hp)],
                       linewidths=0.25, edgecolors="black", alpha=0.85, zorder=3 if cat != "NS" else 2)
        else:
            ax.scatter(sub.l2fc, sub.nl2p, s=1.6, c=VOLC[(cat, hp)],
                       linewidths=0, alpha=0.85, zorder=1 if cat == "NS" else 2)
    ax.set_xlabel("Log2(FC)"); ax.set_ylabel("-Log2(padj)")
    lim = np.nanmax(np.abs(g.l2fc)) * 1.08
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.5, (ymax if ymax else np.nanmax(g.nl2p) * 1.06))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # stacked % bar -----------------------------------------------------------
    m = g.cat == claim
    tab = [[int((m & g.haspeak).sum()), int((~m & g.haspeak).sum())],
           [int((m & ~g.haspeak).sum()), int((~m & ~g.haspeak).sum())]]
    orr, pv = fisher_exact(tab)
    for i, hp in enumerate([True, False]):
        sub = g[g.haspeak == hp]
        n = max(len(sub), 1)
        bottom = 0.0
        for cat in ["UP", "NS", "DOWN"]:
            v = (sub.cat == cat).sum() / n * 100
            axbar.bar(i, v, bottom=bottom, width=0.62, color=VOLC[(cat, hp)],
                      edgecolor="none")
            bottom += v
    axbar.set_xticks([0, 1]); axbar.set_xticklabels(["Peak", "No Peak"])
    axbar.set_ylim(0, 108); axbar.set_yticks([0, 25, 50, 75, 100])
    axbar.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    axbar.set_ylabel("% genes")
    axbar.set_xlim(-0.6, 1.6)
    for sp in ("top", "right"):
        axbar.spines[sp].set_visible(False)
    axbar.text(0.5, 103, stars(pv), ha="center", va="bottom", fontsize=FS["star"])

    if axleg is not None:
        draw_volcano_legend(axleg, window_kb)
    return dict(kzfp=kzfp, claim=claim, n_genes=len(g), n_peak=int(g.haspeak.sum()),
                table=tab, odds_ratio=orr, p=pv, stars=stars(pv))


def draw_volcano_legend(ax, window_kb):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.255, 0.93, "padj", ha="right", va="center", fontsize=FS["leg"])
    ax.text(0.415, 0.93, "<0.05", ha="center", va="center", fontsize=FS["leg"])
    ax.plot([0.295, 0.535], [0.83, 0.83], color="k", lw=0.6)
    ax.text(0.655, 0.93, ">0.05", ha="center", va="center", fontsize=FS["leg"])
    ax.text(0.745, 0.93, "Peak Value (1/10)", ha="left", va="center", fontsize=FS["legs"])
    ax.text(0.375, 0.68, "Down", ha="center", va="center", fontsize=FS["leg"])
    ax.plot([0.298, 0.452], [0.585, 0.585], color="k", lw=0.6)
    ax.text(0.520, 0.68, "Up", ha="center", va="center", fontsize=FS["leg"])
    ax.plot([0.462, 0.582], [0.585, 0.585], color="k", lw=0.6)
    ax.text(0.255, 0.44, "Diff. Exp.", ha="right", va="center", fontsize=FS["leg"])
    ax.text(0.255, 0.20, f"Peak <{window_kb}Kb", ha="right", va="center", fontsize=FS["leg"])
    xs = [0.300, 0.373, 0.446, 0.519, 0.592, 0.665]
    keys = [("DOWN", True), ("DOWN", False), ("UP", True), ("UP", False),
            ("NS", True), ("NS", False)]
    for x, k in zip(xs, keys):
        ax.add_patch(Rectangle((x, 0.36), 0.050, 0.17, facecolor=VOLC[k], edgecolor="none"))
        ax.text(x + 0.025, 0.20, "+" if k[1] else "\u2013", ha="center", va="center", fontsize=FS["leg"])
    for v, y in [(60, 0.62), (40, 0.55), (20, 0.48)]:
        ax.scatter([0.800], [y], s=bubble_area(v), facecolor="#828282",
                   edgecolor="black", linewidths=0.4, alpha=0.75, zorder=3)
    for y, lab in [(0.76, "60"), (0.62, "40"), (0.48, "20")]:
        ax.plot([0.828, 0.880], [0.56, y], color="k", lw=0.4, zorder=2)
        ax.text(0.888, y, lab, ha="left", va="center", fontsize=FS["leg"])
    ax.plot([0.22, 1.0], [0.04, 0.04], color="k", lw=0.6)


# --------------------------------------------------------------- alluvial panel
def _ribbon(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha=1.0, lw=0.0, zorder=1):
    dx = (x1 - x0) * 0.42
    verts = [(x0, y0a), (x0 + dx, y0a), (x1 - dx, y1a), (x1, y1a),
             (x1, y1b), (x1 - dx, y1b), (x0 + dx, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none",
                           alpha=alpha, lw=lw, zorder=zorder))


def alluvial_panel(ax, s, node_w=0.072, gap_frac=0.013, right_labels=True,
                   status_order=("UP", "DOWN")):
    """Three-stage alluvial: Status -> distance bin -> repeat assignation.

    Every flow is tracked on the FULL (Direction, bin, assignation) triple, so a
    ribbon keeps one colour from the Status node to the Assignation node.  The
    stacking key inside the middle node is identical on the incoming and the
    outgoing side, so each band is continuous across the panel.
    """
    s = s.copy()
    dir_order = [d for d in status_order if (s.Direction == d).any()]
    bin_order = [b for b in reversed(BINS) if (s.bin == b).any()]
    # third-stage order: barycentre of the incoming distance bins, so ribbons run
    # as flat as possible and crossings are minimised
    _bo = [b for b in reversed(BINS) if (s.bin == b).any()]
    _bidx = {b: i for i, b in enumerate(_bo)}
    _tmp = s.copy(); _tmp["_i"] = _tmp.bin.map(_bidx)
    _bc = _tmp.groupby("assign2")["_i"].mean()
    asg_order = list(_bc.sort_values(kind="stable").index)
    di = {v: i for i, v in enumerate(dir_order)}
    bi = {v: i for i, v in enumerate(bin_order)}
    ai = {v: i for i, v in enumerate(asg_order)}
    FL = [(d, b, a, int(h)) for (d, b, a), h
          in s.groupby(["Direction", "bin", "assign2"]).size().items() if h > 0]
    n = len(s)
    gap = n * gap_frac
    levels = [dir_order, bin_order, asg_order]

    node = {}
    for si, lv in enumerate(levels):
        y = 0.0
        for l in lv:
            h = sum(f[3] for f in FL if f[si] == l)
            node[(si, l)] = (y, y + h)
            y += h + gap

    # offsets inside each node, one deterministic key per node ----------------
    def offsets(si, key):
        o = {}
        for l in levels[si]:
            y = node[(si, l)][0]
            for f in sorted([f for f in FL if f[si] == l], key=key):
                o[f[:3]] = y
                y += f[3]
        return o
    offA = offsets(0, lambda f: (bi[f[1]], ai[f[2]]))
    offB = offsets(1, lambda f: (di[f[0]], ai[f[2]]))
    offC = offsets(2, lambda f: (di[f[0]], bi[f[1]]))

    xs = [0.0, 0.5, 1.0]
    # large flows first so that a crossing occludes rather than blends colours
    for zi, (d, b, a, h) in enumerate(sorted(FL, key=lambda f: -f[3])):
        k = (d, b, a)
        col = ALLUV[(d, b)]
        for i in (0, 1):
            _ribbon(ax, xs[i] + node_w / 2, xs[i + 1] - node_w / 2,
                    (offA if i == 0 else offB)[k], (offA if i == 0 else offB)[k] + h,
                    (offB if i == 0 else offC)[k], (offB if i == 0 else offC)[k] + h,
                    col, alpha=0.62, lw=0.0, zorder=1 + zi * 0.01)

    for si, lv in enumerate(levels):
        for l in lv:
            y0, y1 = node[(si, l)]
            ax.add_patch(Rectangle((xs[si] - node_w / 2, y0), node_w, y1 - y0,
                                   facecolor=NODE_GREY, edgecolor="#4D4D4D", lw=0.45, zorder=4))
            if si < 2 and (y1 - y0) > n * 0.07:
                ax.text(xs[si], (y0 + y1) / 2, l, rotation=90, ha="center", va="center",
                        fontsize=FS["node"], zorder=5)
            if si == 2 and right_labels and (y1 - y0) > n * 0.02:
                ax.text(xs[si] + node_w / 2 + 0.03, (y0 + y1) / 2, l,
                        ha="left", va="center", fontsize=FS["node"], zorder=5)
    ax.set_xlim(-0.12, 1.42)
    ax.set_ylim(-n * 0.02, max(v[1] for v in node.values()) * 1.02)
    ax.set_xticks(xs)
    ax.set_xticklabels(["Status", "Distance TSS - Peak", "Assignation"])
    ax.set_ylabel("Number of genes")
    ax.tick_params(axis="x", length=0)
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    return node


def alluvial_legend(ax, ncol=3):
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    hs = [Rectangle((0, 0), 1, 1, facecolor=ALLUV[(d, b)], edgecolor="none")
          for d in ("DOWN", "UP") for b in BINS]
    ls = [f"{d}; {b}" for d in ("DOWN", "UP") for b in BINS]
    order = [0, 3, 1, 4, 2, 5] if ncol == 3 else list(range(6))
    ax.legend([hs[i] for i in order], [ls[i] for i in order], ncol=ncol, frameon=False,
              loc="upper left", bbox_to_anchor=(0.06, 1.25), handlelength=0.9, handleheight=0.9, handletextpad=0.35,
              columnspacing=1.0, labelspacing=0.35, fontsize=FS["legs"])
