
"""Figure 3 assembly for PLOS Genetics (7.48 x 8.70 in, all text >= 8 pt).

Panels are lettered in reading order, down the left column then down the right.

Left  A volcano with the peak-proximity bar | B alluvial | C targets across cell
      types | D germ-cell specificity | E transcriptome-wide volcano
Right F MAGEA cluster occupancy | G MAGEA6 5' end | H 5' alignment at
      the motif | I MAGEA ages | J spermatogenesis stages | K model

Panel K states the evolutionary model the ages in I support: tandem duplication
of an ancestral MAGEA6 carried its GAGGCA element into every descendant copy,
so one ZNF257 recognition sequence now sits at the 5' end of each copy.
Feature placement audited in Fig3H_motif_gene_feature.csv: the motif is in exon 1
in 23 of 27 isoforms, the peak summit in intron 1 for MAGEA6 and MAGEA2B.

Panel K, the model, is drawn natively (model_fig3.py) rather than lifted from
the published artwork, whose lettering is about 4.5 pt. It sits last in the
right column, below J.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram
import rnaseq_panels as RP
from plos_export import PK  # pins output to the figure's declared size

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
RP.FS.update(star=8.0, leg=8.0, legs=8.0, node=8.0)

OUT = "out/fig3"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 8.70

CMAP_HM = LinearSegmentedColormap.from_list(
    "kzfp", ["#0000CD", "#8A7DD1", "#C9C3E3", "#E8E8E8", "#F0C0B4", "#E06B52", "#CC1F0A"])
C_MAG, C_ZNF = "#9B2FAE", "#D7263D"
C_DOWN = "#377EB8"
#: status of every MAGE gene, derived from the DE table by mage_status()
MAGEA_DOWN = ["MAGEA6", "MAGEA2B", "MAGEA12", "MAGEA3"]
C_UNTESTED = "#FFFFFF"


def mage_status(de):
    """DOWN / UP / UNCHANGED / 'not tested' per MAGE gene, from the DE table."""
    st = {}
    for _, r in de[de.symbol_de.astype(str).str.startswith("MAGE")].iterrows():
        st[r.symbol_de] = r.Direction
    return st


def mage_face(gene, st):
    d = st.get(gene)
    if d == "DOWN":
        return C_MAG, "black", None
    if d is None:
        return C_UNTESTED, "#6E6E6E", "////"
    return "#B8B8B8", "black", None
NC = {"A": "#33A02C", "C": "#1F78B4", "G": "#FF7F00", "T": "#E31A1C",
      "-": "#DDDDDD", "N": "#999999"}
LX = 0.092
RX = 0.555
SHORT_BIN = {"<1Kb": "<1", "1 to 10Kb": "1-10", "10 to 50Kb": "10-50"}
LEG_BIN = {"<1Kb": "<1 kb", "1 to 10Kb": "1-10 kb", "10 to 50Kb": "10-50 kb"}


def LT(F3, x, y, s):
    F3.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def hspread(vals, lo, hi, gap):
    """Push labels apart along one axis, keeping their order."""
    o = np.argsort(vals)
    out = np.array(vals, float)
    cur = lo
    for i in o:
        out[i] = max(out[i], cur)
        cur = out[i] + gap
    if out.max() > hi:
        out = out - (out.max() - hi)
    return out


# ------------------------------------------------------------------ A
def panel_A(F3, g, ymax):
    ax = F3.add_axes([LX, 0.850, 0.168, 0.096])
    axb = F3.add_axes([LX + 0.206, 0.850, 0.054, 0.096])
    s = RP.volcano_panel(ax, axb, g, "ZNF257", 1, axleg=None, ymax=ymax, claim="DOWN")
    axb.set_xticklabels(["Peak", "No peak"], rotation=90, fontsize=8)
    axb.set_yticks([0, 50, 100])
    axb.set_yticklabels(["0", "50", "100"], fontsize=8)
    axb.set_ylabel("% genes", fontsize=8, labelpad=1.0)
    axb.tick_params(axis="y", pad=1.0)
    ax.xaxis.labelpad = 1.0
    # one-column key in the inter-column gap, bubble sizes folded in (ann. 1)
    keys = [("DOWN", True), ("DOWN", False), ("UP", True),
            ("UP", False), ("NS", True), ("NS", False)]
    labs = ["Down, peak", "Down, no peak", "Up, peak",
            "Up, no peak", "n.s., peak", "n.s., no peak"]
    hs = [Rectangle((0, 0), 1, 1, facecolor=RP.VOLC[k], edgecolor="none") for k in keys]
    for v in (20, 40, 60):
        hs.append(mlines.Line2D([], [], marker="o", linestyle="none", color="#828282",
                                markeredgecolor="black", markeredgewidth=0.4,
                                markersize=np.sqrt(RP.bubble_area(v)) * 0.75))
        labs.append(f"score {v}")
    ax.legend(hs, labs, ncol=1, frameon=False, loc="upper left",
              bbox_to_anchor=(0.352, 0.952), bbox_transform=F3.transFigure,
              handlelength=1.0, handleheight=1.0,
              handletextpad=0.4, labelspacing=0.30, fontsize=8)
    LT(F3, 0.010, 0.996, "A")
    return s


# ------------------------------------------------------------------ B
def panel_B(F3, s):
    """Alluvial of status x TSS-peak distance x repeat assignation.

    The legend is placed from the measured right edge of the assignation labels,
    so it cannot ride on top of them whatever the label text turns out to be.
    """
    ax = F3.add_axes([LX, 0.668, 0.228, 0.096])
    RP.alluvial_panel(ax, s, node_w=0.085, gap_frac=0.013)
    ax.set_xticklabels(["Status", "Distance\nTSS - peak (kb)", "Assignation"])
    for t in ax.texts:                                 # annotation 6
        if t.get_rotation() == 90 and t.get_text() in SHORT_BIN:
            t.set_text(SHORT_BIN[t.get_text()])
        elif "_" in t.get_text():
            t.set_text(t.get_text().replace("_", " "))
    hs = [Rectangle((0, 0), 1, 1, facecolor=RP.ALLUV[(d, b)], edgecolor="none")
          for d in ("DOWN", "UP") for b in RP.BINS]
    ls = [f"{d[0]}{d[1:].lower()}, {LEG_BIN[b]}" for d in ("DOWN", "UP") for b in RP.BINS]
    F3.canvas.draw()
    inv = F3.transFigure.inverted()
    right = max((inv.transform(t.get_window_extent(F3.canvas.get_renderer()).corners())
                 [:, 0].max() for t in ax.texts if t.get_rotation() == 0), default=0.36)
    ax.legend(hs, ls, ncol=1, frameon=False, loc="upper left",
              bbox_to_anchor=(right + 0.010, 0.772), bbox_transform=F3.transFigure,
              handlelength=1.0, handleheight=1.0,
              handletextpad=0.4, labelspacing=0.30, fontsize=8)
    LT(F3, 0.010, 0.806, "B")


# ------------------------------------------------------------------ C
def panel_C(F3, sc, targets, gene="ZNF257", vmax=3):
    SHORT_CT = {"Distal tubular cells": "Distal tubular",
                "Basal respiratory cells": "Basal respiratory",
                "Cone photoreceptor cells": "Cone photoreceptor"}
    r = sc[sc["Gene name"] == gene][["Cell type", "nTPM"]].copy()
    r = r.sort_values(["nTPM", "Cell type"], ascending=[False, True]).reset_index(drop=True)
    hi = r.head(3)
    lo = r.sort_values(["nTPM", "Cell type"], ascending=[True, True]).head(3)
    lo = lo.sort_values(["nTPM", "Cell type"], ascending=[False, True])
    sel = pd.concat([hi, lo]).drop_duplicates("Cell type")
    cols = sel["Cell type"].tolist()
    barvals = dict(zip(sel["Cell type"], sel.nTPM))
    m = sc[sc["Gene name"].isin(targets) & sc["Cell type"].isin(cols)]
    W = m.pivot_table(index="Gene name", columns="Cell type", values="nTPM", aggfunc="mean")
    L = np.log(W.reindex(columns=cols).dropna(how="all") + 1)
    C = L.sub(L.mean(axis=1), axis=0)
    Z = linkage(C.values, method="complete", metric="euclidean")
    C = C.iloc[dendrogram(Z, no_plot=True)["leaves"]]

    x0, w = 0.168, 0.152
    axbar = F3.add_axes([x0, 0.580, w, 0.019])
    axhm = F3.add_axes([x0, 0.498, w, 0.078])
    axbar.bar(range(len(cols)), [barvals[c] for c in cols], width=0.70,
              color="#4D4D4D", edgecolor="none")
    axbar.set_xlim(-0.5, len(cols) - 0.5)
    axbar.set_xticks([])
    axbar.set_ylabel("$ZNF257$\n(nTPM)", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbar.spines[["top", "right"]].set_visible(False)

    axhm.imshow(C.values, aspect="auto", cmap=CMAP_HM, vmin=-vmax, vmax=vmax,
                interpolation="nearest")
    axhm.set_xticks(range(len(cols)))
    axhm.set_xticklabels([SHORT_CT.get(c, c) for c in cols], rotation=45, fontsize=8,
                         va="top", rotation_mode="anchor", ha="right")   # annotation 5
    axhm.tick_params(axis="x", pad=1.5, length=1.8)
    axhm.set_yticks([])
    axhm.set_ylabel("Targets", fontsize=8, labelpad=7)
    axbar.text(1.0, 1.08, f"n = {len(C)} targets", transform=axbar.transAxes,
               ha="right", va="bottom", fontsize=8)
    for sp in axhm.spines.values():
        sp.set_visible(False)
    # MAGEA rows called out to the right, spread so every name is legible (ann. 2)
    mag = sorted([(i, gn) for i, gn in enumerate(C.index) if gn.startswith("MAGE")])
    n = len(C)
    ys = np.linspace(n * 0.06, n * 0.94, len(mag)) if len(mag) > 1 else [n * 0.5]
    for (i, gn), yy in zip(mag, ys):
        axhm.annotate(gn, xy=(len(cols) - 0.4, i), xytext=(len(cols) + 0.9, yy),
                      fontsize=8, style="italic", color=C_MAG, va="center", ha="left",
                      annotation_clip=False,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color=C_MAG,
                                      shrinkA=0, shrinkB=1))
    # colour bar to the LEFT of the heatmap (annotation 8)
    cax = F3.add_axes([0.126, 0.512, 0.009, 0.042])
    cb = F3.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax), cmap=CMAP_HM),
                     cax=cax, ticks=[-vmax, 0, vmax])
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.yaxis.set_ticks_position("left")
    cax.set_ylabel("Centred\nlog(nTPM+1)", fontsize=8, labelpad=2, linespacing=1.1)
    cax.yaxis.set_label_position("left")
    LT(F3, 0.010, 0.640, "C")
    return cols


# ------------------------------------------------------------------ D
def panel_D(F3, D3, show):
    ax = F3.add_axes([LX + 0.040, 0.286, 0.250, 0.090])
    pv = ax.violinplot([D3.lg_germ, D3.lg_other], positions=[0, 1], widths=0.78,
                       showextrema=False, showmedians=False)
    for i, b in enumerate(pv["bodies"]):
        b.set_facecolor("#D9D9D9" if i == 0 else "#BFBFBF")
        b.set_edgecolor("black")
        b.set_linewidth(0.5)
        b.set_alpha(1)
    for i, v in enumerate([D3.lg_germ, D3.lg_other]):
        ax.hlines(v.median(), i - 0.19, i + 0.19, color="black", lw=1.1, zorder=4)
    for gn in show:
        col = C_ZNF if gn == "ZNF257" else C_MAG
        yg, yo = D3.loc[gn, "lg_germ"], D3.loc[gn, "lg_other"]
        ax.plot([0, 1], [yg, yo], color=col, lw=0.8, alpha=0.85, zorder=5)
        ax.scatter([0, 1], [yg, yo], s=20, color=col, edgecolor="black",
                   linewidth=0.4, zorder=6)
    ordered = sorted(show, key=lambda gn: -D3.loc[gn, "lg_germ"])
    ytxt = np.linspace(D3.lg_germ.max() * 0.72, D3.lg_germ.max() * 0.22, len(ordered))
    for gn, yy in zip(ordered, ytxt):
        col = C_ZNF if gn == "ZNF257" else C_MAG
        ax.text(0.5, yy, gn, fontsize=8, fontweight="bold", style="italic",
                color=col, ha="center", va="center", zorder=7,
                bbox=dict(boxstyle="square,pad=0.10", fc="white", ec="none", alpha=0.85))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Spermatogonia +\nspermatocytes", "All other\ncell types"], fontsize=8)
    ax.set_ylabel("log$_{10}$(nTPM + 1)", fontsize=8)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(-0.12, D3.lg_germ.max() * 1.06)
    ax.text(0.5, 1.02, f"n = {len(D3):,} genes; median {D3.lg_germ.median():.2f} vs "
            f"{D3.lg_other.median():.2f}; Wilcoxon p < 1e-300",
            transform=ax.transAxes, fontsize=8, ha="center", va="bottom",
            color="#5A5A5A")
    ax.spines[["top", "right"]].set_visible(False)
    LT(F3, 0.010, 0.404, "D")


# ------------------------------------------------------------------ G
def panel_I2(F3, age, st):
    ax = F3.add_axes([RX + 0.100, 0.462, 0.192, 0.068])
    age = sorted(age, key=lambda t: t[1])
    y = np.arange(len(age))
    faces = [mage_face(gn, st) for gn, _ in age]
    ax.barh(y, [v for _, v in age], color=[f[0] for f in faces],
            edgecolor=[f[1] for f in faces], lw=0.5, height=0.68)
    for bar, f in zip(ax.patches, faces):
        if f[2]:
            bar.set_hatch(f[2])
    ax.set_yticks(y)
    ax.set_yticklabels([gn for gn, _ in age], fontsize=8, fontstyle="italic")
    for i, (gn, v) in enumerate(age):
        ax.text(v + 2.0, i, f"{v:g}", va="center", ha="left", fontsize=8)
    ax.set_xlabel("Evolutionary age (Myr)", labelpad=1.0)
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 50, 100])
    ax.tick_params(length=1.8, pad=1.2)
    ax.spines[["top", "right"]].set_visible(False)
    hs = [Rectangle((0, 0), 1, 1, facecolor=C_MAG, edgecolor="black", lw=0.5),
          Rectangle((0, 0), 1, 1, facecolor=C_UNTESTED, edgecolor="#6E6E6E", lw=0.5,
                    hatch="////")]
    ax.legend(hs, ["repressed", "not tested"], fontsize=8, frameon=False,
              loc="upper left", bbox_to_anchor=(1.01, 1.06), handlelength=1.0,
              handleheight=1.0, handletextpad=0.4, labelspacing=0.22,
              title="On ZNF257 OE", title_fontsize=8, alignment="left")
    LT(F3, 0.512, 0.548, "I")


# ------------------------------------------------------------------ E
def panel_E(F3, de, label_genes):
    ax = F3.add_axes([LX + 0.030, 0.062, 0.262, 0.112])
    d = de.copy()
    d["fc"] = pd.to_numeric(d.foldChange, errors="coerce")
    d["padj"] = pd.to_numeric(d.padj, errors="coerce")
    d = d[d.fc.notna() & d.padj.notna()]
    d["l2fc"] = np.sign(d.fc) * np.log2(np.abs(d.fc))
    d["nl10p"] = -np.log10(d.padj.clip(lower=1e-300))
    COL = {"UP": "#E41A1C", "UNCHANGED": "#999999", "DOWN": C_DOWN}
    for k in ("UNCHANGED", "UP", "DOWN"):
        sub = d[d.Direction == k]
        ax.scatter(sub.l2fc, sub.nl10p, s=1.8, c=COL[k], linewidths=0, alpha=0.8,
                   rasterized=True, zorder=2 if k == "UNCHANGED" else 3,
                   label=f"{k.capitalize()} ({len(sub):,})")
    lab = d[d.symbol_de.isin(label_genes)].sort_values("padj").reset_index(drop=True)
    lim = max(np.nanpercentile(np.abs(d.l2fc), 99.9) * 1.12,
              float(np.abs(lab.l2fc).max()) * 1.10)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.5, d.nl10p.max() * 1.30)
    ax.scatter(lab.l2fc, lab.nl10p, s=34, facecolor=C_DOWN, edgecolor="black",
               linewidths=0.6, zorder=7)
    ty = [0.985, 0.905, 0.825]
    for i in range(len(lab)):
        ax.annotate(lab.symbol_de.iloc[i], xy=(lab.l2fc.iloc[i], lab.nl10p.iloc[i]),
                    xycoords="data", xytext=(0.030, ty[i]), textcoords="axes fraction",
                    fontsize=8, style="italic", fontweight="bold", color=C_DOWN,
                    ha="left", va="center", zorder=9,
                    arrowprops=dict(arrowstyle="-", lw=0.7, color=C_DOWN,
                                    shrinkA=2, shrinkB=4,
                                    connectionstyle="arc3,rad=0.12"))
    ax.set_xlabel("Log2(FC)", labelpad=1.0)
    ax.set_ylabel("-Log10(padj)")
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.04, 1.04),
              handletextpad=0.25, labelspacing=0.22, markerscale=3.2, fontsize=8,
              borderaxespad=0.1)
    ax.spines[["top", "right"]].set_visible(False)
    LT(F3, 0.010, 0.196, "E")
    return lab[["symbol_de", "l2fc", "padj"]]


# ------------------------------------------------------------------ F
def panel_F(F3, bw, bg, MB, pk, a, b, st, nb=1400):
    ax = F3.add_axes([RX + 0.052, 0.884, 0.298, 0.048])
    axg = F3.add_axes([RX + 0.052, 0.846, 0.298, 0.034])
    xs = np.linspace(a, b, nb)
    v = np.nan_to_num(bw.values("chrX", a, b, nbins=nb, agg="mean")) / bg
    ax.fill_between(xs, 0, v, color="#2B2B2B", lw=0, rasterized=True)
    ax.set_xlim(a, b)
    ax.set_ylim(0, v.max() * 1.22)
    ax.set_xticks([])
    ax.set_ylabel("ZNF257 (fold)", fontsize=8, labelpad=2)
    ax.set_title(f"MAGEA cluster   chrX:{a/1e6:.3f}\u2013{b/1e6:.3f} Mb (hg19)",
                 fontsize=8, pad=2.0)
    ax.spines[["top", "right"]].set_visible(False)
    # peak summits, not peak starts (annotation 7)
    for p in pk.itertuples():
        ax.plot([p.summit], [v.max() * 1.12], marker="v", ms=3.4, color="#B00000",
                clip_on=False, zorder=5)
    span = b - a
    genes = MB[(MB[0] == "chrX") & MB[1].between(a, b)].sort_values(1).reset_index(drop=True)
    mids = ((genes[1] + genes[2]) / 2).values
    for i, r_ in genes.iterrows():
        fc, ec, ha_ = mage_face(r_[3], st)
        col = fc if fc != C_UNTESTED else "#6E6E6E"
        axg.add_patch(Rectangle((r_[1], 1.80), max(r_[2] - r_[1], span * 0.004), 0.40,
                                facecolor=fc, edgecolor=ec, lw=0.5, hatch=ha_))
        yy = 1.2 if i % 2 == 0 else 0.4
        axg.annotate(r_[3], xy=(mids[i], 1.78), xytext=(mids[i], yy),
                     ha="center", va="center", fontsize=8, style="italic", color=col,
                     arrowprops=dict(arrowstyle="-", lw=0.5, color=col,
                                     shrinkA=1, shrinkB=1))
    axg.set_xlim(a, b)
    axg.set_ylim(0, 2.4)
    axg.set_yticks([])
    axg.set_xticks([])
    axg.set_xticks([]); axg.set_yticks([])
    axg.axis("off")
    axg.text(0.0, -0.08, f"{len(pk)} ZNF257 peaks (\u25bc)",
             transform=axg.transAxes, ha="left", va="top", fontsize=8, color="#B00000")
    LT(F3, 0.512, 0.952, "F")


# ------------------------------------------------------------------ H
def panel_G2(F3, bw, bg, MB, summit, nb=900, half=900):
    gi = MB[(MB[0] == "chrX") & (MB[3] == "MAGEA6")].iloc[0]
    tss = int(gi[1]) if gi[5] == "+" else int(gi[2])
    a, b = tss - half, tss + half
    ax = F3.add_axes([RX + 0.052, 0.746, 0.298, 0.048])
    axg = F3.add_axes([RX + 0.052, 0.718, 0.298, 0.022])
    xs = np.linspace(a, b, nb)
    v = np.nan_to_num(bw.values("chrX", a, b, nbins=nb, agg="mean")) / bg
    ax.fill_between(xs, 0, v, color="#2B2B2B", lw=0, rasterized=True)
    ax.plot([summit], [v.max() * 1.12], marker="v", ms=3.4, color="#B00000",
            clip_on=False, zorder=5)
    ax.set_xlim(a, b)
    ax.set_ylim(0, v.max() * 1.24)
    ax.set_xticks([])
    ax.set_ylabel("ZNF257 (fold)", fontsize=8, labelpad=2)
    ax.set_title(f"$MAGEA6$ 5′ end, peak {v.max():.0f}x genome mean",
                 fontsize=8, pad=6.0, color="#B00000", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    span = b - a
    axg.add_patch(Rectangle((int(gi[1]), 1.05), int(gi[2]) - int(gi[1]), 0.45,
                            facecolor=C_MAG, edgecolor="none"))
    axg.plot([tss, tss], [0.90, 1.65], lw=1.1, color=C_MAG, solid_capstyle="butt")
    axg.text(tss - span * 0.020, 1.62, "MAGEA6", ha="right", va="center",
             fontsize=8, style="italic", color=C_MAG)
    axg.set_xlim(a, b)
    axg.set_ylim(0, 1.80)
    axg.set_yticks([])
    axg.set_xticks([a, tss, b])
    axg.set_xticklabels([f"-{half/1000:.1f} kb", "TSS", f"+{half/1000:.1f} kb"],
                        fontsize=8)
    axg.tick_params(axis="x", length=1.8, pad=1.0)
    axg.spines[["top", "right", "left"]].set_visible(False)
    LT(F3, 0.512, 0.818, "G")
    return float(v.max()), (a, b)


# ------------------------------------------------------------------ I
def panel_H2(F3, AL, ident, order, mstart, half=11):
    """Promoter alignment restricted so an 8 pt base fits inside its cell."""
    c0, c1 = mstart - half, mstart + half + len("GAGGCA")
    Wn = c1 - c0
    axA = F3.add_axes([RX + 0.066, 0.604, 0.288, 0.058])
    axB = F3.add_axes([RX + 0.066, 0.580, 0.288, 0.020])
    for r, name in enumerate(order):
        s = AL[name][c0:c1]
        for j, chb in enumerate(s):
            axA.add_patch(Rectangle((j, len(order) - 1 - r), 1, 1,
                                    facecolor=NC.get(chb, "#999999"), edgecolor="none"))
            axA.text(j + 0.5, len(order) - 1 - r + 0.5, chb, ha="center", va="center",
                     fontsize=8, color="white", fontweight="bold")
    axA.set_xlim(0, Wn)
    axA.set_ylim(0, len(order))
    axA.set_yticks([len(order) - 1 - r + 0.5 for r in range(len(order))])
    axA.set_yticklabels([o.replace("MAGEA", "A") for o in order], fontsize=8,
                        style="italic")
    axA.set_ylabel("$MAGEA$", fontsize=8, labelpad=2.0)
    axA.set_xticks([])
    axA.tick_params(axis="y", length=0, pad=2)
    for sp in axA.spines.values():
        sp.set_linewidth(0.5)
    axB.bar(np.arange(Wn) + 0.5, ident[c0:c1], width=1.0, color="#404040", edgecolor="none")
    axB.set_xlim(0, Wn)
    axB.set_ylim(0, 1.05)
    axB.set_yticks([0, 1])
    axB.set_yticklabels(["0", "1"], fontsize=8)
    axB.set_ylabel("Identity", fontsize=8, labelpad=1.5)
    axB.set_xticks([0.5, half + 3, Wn - 0.5])
    axB.set_xticklabels([f"-{half}", "0", f"+{half}"], fontsize=8)
    axB.set_xlabel("Position relative to GAGGCA (bp)", labelpad=1.0)
    for sp in ("top", "right", "bottom"):
        axB.spines[sp].set_visible(False)
    mx = mstart - c0
    axA.add_patch(Rectangle((mx, 0), len("GAGGCA"), len(order), facecolor="none",
                            edgecolor="red", linewidth=1.2, zorder=6, clip_on=False))
    axA.text(mx + 3, len(order) + 0.14, "GAGGCA", ha="center", va="bottom", fontsize=8,
             color="red", fontweight="bold")
    LT(F3, 0.512, 0.682, "H")
    return Wn


# ------------------------------------------------------------------ K
def panel_K(F3, sc, krows, stages):
    genes = [g for _, gl in krows for g in gl]
    sub = sc[sc["Gene name"].isin(genes) & sc["Cell type"].isin(stages)]
    M = sub.pivot_table(index="Gene name", columns="Cell type", values="nTPM",
                        aggfunc="mean").reindex(columns=stages)
    Lg = np.log10(M + 1)
    Zs = Lg.sub(Lg.mean(axis=1), axis=0).div(Lg.std(axis=1).replace(0, np.nan), axis=0)
    rows = [g for _, gl in krows for g in gl]
    Zs = Zs.reindex(rows)
    STG = ["Spg", "Spc", "Early\nspd", "Late\nspd"]
    axh = F3.add_axes([RX + 0.090, 0.252, 0.116, 0.148])
    axl = F3.add_axes([RX + 0.268, 0.252, 0.082, 0.148])
    axh.imshow(Zs.values, aspect="auto", cmap=CMAP_HM, vmin=-1.6, vmax=1.6,
               interpolation="nearest")
    axh.set_xticks(range(len(stages)))
    axh.set_xticklabels([s.replace("\n", " ") for s in STG], fontsize=8,
                        rotation=45, ha="right", va="top", rotation_mode="anchor")
    axh.set_yticks(range(len(rows)))
    axh.set_yticklabels(rows, fontsize=8, style="italic")
    axh.tick_params(axis="y", length=0, pad=1.5)
    axh.tick_params(axis="x", length=1.8, pad=1.5)
    for sp in axh.spines.values():
        sp.set_visible(False)
    y0 = 0
    for lab, gl in krows:
        axh.plot([-0.75, -0.75], [y0 - 0.40, y0 + len(gl) - 0.60], lw=2.0,
                 color=C_MAG if lab in ("MAGEs", "ZNF257") else "#8C8C8C",
                 clip_on=False, solid_capstyle="butt")
        y0 += len(gl)
    # vertical colour bar to the LEFT of the heatmap, clear of the gene names
    # (annotation 8, same treatment as panel C)
    cax = F3.add_axes([RX - 0.030, 0.286, 0.009, 0.056])
    cb = F3.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -1.6, 1.6), cmap=CMAP_HM),
                     cax=cax, ticks=[-1.6, 0, 1.6])
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.yaxis.set_ticks_position("left")
    cax.set_ylabel("Row z-score", fontsize=8, labelpad=2)
    cax.yaxis.set_label_position("left")
    for gn, col in [("MAGEA3", C_MAG), ("MAGEA6", "#C77CFF"), ("ZNF257", C_ZNF)]:
        axl.plot(range(len(stages)), Lg.loc[gn, stages].values, marker="o", ms=3.0,
                 lw=1.0, color=col, label=gn)
    axl.set_xticks(range(len(stages)))
    axl.set_xticklabels([s.replace("\n", " ") for s in STG], fontsize=8,
                        rotation=90)
    axl.set_ylabel("log$_{10}$(nTPM+1)", fontsize=8, labelpad=0.5)
    axl.tick_params(length=1.8, pad=1.5)
    # legend above the axes so it cannot sit on the curves (annotation 9)
    axl.legend(frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(0.00, 1.00),
               handlelength=1.0, handletextpad=0.35, labelspacing=0.22, ncol=1)
    axl.spines[["top", "right"]].set_visible(False)
    LT(F3, 0.512, 0.428, "J")
    return M


def build(g, s, sc, targets, D3, show, de, label_genes, bw, bg, MB, pk, fwin,
          summit, AL, ident, order, mstart, age, krows, stages, ymaxA):
    """Panels run in reading order: A-E down the left, F-K down the right."""
    import model_fig3 as M3
    plt.close("all")
    st = mage_status(de)
    F3 = plt.figure(figsize=(SHEET_W, SHEET_H))
    stA = panel_A(F3, g, ymaxA)
    panel_B(F3, s)
    cols = panel_C(F3, sc, targets)
    panel_D(F3, D3, show)
    labE = panel_E(F3, de, label_genes)
    panel_F(F3, bw, bg, MB, pk, *fwin, st)
    hmax, hwin = panel_G2(F3, bw, bg, MB, summit)
    nI = panel_H2(F3, AL, ident, order, mstart)
    panel_I2(F3, age, st)
    MK = panel_K(F3, sc, krows, stages)
    axK = F3.add_axes([RX + 0.020, 0.030, 0.418, 0.158])
    M3.draw_model(axK)
    LT(F3, 0.512, 0.200, "K")
    return F3, stA, cols, labE, hmax, MK, nI, st


def export(F3, out=OUT):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F3.savefig(f"{out}/Figure_3_PLOS.png", dpi=600, **PK(F3))
    F3.savefig(f"{out}/Figure_3_300.png", dpi=300, **PK(F3))
    F3.savefig(f"{out}/Figure_3_PLOS.svg", **PK(F3))
    _im = _PILI.open(f"{out}/Figure_3_300.png").convert("RGB")
    _im.save(f"{out}/Figure_3_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
