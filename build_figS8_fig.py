
"""S8 Fig - selection on the ZNF18 and ZNF498 coding sequences, and the
SCAN-deletion controls (published Fig. S6, page 1).

A  per-site dN/dS (beta/alpha) along each protein, SCAN / KRAB / ZF shaded
B  dN/dS by domain group, as published (ZF, Between_ZF, Linker, SCAN, KRAB)
C  per-domain fraction of sites under SIGNIFICANT purifying selection
D  HA signal normalised on actin, full length vs dSCAN
E  replicate concordance of the SCAN-dependence measurement

Panel C exists because the published panel B cannot support the claim it is
cited for: beta/alpha is exactly 0 at 28% (ZNF18) and 48% (ZNF498) of sites,
so the medians of ZF, Between_ZF, SCAN and (ZNF498) KRAB are all 0.00 and the
only significant contrasts are Linker and KRAB being HIGHER than the rest.
The per-site likelihood-ratio test in the same table does separate the
domains, and that is what C shows.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from scipy import stats as sst
from plos_export import PK  # pins output to the figure's declared size

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

OUT = "out/figS8"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 8.70
KZS = ("ZNF18", "ZNF498")
GRPS = ["ZF", "Between_ZF", "Linker", "SCAN", "KRAB"]
GLAB = {"ZF": "ZF", "Between_ZF": "Between\nZF", "Linker": "Linker",
        "SCAN": "SCAN", "KRAB": "KRAB"}
C_DOM = {"SCAN": "#2E7D32", "KRAB": "#C2185B", "ZF": "#1F6FB4"}
C_KZ = {"ZNF18": "#1F4E9C", "ZNF498": "#E07B39"}
YCAP = 3.6


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def _spans(d):
    out = {}
    for g in ("SCAN", "KRAB"):
        p = d.pos[d.domains_uniprot == g]
        if len(p):
            out[g] = [(int(p.min()), int(p.max()))]
    zf = []
    for g in sorted({x for x in d.domains_uniprot.dropna() if x.startswith("ZF")}):
        p = d.pos[d.domains_uniprot == g]
        zf.append((int(p.min()), int(p.max())))
    if zf:
        out["ZF"] = zf
    return out


def panel_A(F, DN, rects):
    for kz, rect in zip(KZS, rects):
        d = DN[kz]
        ax = F.add_axes(rect)
        for g, blocks in _spans(d).items():
            for a, b in blocks:
                ax.add_patch(Rectangle((a - 0.5, 0), b - a + 1, YCAP,
                                       facecolor=C_DOM[g], alpha=0.16,
                                       edgecolor="none", zorder=0))
        y = d.ba.replace([np.inf, -np.inf], np.nan)
        ok = y.notna()
        hi = ok & (y > YCAP)
        ax.scatter(d.pos[ok & ~hi], y[ok & ~hi], s=2.2, color="#1A1A1A",
                   linewidths=0, zorder=3, rasterized=True)
        ax.scatter(d.pos[hi], np.full(int(hi.sum()), YCAP * 0.985), s=5,
                   marker="^", color="#1A1A1A", linewidths=0, zorder=3,
                   clip_on=False)
        ax.set_xlim(0, len(d) + 1)
        ax.set_ylim(0, YCAP)
        ax.set_ylabel("$d_N/d_S$", labelpad=1.5)
        ax.set_yticks([0, 1, 2, 3])
        ax.tick_params(length=1.8, pad=1.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(0.0, 1.02, f"${kz}$", transform=ax.transAxes, ha="left",
                va="bottom", fontsize=8, fontweight="bold")
        nz = int(ok.sum())
        ax.text(1.0, 1.02, f"{len(d)} codons; {len(d)-nz} with "
                f"$d_S$ = 0 not plotted; {int(hi.sum())} above axis",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color="#6E6E6E")
        if kz == KZS[0]:
            ax.legend(handles=[Patch(facecolor=C_DOM[g], alpha=0.5,
                                     edgecolor="none", label=g)
                               for g in ("SCAN", "KRAB", "ZF")],
                      loc="upper left", bbox_to_anchor=(0.012, 0.98), ncol=3,
                      frameon=False, handlelength=1.1, handleheight=0.9,
                      handletextpad=0.4, columnspacing=1.0, fontsize=8)
        else:
            ax.set_xlabel("Amino-acid position", labelpad=1.0)
    LT(F, 0.012, 0.992, "A")


def panel_B(F, DN, MW, rects):
    for kz, rect in zip(KZS, rects):
        d = DN[kz]
        ax = F.add_axes(rect)
        y = d.ba.replace([np.inf, -np.inf], np.nan)
        data = [y[d.grp == g].dropna().values for g in GRPS]
        vp = ax.violinplot(data, positions=range(len(GRPS)), widths=0.86,
                           showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor("#C8C8C8")
            b.set_edgecolor("#4D4D4D")
            b.set_linewidth(0.5)
            b.set_alpha(1.0)
        ax.boxplot(data, positions=range(len(GRPS)), widths=0.14,
                   showfliers=False, whis=(10, 90), patch_artist=True,
                   medianprops=dict(color="black", lw=0.8),
                   boxprops=dict(facecolor="white", edgecolor="black", lw=0.5),
                   whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
        ax.set_xticks(range(len(GRPS)))
        ax.set_xticklabels([GLAB[g] for g in GRPS], fontsize=8,
                           linespacing=1.05)
        ax.set_xlim(-0.6, len(GRPS) - 0.4)
        ax.set_ylim(-0.15, 2.55)
        ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
        ax.set_ylabel("$d_N/d_S$", labelpad=1.5)
        ax.tick_params(length=1.8, pad=1.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(0, color="#B5B5B5", lw=0.5, zorder=0)
        ax.text(0.0, 1.04, f"${kz}$", transform=ax.transAxes, ha="left",
                va="bottom", fontsize=8, fontweight="bold")
        med = [f"{np.median(v):.2f}" for v in data]
        ax.text(1.0, 1.04, "medians " + " / ".join(med), transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="#6E6E6E")
        # the one contrast the panel is cited for
        sub = MW[(MW.kzfp == kz) & (MW.g1 == "ZF") & (MW.g2 == "SCAN")]
        if len(sub):
            p = float(sub.padj.iloc[0])
            txt = "n.s." if p >= 0.05 else f"q = {p:.2g}"
            # bracket with end ticks; the label sits clear above the rule
            yb, tick = 1.92, 0.08
            ax.plot([0, 0, 3, 3], [yb - tick, yb, yb, yb - tick],
                    lw=0.6, color="black", solid_joinstyle="miter")
            ax.text(1.5, yb + 0.09, f"ZF vs SCAN  {txt}", ha="center",
                    va="bottom", fontsize=8)
    LT(F, 0.012, 0.620, "B")


def panel_C(F, SEL8, rect):
    ax = F.add_axes(rect)
    w = 0.38
    x = np.arange(len(GRPS))
    for i, kz in enumerate(KZS):
        s = SEL8[SEL8.kzfp == kz].set_index("domain").reindex(GRPS)
        ax.bar(x + (i - 0.5) * w, s.pct_purifying, width=w, color=C_KZ[kz],
               edgecolor="none", label=f"${kz}$")
        for xx, v, n in zip(x + (i - 0.5) * w, s.pct_purifying, s.n):
            ax.text(xx, v + 1.4, f"{v:.0f}", ha="center", va="bottom",
                    fontsize=8, color=C_KZ[kz])
    ax.set_xticks(x)
    ax.set_xticklabels([GLAB[g] for g in GRPS], fontsize=8, linespacing=1.05)
    ax.set_xlim(-0.6, len(GRPS) - 0.4)
    ax.set_ylim(0, 104)
    ax.set_ylabel("Codons under significant\npurifying selection (%)",
                  linespacing=1.15, labelpad=1.5)
    ax.tick_params(length=1.8, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper right", bbox_to_anchor=(1.02, 1.14), ncol=2,
              frameon=False, handlelength=1.1, handleheight=0.9,
              handletextpad=0.4, columnspacing=1.2, fontsize=8)
    ax.set_xlabel("per-site LRT, $d_N$ < $d_S$, P < 0.05", fontsize=8,
                  color="#6E6E6E", labelpad=2.0)
    LT(F, 0.012, 0.388, "C")


def panel_D(F, WB, rect):
    ax = F.add_axes(rect)
    lab, val, col = [], [], []
    for kz in KZS:
        fl = float(WB[(WB.kzfp == kz) & (WB.form == "FL")].norm.iloc[0])
        ds = float(WB[(WB.kzfp == kz) & (WB.form == "dSCAN")].norm.iloc[0])
        lab += ["FL", "$\\Delta$SCAN"]
        val += [fl, ds]
        col += [C_KZ[kz], C_KZ[kz]]
    ax.bar(range(len(val)), val, width=0.62, color=col, edgecolor="none",
           alpha=1.0)
    for i, (v, c) in enumerate(zip(val, col)):
        if i % 2:
            ax.patches[i].set_alpha(0.45)
        ax.text(i, v + 0.035, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    for i, kz in zip((0, 2), KZS):
        pct = 100 * val[i + 1] / val[i]
        ax.text(i + 0.5, max(val[i], val[i + 1]) + 0.20, f"{pct:.0f}%",
                ha="center", va="bottom", fontsize=8, color="#6E6E6E")
        ax.text(i + 0.5, -0.46, f"${kz}$", ha="center", va="top", fontsize=8,
                fontweight="bold", color=C_KZ[kz])
    ax.set_xticks(range(len(val)))
    ax.set_xticklabels(lab, fontsize=8)
    ax.set_xlim(-0.6, len(val) - 0.4)
    ax.set_ylim(0, 2.0)
    ax.set_ylabel("HA signal\n(normalised on actin)", linespacing=1.15,
                  labelpad=1.5)
    ax.tick_params(length=1.8, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(1.0, 1.02, "one blot per pair; no replicate",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color="#6E6E6E")
    LT(F, 0.512, 0.388, "D")


def panel_E(F, R2, rect, cap=7.0):
    """Replicate concordance: rep1 vs rep2 log2 FL/dSCAN for the 27 pairs."""
    ax = F.add_axes(rect)
    d = R2.copy()
    inf = ~np.isfinite(d.rep2_l2) & d.rep2_l2.notna()
    fin = d[np.isfinite(d.rep2_l2)]
    lim = max(float(fin.rep1_l2.max()), float(fin.rep2_l2.max())) * 1.10
    ax.plot([-1, lim], [-1, lim], lw=0.6, ls="--", color="#9A9A9A", zorder=1)
    ax.axhline(0, color="#C8C8C8", lw=0.5, zorder=0)
    ax.axvline(0, color="#C8C8C8", lw=0.5, zorder=0)
    for kz in KZS:
        s = fin[fin.bait == kz]
        ax.scatter(s.rep1_l2, s.rep2_l2, s=22, facecolor=C_KZ[kz],
                   edgecolor="black", linewidths=0.4, zorder=4, label="$" + kz + "$")
    # pairs detected in FL only in rep2: complete loss on dSCAN, drawn on the cap line
    io = d[inf]
    for kz in KZS:
        s = io[io.bait == kz]
        if len(s):
            ax.scatter(s.rep1_l2, np.full(len(s), lim * 0.96), s=22, marker="^",
                       facecolor=C_KZ[kz], edgecolor="black", linewidths=0.4,
                       zorder=4, clip_on=False)
    bad = fin[fin.rep2_l2 < 0].sort_values("rep1_l2")
    ax.scatter(bad.rep1_l2, bad.rep2_l2, s=38, facecolor="none",
               edgecolor="#8A8A8A", linewidths=0.8, zorder=5)
    ax.text(0.985, 0.30, "below zero, not reproduced (all tier D):\n"
            + ", ".join(bad.gene), transform=ax.transAxes, ha="right",
            va="top", fontsize=8, color="#8A8A8A", linespacing=1.25,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8), zorder=6)
    ax.set_xlabel("Replicate 1, log$_2$ full length / $\\Delta$SCAN", labelpad=1.0)
    ax.set_ylabel("Replicate 2,\nlog$_2$ FL / $\\Delta$SCAN", linespacing=1.15,
                  labelpad=1.5)
    ax.set_xlim(-0.6, lim)
    ax.set_ylim(-1.6, lim)
    ax.tick_params(length=1.8, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    # legend below the axes: the upper-left corner carries data points
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.26),
              ncol=3, handletextpad=0.3, columnspacing=1.1, labelspacing=0.22,
              markerscale=1.0, fontsize=8, borderaxespad=0.0)
    n_ok = int((fin.rep2_l2 > 0).sum()) + int(inf.sum())
    rho = sst.spearmanr(fin.rep1_l2, fin.rep2_l2)
    ax.text(1.0, 1.02, "{} of {} pairs keep the SCAN-dependent direction; "
            "Spearman $\\rho$ = {:.2f}, P = {:.0e}".format(n_ok, len(d), rho.statistic,
                                                          rho.pvalue),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color="#6E6E6E")
    ax.text(0.02, 0.62, "$\\blacktriangle$ lost entirely on $\\Delta$SCAN",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=8,
            color="#6E6E6E", bbox=dict(facecolor="white", edgecolor="none", pad=0.8), zorder=6)
    LT(F, 0.012, rect[1] + rect[3] + 0.034, "E")


def build(DN, MW, SEL8, WB, R2):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    panel_A(F, DN, ([0.105, 0.808, 0.845, 0.136],
                    [0.105, 0.636, 0.845, 0.136]))
    panel_B(F, DN, MW, ([0.105, 0.444, 0.370, 0.134],
                        [0.590, 0.444, 0.370, 0.134]))
    panel_C(F, SEL8, [0.105, 0.236, 0.330, 0.110])
    panel_D(F, WB, [0.600, 0.248, 0.230, 0.098])
    panel_E(F, R2, [0.185, 0.046, 0.600, 0.118])
    return F


def export(F, out=OUT, name="Figure_S8"):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F.savefig(f"{out}/{name}_PLOS.png", dpi=600, **PK(F))
    F.savefig(f"{out}/{name}_300.png", dpi=300, **PK(F))
    F.savefig(f"{out}/{name}_PLOS.svg", **PK(F))
    _im = _PILI.open(f"{out}/{name}_300.png").convert("RGB")
    _im.save(f"{out}/{name}_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
