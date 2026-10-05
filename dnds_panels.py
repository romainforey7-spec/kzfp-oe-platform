"""SUPERSEDED - not used by any delivered figure.

An earlier version of the per-residue dN/dS panels. Superseded by the panel functions inside build_figS8_fig.py, which draw the same content at the 8 pt floor; this file still sets type at 4.6 to 6.0 pt.
No builder in this repository imports it. Kept for the record; running it
would produce output below the 8 pt legibility floor.
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import mannwhitneyu
from itertools import combinations

DOMCOL = {"SCAN": "#D95F02", "KRAB": "#7570B3", "ZF": "#1B9E77",
          "DUF3669": "#E7298A", "Between_ZF": "#BDBDBD", "Linker": "#F0F0F0"}
GROUP_ORDER = ["ZF", "Between_ZF", "Linker"]


def prepare(df):
    d = df.rename(columns={"Unnamed: 0": "pos"}).sort_values("pos").reset_index(drop=True)
    d["beta_alpha"] = np.where(d.alpha > 0, d.beta / d.alpha, np.nan)
    dom = d.get("domains_uniprot")
    d["dom"] = dom.where(dom.notna() & (dom.astype(str).str.strip() != ""), np.nan)
    d["is_ZF"] = d.dom.fillna("").str.contains("ZF", case=False)
    d["domain_group"] = np.where(d.is_ZF, "ZF", np.where(d.dom.notna(), d.dom, "Linker"))
    zf = d.pos[d.is_ZF]
    if len(zf):
        inside = d.dom.isna() & d.pos.between(zf.min(), zf.max())
        d.loc[inside, "domain_group"] = "Between_ZF"
    d["dom_merged"] = np.where(d.is_ZF, "ZF", d.dom)
    return d


def domain_blocks(d):
    """Contiguous runs of a non-missing merged domain label."""
    out, cur, start = [], None, None
    vals = d.dom_merged.values; pos = d.pos.values
    for i, v in enumerate(vals):
        v = None if (v is None or (isinstance(v, float) and np.isnan(v))) else v
        if v != cur:
            if cur is not None:
                out.append((cur, pos[start], pos[i - 1]))
            cur, start = v, i
    if cur is not None:
        out.append((cur, pos[start], pos[-1]))
    return out


def dot_panel(ax, d, label, ymax=None):
    for name, s, e in domain_blocks(d):
        base = "ZF" if str(name).upper().startswith("ZF") else name
        ax.axvspan(s - 0.5, e + 0.5, color=DOMCOL.get(base, "#CCCCCC"), alpha=0.22,
                   lw=0, zorder=0)
        ax.text((s + e) / 2, 1.015, base, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=4.6, rotation=0, clip_on=False)
    ax.plot(d.pos, d.beta_alpha, "o", ms=1.5, color="black", mew=0, zorder=3)
    ax.axhline(1.0, color="#B00000", lw=0.5, ls=(0, (3, 2)), zorder=2)
    ax.set_xlim(d.pos.min() - 2, d.pos.max() + 2)
    ax.set_ylim(0, ymax if ymax else np.nanpercentile(d.beta_alpha, 99.5) * 1.12)
    ax.set_xlabel("Amino acid position", fontsize=5.6)
    ax.set_ylabel(r"dN/dS ($\beta/\alpha$)", fontsize=5.6)
    ax.tick_params(labelsize=5.2, length=1.8, pad=1)
    ax.text(0.995, 0.96, label, transform=ax.transAxes, ha="right", va="top",
            fontsize=6.0, style="italic")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def violin_panel(ax, d, label, ylim=(-1, 2)):
    groups = [g for g in GROUP_ORDER if (d.domain_group == g).any()]
    groups += [g for g in d.domain_group.unique() if g not in GROUP_ORDER]
    data = [d.beta_alpha[(d.domain_group == g)].dropna().values for g in groups]
    vp = ax.violinplot(data, positions=range(len(groups)), widths=1.05,
                       showextrema=False, showmedians=False)
    for b in vp["bodies"]:
        b.set_facecolor("#B3B3B3"); b.set_edgecolor("black"); b.set_lw(0.4); b.set_alpha(0.85)
    bp = ax.boxplot(data, positions=range(len(groups)), widths=0.10, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.6),
                    boxprops=dict(facecolor="white", edgecolor="black", lw=0.4),
                    whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=5.2)
    ax.set_ylim(*ylim); ax.set_ylabel(r"dN/dS ($\beta/\alpha$)", fontsize=5.6)
    ax.tick_params(labelsize=5.2, length=1.8, pad=1)
    ax.text(0.99, 0.97, label, transform=ax.transAxes, ha="right", va="top",
            fontsize=6.0, style="italic")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return groups, data


def pairwise_stats(d, gene):
    groups = [g for g in GROUP_ORDER if (d.domain_group == g).any()]
    groups += [g for g in d.domain_group.unique() if g not in GROUP_ORDER]
    rows = []
    for a, b in combinations(groups, 2):
        xa = d.beta_alpha[d.domain_group == a].dropna()
        xb = d.beta_alpha[d.domain_group == b].dropna()
        if len(xa) < 3 or len(xb) < 3:
            continue
        u, pv = mannwhitneyu(xa, xb, alternative="two-sided")
        rows.append(dict(gene=gene, group1=a, group2=b, n1=len(xa), n2=len(xb),
                         median1=round(float(xa.median()), 4), median2=round(float(xb.median()), 4),
                         U=float(u), p=pv))
    s = pd.DataFrame(rows)
    if len(s):                                   # Benjamini-Hochberg
        o = np.argsort(s.p.values); m = len(s)
        adj = np.empty(m); prev = 1.0
        for rank, i in enumerate(o[::-1]):
            prev = min(prev, s.p.values[i] * m / (m - rank))
            adj[i] = prev
        s["p_adj_BH"] = adj
        s["signif"] = np.where(s.p_adj_BH > 0.05, "ns",
                        np.where(s.p_adj_BH > 0.01, "*",
                          np.where(s.p_adj_BH > 0.001, "**", "***")))
    return s
