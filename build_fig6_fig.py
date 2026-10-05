"""Figure 6 assembly for PLOS Genetics (7.48 x 8.70 in, all text >= 8 pt).

A  anti-HA / anti-actin blots: GFP, and full length vs dSCAN for both baits.
   Three content-tight crops of the carry-over experimental block, laid in one
   row; the published legend row is redrawn here at 8 pt with "Full lenght"
   corrected to "Full length".
B  full-length vs dSCAN proliferation series (ZNF498 left, ZNF18 right), the
   third crop of the same block, with its own letter.
C  ZNF18 recovered partners, enrichment over GFP (11 partners)
D  ZNF498 recovered partners, enrichment over GFP (12 partners)
E  partners shared between the two baits (5, of which 3 carry a SCAN domain)
F  InterPro domain enrichment among the 18 unique partners of the REPLICATED
   set (11 for ZNF18 and 12 for ZNF498 over 23 bait-partner pairs, 5 shared).
   The 27-pair / 22-protein discovery set is retained in S1 Table.
G  SCAN-domain share of each bait's partner list, as counts
H  ZNF18 and ZNF498 placed in the documented SCAN heterodimerisation network:
   BioGRID multi-validated SCAN-SCAN edges in grey, the edges added here in
   colour. Built from UniProt InterPro IPR003309 (human, reviewed) intersected
   with BIOGRID-MV-Physical-5.0.261.

Panel letters run left to right then top to bottom; A..H are drawn by
panel_A, panel_BC (C and D), panel_D (E), panel_E (F), panel_scanfrac (G)
and panel_network (H) - the function names predate the current lettering and
are kept so that the call sites in other scripts do not break.

Removed from the submitted version: the Reactome enrichment panel (with 11 and
12 partners per bait no pathway test is meaningful) and the log2 FL/dSCAN
against log2 over GFP scatter, whose content is carried by the tier colours in
C and D and by S8E Fig.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Rectangle, Circle
from plos_export import PK  # pins output to the figure's declared size

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

OUT = "out/fig6"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 8.70
C_GFP, C_FL, C_DS = "#7F7F7F", "#1F4E9C", "#E07B39"
C_SCAN = "#8B1A8B"
TIER_FC = {"A (high)": "#1F4E9C", "B (medium)": "#5B86C4",
           "C (cross-bait rescue)": "#9DB9DC", "D (detection / two-part)": "#D6D6D6"}
TIER_SHORT = {"A (high)": "A", "B (medium)": "B",
              "C (cross-bait rescue)": "C", "D (detection / two-part)": "D"}


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def panel_A(F, imgs):
    """Blots and the DeltaSCAN viability series, each cropped to its content."""
    gfp, pairs, curves = imgs
    for im, rect in ((gfp,   [0.075, 0.862, 0.106, 0.100]),
                     (pairs, [0.200, 0.862, 0.430, 0.100]),
                     (curves,[0.700, 0.862, 0.252, 0.100])):
        ax = F.add_axes(rect)
        ax.imshow(im)
        ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")
    axl = F.add_axes([0.700, 0.962, 0.252, 0.001])
    axl.set_xticks([]); axl.set_yticks([]); axl.axis("off")
    hs = [mlines.Line2D([], [], color=c, lw=1.6) for c in (C_GFP, C_FL, C_DS)]
    axl.legend(hs, ["GFP", "Full length", "$\\Delta$SCAN"], ncol=3, frameon=False,
               loc="lower left", bbox_to_anchor=(-0.02, 0.0), handlelength=1.3,
               handletextpad=0.4, columnspacing=1.2, fontsize=8)
    LT(F, 0.012, 0.986, "A")
    LT(F, 0.660, 0.986, "B")
    return None


def _partner_bars(F, rect, d, bait, xmax):
    ax = F.add_axes(rect)
    d = d.sort_values("log2 vs GFP", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    cols = [TIER_FC[t] for t in d.Tier]
    ax.barh(y, d["log2 vs GFP"], xerr=d.SE, height=0.70, color=cols,
            edgecolor="black", lw=0.4,
            error_kw=dict(elinewidth=0.6, capsize=1.2, capthick=0.6, ecolor="#3A3A3A"))
    labs = [g + (" *" if s else "") for g, s in zip(d.Gene, d["SCAN-domain protein"])]
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=8, style="italic")
    ax.tick_params(axis="y", length=0, pad=1.5)
    for tl, s in zip(ax.get_yticklabels(), d["SCAN-domain protein"]):
        if s:
            tl.set_color(C_SCAN)
    ax.set_xlabel("Log$_2$ enrichment over GFP", labelpad=1.0)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.set_title(f"{bait}  ({len(d)} partners)", fontsize=8, pad=2.0)
    ax.spines[["top", "right"]].set_visible(False)
    return ax, d


def panel_BC(F, SI):
    xmax = float((SI["log2 vs GFP"] + SI.SE).max()) * 1.10
    axB, dB = _partner_bars(F, [0.150, 0.592, 0.255, 0.200], SI[SI.Bait == "ZNF18"],
                            "ZNF18", xmax)
    axC, dC = _partner_bars(F, [0.600, 0.592, 0.255, 0.200], SI[SI.Bait == "ZNF498"],
                            "ZNF498", xmax)
    hs = [Rectangle((0, 0), 1, 1, facecolor=TIER_FC[t], edgecolor="black", lw=0.4)
          for t in TIER_FC]
    axB.legend(hs, [f"tier {TIER_SHORT[t]}" for t in TIER_FC], ncol=4, frameon=False,
               loc="lower left", bbox_to_anchor=(0.0, 1.06), handlelength=1.0,
               handleheight=1.0, handletextpad=0.4, columnspacing=1.0, fontsize=8)
    axC.text(1.0, 1.06, "* SCAN-domain protein", transform=axC.transAxes,
             ha="right", va="bottom", fontsize=8, color=C_SCAN)
    LT(F, 0.012, 0.840, "C")
    LT(F, 0.500, 0.840, "D")
    return dB, dC


def panel_D(F, SI):
    """Shared versus bait-specific partners."""
    ax = F.add_axes([0.075, 0.402, 0.390, 0.130])
    sh = sorted(SI.loc[SI["Shared between baits"], "Gene"].unique())
    only18 = sorted(set(SI[SI.Bait == "ZNF18"].Gene) - set(sh))
    only498 = sorted(set(SI[SI.Bait == "ZNF498"].Gene) - set(sh))
    scan = set(SI.loc[SI["SCAN-domain protein"], "Gene"])
    for x0, x1, col, lab in [(0.01, 0.49, C_FL, "ZNF18"),
                             (0.30, 0.99, C_DS, "ZNF498")]:
        ax.add_patch(Rectangle((x0, 0.10), x1 - x0, 0.76, facecolor=col, alpha=0.14,
                               edgecolor=col, lw=0.8))
        ax.text(x0 + 0.03 if lab == "ZNF18" else x1 - 0.03, 0.885, lab,
                ha="left" if lab == "ZNF18" else "right", va="bottom", fontsize=8,
                color=col, fontweight="bold")
    def _col(gs):
        return "\n".join(g + (" *" if g in scan else "") for g in gs)
    ax.text(0.125, 0.80, _col(only18), ha="center", va="top", fontsize=8,
            style="italic", linespacing=1.22)
    ax.text(0.390, 0.80, _col(sh), ha="center", va="top", fontsize=8,
            style="italic", fontweight="bold", linespacing=1.22)
    half = (len(only498) + 1) // 2
    ax.text(0.660, 0.80, _col(only498[:half]), ha="center", va="top", fontsize=8,
            style="italic", linespacing=1.22)
    ax.text(0.890, 0.80, _col(only498[half:]), ha="center", va="top", fontsize=8,
            style="italic", linespacing=1.22)
    ax.text(0.5, 0.015, f"{len(sh)} shared ({sum(1 for g in sh if g in scan)} "
                        f"SCAN-domain), {len(only18)} ZNF18-only, "
                        f"{len(only498)} ZNF498-only", ha="center", va="bottom",
            fontsize=8, color="#5A5A5A")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, 0.012, 0.546, "E")
    return sh, only18, only498


def panel_E(F, ip, se):
    ax = F.add_axes([0.620, 0.464, 0.290, 0.054])
    d = ip[(ip.FDR <= 0.05) & (ip.n >= 2)].sort_values(
        "FDR", ascending=False).reset_index(drop=True)
    y = np.arange(len(d))
    cols = [C_SCAN if "SCAN" in t else "#8C8C8C" for t in d.Term]
    ax.barh(y, -np.log10(d.FDR), height=0.66, color=cols, edgecolor="black", lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([t if len(t) < 34 else t[:31] + "..." for t in d.Term],
                       fontsize=8)
    ax.tick_params(axis="y", length=0, pad=1.5)
    for i, r in d.iterrows():
        ax.text(-np.log10(r.FDR) + 0.15, i, f"{int(r.n)}", va="center", ha="left",
                fontsize=8)
    ax.axvline(-np.log10(0.05), color="#B00000", lw=0.8, ls="--")
    ax.set_xlabel("-Log$_{10}$(FDR)   (bar label = proteins)", labelpad=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    r0 = se.iloc[0]
    LT(F, 0.500, 0.546, "F")
    return d


def panel_F(F, SI):
    """Enrichment over GFP against SCAN-dependence, one point per partner."""
    ax = F.add_axes([0.110, 0.148, 0.290, 0.096])
    scan = SI["SCAN-domain protein"].values
    for b, col in (("ZNF18", C_FL), ("ZNF498", C_DS)):
        m = SI.Bait == b
        ax.scatter(SI.loc[m, "log2 vs GFP"], SI.loc[m, "log2 FL/dSCAN"],
                   s=np.where(SI.loc[m, "SCAN-domain protein"], 42, 20),
                   facecolor=col,
                   edgecolor=np.where(SI.loc[m, "SCAN-domain protein"],
                                      C_SCAN, "black"),
                   linewidths=np.where(SI.loc[m, "SCAN-domain protein"], 1.1, 0.4),
                   label=b, zorder=4)
    lim = float(max(SI["log2 vs GFP"].max(), SI["log2 FL/dSCAN"].max())) * 1.12
    ax.plot([0, lim], [0, lim], lw=0.7, ls="--", color="#999999", zorder=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    for g in ("SCAND1",):
        for r in SI[SI.Gene == g].itertuples():
            ax.annotate(g, xy=(getattr(r, "_9"), getattr(r, "_11")),
                        xytext=(6, -6), textcoords="offset points", fontsize=8,
                        style="italic", color=C_SCAN if r._17 else "#3A3A3A",
                        zorder=6)
    ax.set_xlabel("Log$_2$ over GFP", labelpad=1.0)
    ax.set_ylabel("Log$_2$ FL / $\\Delta$SCAN", labelpad=1.5)
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.02, 0.02),
              handletextpad=0.25, labelspacing=0.20, markerscale=1.0, fontsize=8)
    ax.text(0.03, 0.97, "ringed = SCAN domain", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color=C_SCAN)
    ax.spines[["top", "right"]].set_visible(False)
    LT(F, 0.012, 0.260, "G")
    return SI


def panel_G(F, SI, sg):
    """Interactome: AP-MS edges with the STRING prior drawn underneath."""
    ax = F.add_axes([0.470, 0.030, 0.500, 0.212])
    baits = ["ZNF18", "ZNF498"]
    sh = sorted(SI.loc[SI["Shared between baits"], "Gene"].unique())
    o18 = sorted(set(SI[SI.Bait == "ZNF18"].Gene) - set(sh))
    o498 = sorted(set(SI[SI.Bait == "ZNF498"].Gene) - set(sh))
    scan = set(SI.loc[SI["SCAN-domain protein"], "Gene"])
    pos = {"ZNF18": (0.22, 0.88), "ZNF498": (0.78, 0.88)}
    row1 = o18 + sh                      # 11 nodes
    row2 = o498                          # 11 nodes
    for i, g in enumerate(row1):
        pos[g] = (0.035 + i * (0.93 / max(len(row1) - 1, 1)), 0.655)
    for i, g in enumerate(row2):
        pos[g] = (0.035 + i * (0.93 / max(len(row2) - 1, 1)), 0.315)
    nstr = 0
    for r in sg.itertuples():
        a, b = r.preferredName_A, r.preferredName_B
        if a in pos and b in pos:
            ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], lw=1.3,
                    color="#C9A227", alpha=0.70, zorder=1)
            nstr += 1
    napms = 0
    for r in SI.itertuples():
        a, b = r.Bait, r.Gene
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], lw=0.5,
                color="#9A9A9A", alpha=0.80, zorder=2)
        napms += 1
    for g, (xx, yy) in pos.items():
        isb = g in baits
        ax.scatter([xx], [yy], s=150 if isb else 44,
                   facecolor=(C_FL if g == "ZNF18" else C_DS) if isb
                   else (C_SCAN if g in scan else "#CFCFCF"),
                   edgecolor="black", linewidths=0.5, zorder=4)
        if isb:
            ax.text(xx, yy + 0.045, g, ha="center", va="bottom", fontsize=8,
                    style="italic", fontweight="bold", zorder=5)
        else:
            ax.text(xx, yy - 0.040, g, ha="center", va="top", fontsize=8,
                    style="italic", rotation=90,
                    color=C_SCAN if g in scan else "black", zorder=5)
    hs = [mlines.Line2D([], [], color="#9A9A9A", lw=0.9),
          mlines.Line2D([], [], color="#C9A227", lw=1.5),
          mlines.Line2D([], [], marker="o", linestyle="none", color=C_SCAN,
                        markeredgecolor="black", markeredgewidth=0.5, markersize=4)]
    ax.legend(hs, [f"AP-MS ({napms})", f"STRING prior ({nstr})", "SCAN domain"],
              ncol=3, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, 1.00),
              handlelength=1.4, handletextpad=0.4, columnspacing=1.2, fontsize=8)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0.0, 1.00)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, 0.432, 0.266, "H")
    return napms, nstr



# ---------------------------------------------------------------- new panels
def panel_scanfrac(F, SI, se, rect, letter_xy):
    """SCAN-domain share of each bait's replicated partner list."""
    ax = F.add_axes(rect)
    baits = ["ZNF18", "ZNF498"]
    y = np.arange(len(baits))
    ns, no = [], []
    for b in baits:
        d = SI[SI.Bait == b]
        ns.append(int(d["SCAN-domain protein"].sum()))
        no.append(int(len(d) - d["SCAN-domain protein"].sum()))
    ax.barh(y, ns, height=0.58, color=C_SCAN, edgecolor="black", lw=0.4,
            label="SCAN domain")
    ax.barh(y, no, left=ns, height=0.58, color="#D5D5D5", edgecolor="black",
            lw=0.4, label="other")
    for i, (a, b_) in enumerate(zip(ns, no)):
        ax.text(a / 2, i, str(a), ha="center", va="center", fontsize=8,
                color="white", fontweight="bold")
        ax.text(a + b_ / 2, i, str(b_), ha="center", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"${b}$" for b in baits], fontsize=8)
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.set_xlim(0, max(np.array(ns) + np.array(no)) * 1.10)
    ax.spines[["top", "right"]].set_visible(False)
    r0 = se.iloc[0]
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5),
              ncol=1, handlelength=1.0, handleheight=1.0, handletextpad=0.4,
              labelspacing=0.35, fontsize=8)
    ax.set_xlabel("Replicated partners", labelpad=1.0)
    F.text(0.5, rect[1] - 0.040, f"{int(r0.partners_scan)}/"
           f"{int(r0.partners_total)} replicated partners carry a SCAN domain "
           f"against {int(r0.bg_scan)}/{int(r0.bg_total)} of the detected proteome: "
           f"{r0.fold:.0f}-fold, one-sided hypergeometric P = {r0.p:.0e}",
           ha="center", va="top", fontsize=8, color="#5A5A5A")
    LT(F, *letter_xy, "G")
    return ns, no


def panel_network(F, SI, known_pairs, scan_human, rect, letter_xy, seed=4):
    """ZNF18 and ZNF498 placed in the documented SCAN heterodimerisation network."""
    import networkx as nx
    baits = ["ZNF18", "ZNF498"]
    ours = {g.upper() for g in SI.Gene} & set(scan_human)
    members = {x for p in known_pairs for x in p} | ours | set(baits)
    G = nx.Graph()
    G.add_nodes_from(members)
    for p in known_pairs:
        a, b_ = tuple(p)
        G.add_edge(a, b_, kind="known")
    for r in SI.itertuples():
        if r.Gene.upper() in ours:
            G.add_edge(r.Bait, r.Gene.upper(), kind="new")
    ax = F.add_axes(rect)
    pos = nx.spring_layout(G, seed=seed, k=2.2, iterations=800)
    # the panel is a wide strip; stretch the layout to fill it so labels separate
    px = np.array([p[0] for p in pos.values()])
    py = np.array([p[1] for p in pos.values()])
    px = (px - px.min()) / (px.max() - px.min())
    py = (py - py.min()) / (py.max() - py.min())
    pos = {n: (px[i], py[i]) for i, n in enumerate(pos)}
    kn = [(u, v) for u, v, d in G.edges(data=True) if d["kind"] == "known"]
    nw = [(u, v) for u, v, d in G.edges(data=True) if d["kind"] == "new"]
    nx.draw_networkx_edges(G, pos, edgelist=kn, edge_color="#DCDCDC", width=0.8,
                           ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=[e for e in nw if "ZNF18" in e],
                           edge_color=C_FL, width=1.8, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=[e for e in nw if "ZNF498" in e],
                           edge_color=C_DS, width=1.8, ax=ax)
    cols, sizes = [], []
    for n in G.nodes():
        if n in baits:
            cols.append("#111111"); sizes.append(300)
        elif n == "SCAND1":
            cols.append(C_SCAN); sizes.append(220)
        elif n in ours:
            cols.append("#3E8E62"); sizes.append(165)
        else:
            cols.append("#EAEAEA"); sizes.append(95)
    nx.draw_networkx_nodes(G, pos, node_color=cols, node_size=sizes,
                           edgecolors="white", linewidths=0.8, ax=ax)
    # labels pushed radially outward from the graph centroid, in POINTS,
    # so each one clears its own marker and the crowded core
    r_pt = {n: (np.sqrt(s / np.pi) + 3.0) for n, s in zip(G.nodes(), sizes)}
    cx0 = np.mean([p[0] for p in pos.values()])
    cy0 = np.mean([p[1] for p in pos.values()])
    for n, (x, yv) in pos.items():
        key = n in baits or n == "SCAND1" or n in ours
        if not key:
            continue                      # context nodes stay unlabelled
        vx, vy = x - cx0, (yv - cy0) * 0.55
        nr = np.hypot(vx, vy) or 1.0
        ux, uy = vx / nr, vy / nr
        d = r_pt[n] + 2.0
        ax.annotate(n, xy=(x, yv), xytext=(ux * d, uy * d),
                    textcoords="offset points", fontsize=8,
                    ha="left" if ux > 0.25 else ("right" if ux < -0.25 else "center"),
                    va="bottom" if uy > 0.25 else ("top" if uy < -0.25 else "center"),
                    style="italic",
                    color=("#111111" if n in baits else
                           C_SCAN if n == "SCAND1" else "#2E6B4F"),
                    fontweight="bold", zorder=6,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.12, 1.34)
    ax.set_axis_off()
    hs = [mlines.Line2D([], [], color=C_FL, lw=1.8),
          mlines.Line2D([], [], color=C_DS, lw=1.8),
          mlines.Line2D([], [], color="#DCDCDC", lw=1.2),
          mlines.Line2D([], [], marker="o", ls="", ms=5, color="#111111"),
          mlines.Line2D([], [], marker="o", ls="", ms=4.5, color=C_SCAN),
          mlines.Line2D([], [], marker="o", ls="", ms=4, color="#3E8E62"),
          mlines.Line2D([], [], marker="o", ls="", ms=3.5, color="#D0D0D0")]
    n_ctx = G.number_of_nodes() - len(ours | set(baits))
    ax.legend(hs, ["new: ZNF18 (this study)", "new: ZNF498 (this study)",
                   "known (BioGRID multi-validated)", "bait", "SCAND1",
                   "partner found here",
                   f"{n_ctx} further network members"],
              frameon=False, fontsize=8, loc="lower center",
              bbox_to_anchor=(0.5, 0.80), ncol=4, handlelength=1.3,
              handletextpad=0.4, columnspacing=1.2, labelspacing=0.25)
    LT(F, *letter_xy, "H")
    return G.number_of_nodes(), len(kn), len(nw)


def build(SI, ip, se, sg, imgs, known_pairs=None, scan_human=None):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    panel_A(F, imgs)
    dB, dC = panel_BC(F, SI)
    sh, o18, o498 = panel_D(F, SI)
    de = panel_E(F, ip, se)
    ns, no = panel_scanfrac(F, SI, se, [0.620, 0.376, 0.195, 0.036],
                            (0.500, 0.412))
    nn, nk, nw = panel_network(F, SI, known_pairs, scan_human,
                               [0.030, 0.012, 0.940, 0.300], (0.012, 0.316))
    return F, dict(shared=sh, only18=o18, only498=o498, interpro=de,
                   scan_counts=dict(zip(["ZNF18", "ZNF498"], ns)),
                   net_nodes=nn, net_known=nk, net_new=nw)


def export(F, out=OUT):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F.savefig(f"{out}/Figure_6_PLOS.png", dpi=600, **PK(F))
    F.savefig(f"{out}/Figure_6_300.png", dpi=300, **PK(F))
    F.savefig(f"{out}/Figure_6_PLOS.svg", **PK(F))
    _im = _PILI.open(f"{out}/Figure_6_300.png").convert("RGB")
    _im.save(f"{out}/Figure_6_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
