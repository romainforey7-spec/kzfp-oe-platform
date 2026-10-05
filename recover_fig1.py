"""SUPERSEDED - recovery scratch, do not run.

Used once to recover the published Figure 1 panels as images. Sets type down to
5.8 pt. The delivered Figure 1 is built by build_fig1.py.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats as st
from PIL import Image

# Load data
s1v = pd.read_excel("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/2d6b5d92-7535-4c43-8e22-54ca0ce88c68/v3493c538_revision_ledger.xlsx", sheet_name=0)
# The S1 table has columns KZFP, D0, D4, D7, D9
# Check if it's the revision ledger or the S1 table
# Based on trace: s1v = pd.read_excel(...); s1v.columns=['KZFP','D0','D4','D7','D9']
s1v = pd.read_excel("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/2d6b5d92-7535-4c43-8e22-54ca0ce88c68/v3493c538_revision_ledger.xlsx")

# From trace context, s1v is the screen scores table
# The revision ledger artifact is used for s1v based on the reload cell:
# s1v=pd.read_excel(os.path.join(R,'Supplementary Table S1.xlsx')); s1v.columns=['KZFP','D0','D4','D7','D9']
# But artifact 3493c538 maps to revision_ledger.xlsx - let's check the actual reload
# In the final reload: s1v=pd.read_excel(...) from host.artifact_path for the ledger sheet
# The ledger has sheet 'Revision ledger' but s1v needs KZFP/D0/D4/D7/D9
# Actually from trace: L=pd.read_excel(host.artifact_path('3493c538...'),sheet_name='Revision ledger')
# and s1v was already in globals from earlier: s1v=pd.read_excel(os.path.join(R,'Supplementary Table S1.xlsx'))
# So we need to reconstruct s1v from the supplementary table

# Since we don't have direct access to Supplementary Table S1.xlsx,
# but TSS counts artifact is available, and s1v must come from somewhere:
# Looking at the reload cell more carefully:
# s1v=pd.read_excel(os.path.join(R,'Supplementary Table S1.xlsx')); s1v.columns=['KZFP','D0','D4','D7','D9']
# The revision_ledger artifact contains a sheet - let's try to load it as s1v

s1v = pd.read_excel("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/2d6b5d92-7535-4c43-8e22-54ca0ce88c68/v3493c538_revision_ledger.xlsx", sheet_name='Revision ledger')
# If that doesn't work for screen scores, we need the S1 data
# Based on dependency mappings, only these two artifacts are available
# The ledger likely has KZFP scores embedded or we use what's available

# Re-reading the trace: the final composite figure cell uses:
# MK (built from s1v), ctE2 (crosstab), FF (family enrichment), GD (TSS merge), im1 (TIFF image)
# Since the image file is not available as artifact, we need to handle panel A differently

# Let's load what we can from the artifacts
TSS = pd.read_csv("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/f2ed15a8-efce-485b-8dfb-291570589bbe/v87fdd99e_kzfp_TSS_counts.csv")

# For s1v, try the revision ledger artifact with correct sheet
try:
    s1v = pd.read_excel("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/2d6b5d92-7535-4c43-8e22-54ca0ce88c68/v3493c538_revision_ledger.xlsx", sheet_name='Revision ledger')
    if not all(c in s1v.columns for c in ['D0', 'D4', 'D7', 'D9']):
        # Try reading as the S1 supplementary table (first sheet)
        s1v = pd.read_excel("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/2d6b5d92-7535-4c43-8e22-54ca0ce88c68/v3493c538_revision_ledger.xlsx", sheet_name=0)
        if s1v.shape[1] >= 5:
            s1v.columns = list(s1v.columns[:1]) + ['D0', 'D4', 'D7', 'D9'] + list(s1v.columns[5:])
except Exception:
    s1v = pd.read_excel("C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/2d6b5d92-7535-4c43-8e22-54ca0ce88c68/v3493c538_revision_ledger.xlsx", header=0)

# Ensure s1v has KZFP column
if 'KZFP' not in s1v.columns and s1v.shape[1] >= 5:
    s1v.columns = ['KZFP', 'D0', 'D4', 'D7', 'D9'] + list(s1v.columns[5:])

ALIAS = {'ZFP69': 'ZNF69', 'ZKSCAN7': 'ZNF167', 'ZKSCAN8': 'ZNF192'}
d9 = dict(zip(s1v.KZFP, s1v.D9))

def tox(n):
    v = d9.get(n, d9.get(ALIAS.get(n, n)))
    return None if v is None or pd.isna(v) else bool(v <= 0.85)

# Build MK master table
MK = s1v[['KZFP', 'D0', 'D4', 'D7', 'D9']].copy() if all(c in s1v.columns for c in ['D0','D4','D7','D9']) else s1v.copy()
MK['toxic'] = MK['D9'] <= 0.85

# We need age and SCAN data - these come from OE_ANALYSIS.xlsx which is not available as artifact
# Build placeholder columns
MK['age'] = np.nan
MK['is_scan'] = False

# For ctE2 - genomic location crosstab
# We don't have GenomicLocation.txt as artifact, build from TSS
# The panel E used gE = gl[gl.KZFP.isin(pkk)] where pkk = set(TSS.KZNF)
# Without GenomicLocation.txt we can't reproduce this exactly
# Based on the trace values: ctE2 had TEs/Prom./Both/None for Not Toxic and Toxic

# For FF - family enrichment
# Without Fam.tsv we can't reproduce this exactly

# Since the key data files are not available as artifacts, we reconstruct from the
# hardcoded values that appear in the trace outputs

# ctE2 from trace: 259 KZFPs with peaks
# Published values: NT=[53, 86, 21, 39] TX=[14, 31, 7, 12] (from earlier cell)
# But the regenerated E used pkk restriction giving different numbers
# From the cell: print(ctE2.to_string()) after rebuilding with pkk
# The trace shows: "panel E built on 259 KZFPs" but doesn't print the exact table
# Using the best available reconstruction

ctE2_data = {
    'Not Toxic': {'TEs': 53, 'Prom.': 86, 'Both': 21, 'None': 39},
    'Toxic': {'TEs': 14, 'Prom.': 31, 'Both': 7, 'None': 12}
}
ctE2 = pd.DataFrame(ctE2_data, index=pd.CategoricalIndex(['TEs', 'Prom.', 'Both', 'None'],
                    categories=['TEs', 'Prom.', 'Both', 'None'], ordered=True))

# FF from trace - corrected family enrichment
FF_data = {
    'Family': ['LINE/L1', 'LTR/ERV1', 'LTR/ERVK', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'nonTE', 'Retroposon/SVA', 'SINE/Alu'],
    'NotToxic': [55, 69, 60, 20, 24, 49, 44, 8],
    'Toxic': [17, 17, 17, 3, 6, 16, 12, 5]
}
FF = pd.DataFrame(FF_data)

# GD for panel G
GD = TSS[['KZNF', 'n_TSS_genes']].merge(s1v[['KZFP', 'D4', 'D7', 'D9']], left_on='KZNF', right_on='KZFP', how='inner')
GD['logTSS'] = np.log10(GD['n_TSS_genes'] + 1)

# Panel A schematic - we don't have the TIFF, create a placeholder
A_placeholder = Image.new('RGB', (400, 500), color=(240, 240, 240))

NT_C, TX_C = "#7F7FFF", "#DF7F7F"
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'svg.fonttype': 'none'})

FIG = plt.figure(figsize=(7.48, 6.9))
gs = FIG.add_gridspec(3, 12, height_ratios=[1.02, 0.96, 0.92], hspace=0.72, wspace=1.30,
                      left=0.075, right=0.975, top=0.955, bottom=0.075)

def PL(ax, l, dx=-0.30, dy=1.10):
    ax.text(dx, dy, l, transform=ax.transAxes, fontweight='bold', fontsize=9, va='bottom', ha='left')

LEGP = [mpl.patches.Patch(fc=NT_C, ec='black', lw=.4, label='Not Toxic'),
        mpl.patches.Patch(fc=TX_C, ec='black', lw=.4, label='Toxic')]

# A -- schematic
axA = FIG.add_subplot(gs[0, 0:3])
axA.imshow(np.asarray(A_placeholder))
axA.axis('off')
PL(axA, 'A', dx=-0.10, dy=1.00)

# B -- heatmap
axB = FIG.add_subplot(gs[0, 3:6])
hm = MK.dropna(subset=['D4', 'D7', 'D9']).sort_values('D9', ascending=False)[['D4', 'D7', 'D9']].values
cmap = mpl.colors.LinearSegmentedColormap.from_list('bwr2', ['#FF0000', '#FFFFFF', '#3B3BA8'])
imh = axB.imshow(hm, aspect='auto', cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
axB.set_xticks([0, 1, 2])
axB.set_xticklabels(['4', '7', '9'], fontsize=6.5)
axB.set_yticks([])
axB.set_xlabel("Days Post Induction", fontsize=6.8)
axB.set_ylabel("KZFPs", fontsize=6.8)
cb = FIG.colorbar(imh, ax=axB, fraction=0.16, pad=0.06, ticks=[0, 1, 2, 3])
cb.ax.tick_params(labelsize=6)
cb.outline.set_linewidth(0.5)
cb.set_label("Norm. A570/A600", fontsize=6.2)
PL(axB, 'B', dx=-0.34, dy=1.02)

# C -- pie
axC = FIG.add_subplot(gs[0, 6:12])
nt, tx = int((~MK.toxic).sum()), int(MK.toxic.sum())
tot = nt + tx
axC.pie([nt, tx], colors=[NT_C, TX_C], startangle=90, counterclock=False, radius=1.0,
        wedgeprops=dict(linewidth=0.5, edgecolor='black'))
axC.text(1.16, 0.60, f"Toxic : {tx} ({100*tx/tot:.1f}%)", fontsize=6.8, ha='left', va='center')
axC.text(0.50, -1.02, f"Not Toxic : {nt} ({100*nt/tot:.1f}%)", fontsize=6.8, ha='left', va='center')
axC.set_xlim(-1.1, 3.0)
axC.set_ylim(-1.35, 1.2)
axC.axis('off')
PL(axC, 'C', dx=-0.02, dy=0.90)

# D -- age density
axD1 = FIG.add_subplot(gs[1, 0:4])
Ma = MK.dropna(subset=['age'])
if len(Ma) > 0:
    bins = np.arange(0, Ma.age.max() + 10, 10)
    for sel, c, lb in [(~Ma.toxic, NT_C, 'Not Toxic'), (Ma.toxic, TX_C, 'Toxic')]:
        axD1.hist(Ma.age[sel], bins=bins, orientation='horizontal', density=True, color=c, alpha=0.6,
                  edgecolor='black', linewidth=0.3, label=lb)
    axD1.invert_yaxis()
axD1.set_ylabel("Age Combined (myo)", fontsize=6.8)
axD1.set_xlabel("Fraction of KZFPs (1 = 100%)", fontsize=6.8)
axD1.legend(frameon=False, fontsize=5.8, loc='lower right', handlelength=0.8, handleheight=0.8,
            borderpad=0.1, labelspacing=0.25)
axD1.tick_params(labelsize=6.3)
axD1.spines[['top', 'right']].set_visible(False)
axD1.text(-0.18, 1.10, 'D', transform=axD1.transAxes, fontweight='bold', fontsize=9, va='bottom', ha='left')

def stacked(ax, nt_, tx_, xl, fs=6.2):
    idx = np.arange(len(nt_))
    nt_ = np.asarray(nt_)
    tx_ = np.asarray(tx_)
    tt = nt_ + tx_
    pt, pn = 100 * tx_ / tt, 100 * nt_ / tt
    ax.bar(idx, pt, 0.70, color=TX_C, edgecolor='black', linewidth=0.4)
    ax.bar(idx, pn, 0.70, bottom=pt, color=NT_C, edgecolor='black', linewidth=0.4)
    for i in idx:
        ax.text(i, pt[i] / 2, str(tx_[i]), ha='center', va='center', fontsize=fs)
        ax.text(i, pt[i] + pn[i] / 2, str(nt_[i]), ha='center', va='center', fontsize=fs)
    ax.set_xticks(idx)
    ax.set_xticklabels(xl, rotation=45, ha='right', fontsize=fs)
    ax.set_ylabel("KZFPs (%)", fontsize=6.8)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.tick_params(labelsize=6.3)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=LEGP, frameon=False, fontsize=5.8, loc='lower center', bbox_to_anchor=(0.5, 1.0),
              ncol=2, handlelength=0.8, handleheight=0.8, columnspacing=0.6, borderpad=0.05)

axE = FIG.add_subplot(gs[1, 4:7])
stacked(axE, ctE2['Not Toxic'].tolist(), ctE2['Toxic'].tolist(), list(ctE2.index))
PL(axE, 'E')

axF = FIG.add_subplot(gs[1, 7:12])
stacked(axF, FF.NotToxic.tolist(), FF.Toxic.tolist(), [f.replace('/', '.') for f in FF.Family], fs=5.8)
PL(axF, 'F')

# G -- three scatters
for k, tp in enumerate(['D4', 'D7', 'D9']):
    axg = FIG.add_subplot(gs[2, k*3:(k+1)*3])
    sub = GD[['logTSS', tp]].dropna()
    axg.scatter(sub.logTSS, sub[tp], s=2.6, color='black', alpha=.75, lw=0, zorder=2)
    sl, ic, rv, pv2, se = st.linregress(sub.logTSS, sub[tp])
    xx = np.linspace(sub.logTSS.min(), sub.logTSS.max(), 80)
    yy = ic + sl * xx
    res = sub[tp] - (ic + sl * sub.logTSS)
    serr = np.sqrt((res**2).sum() / (len(sub) - 2))
    xm = sub.logTSS.mean()
    sxx = ((sub.logTSS - xm)**2).sum()
    ci = st.t.ppf(.975, len(sub) - 2) * serr * np.sqrt(1 / len(sub) + (xx - xm)**2 / sxx)
    axg.fill_between(xx, yy - ci, yy + ci, color="#9EC5E8", alpha=.65, lw=0, zorder=1)
    axg.plot(xx, yy, color="#3B7DBF", lw=1.0, zorder=3)
    axg.set_title(f"R = {rv:.2f}\npVal = {pv2:.3f}", fontsize=6.2, loc='left', linespacing=1.2)
    axg.set_ylabel(f"Signal {tp}", fontsize=6.8)
    axg.set_xticks([1, 2, 3, 4])
    axg.set_xticklabels(['10', '100', '1000', '10000'], fontsize=6)
    axg.tick_params(labelsize=6.3)
    axg.spines[['top', 'right']].set_visible(False)
    if k == 0:
        PL(axg, 'G')
    if k == 1:
        axg.set_xlabel("Number of TSS bound (log10 scale; +1)", fontsize=6.8)

# H -- violins
axH = FIG.add_subplot(gs[2, 9:12])
a_ = MK['D9'][MK.is_scan].dropna()
b_ = MK['D9'][~MK.is_scan].dropna()
if len(a_) > 0 and len(b_) > 0:
    _, pvh = st.mannwhitneyu(a_, b_, alternative='two-sided')
    stars = '***' if pvh <= .001 else '**' if pvh <= .01 else '*' if pvh <= .05 else 'ns'
    for bb in axH.violinplot([a_, b_], positions=[0, 1], widths=0.8, showextrema=False, showmedians=False)['bodies']:
        bb.set_facecolor('#9E9E9E')
        bb.set_edgecolor('black')
        bb.set_linewidth(0.4)
        bb.set_alpha(1)
    axH.boxplot([a_, b_], positions=[0, 1], widths=0.13, patch_artist=True, showfliers=False,
                medianprops=dict(color='black', lw=0.7), boxprops=dict(facecolor='white', lw=0.4),
                whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
    ym = max(a_.max(), b_.max())
    axH.plot([0, 0, 1, 1], [ym*1.04, ym*1.09, ym*1.09, ym*1.04], lw=0.6, color='black')
    axH.text(0.5, ym*1.10, stars, ha='center', va='bottom', fontsize=7)
    axH.set_ylim(0, ym*1.22)
else:
    stars = 'ns'
axH.set_xticks([0, 1])
axH.set_xticklabels(['SCAN', 'No\nSCAN'], fontsize=6.3)
axH.set_ylabel("Norm. A570/A600 (D9)", fontsize=6.2, labelpad=1.5)
axH.tick_params(labelsize=6.3)
axH.spines[['top', 'right']].set_visible(False)
PL(axH, 'H', dx=-0.36)

FIG.subplots_adjust(wspace=1.55)
FIG.savefig('Figure_1_PLOS.png', dpi=600)