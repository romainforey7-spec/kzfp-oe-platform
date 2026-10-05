"""SUPERSEDED - do not run.

An earlier draft of the S1 Fig layout. It sets type as small as 3.2 pt (261 text
elements below the 8 pt floor) and it writes to the same filenames as the current
builder, so running it overwrites the compliant figure with a non-compliant one.
The delivered S1 Fig is built by build_figS1_fig.py. Kept only as a record of the
panel arrangement that was tried first.
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, FancyArrowPatch
from PIL import Image
import pypdfium2 as pdfium
from scipy.cluster import hierarchy as hc
from scipy.spatial.distance import pdist

R  = r'C:\Users\forey\Desktop\20260924_For_Revision'
TX = os.path.join(R, '20250110_TE_Fam')
P33 = os.path.join(R, '20250129_SP_CR033_A09-H12_plate3')
OUT = 'out/figS1'
os.makedirs(OUT, exist_ok=True)

mpl.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
                     'xtick.major.size': 2.4, 'ytick.major.size': 2.4,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42})

TOX_THR = 0.85
NT_C, TX_C = '#7F7FFF', '#DF7F7F'
ALIAS_S1 = {'ZFP69': 'ZNF69', 'ZKSCAN7': 'ZNF167', 'ZKSCAN8': 'ZNF192'}

# ---------------------------------------------------------------- source data
s1v = pd.read_excel(os.path.join(R, 'Supplementary Table S1.xlsx'))
s1v.columns = ['KZFP', 'D0', 'D4', 'D7', 'D9']
d9 = dict(zip(s1v.KZFP, s1v.D9))


def tox(n):
    v = d9.get(n, d9.get(ALIAS_S1.get(n, n)))
    return None if v is None or pd.isna(v) else bool(v <= TOX_THR)


def presto_blocks(path, sheet):
    d = pd.read_excel(path, sheet_name=sheet, header=None)
    c0 = d[0].astype(str)
    i570 = [i for i, v in enumerate(c0) if '570' in v]
    i600 = [i for i, v in enumerate(c0) if '600' in v]
    out = []
    for a, b in zip(i570, i600):
        m5 = d.iloc[a + 2:a + 10, 2:14].apply(pd.to_numeric, errors='coerce').values
        m6 = d.iloc[b + 2:b + 10, 2:14].apply(pd.to_numeric, errors='coerce').values
        out.append((m5, m6))
    return out


# panel B: CR033 plate 3, library column 9, D9, 3 h read (matches published panel)
B_WB = os.path.join(P33, 'presto', '20250212_presto_CR033_D9', 'sumup_presto_CR033_D9.xlsx')
m5, m6 = presto_blocks(B_WB, '3h')[0]
rt = m5 / m6
blank_row = np.nanmean([rt[0], rt[7]], axis=0)
WELL = rt - blank_row                                   # 8 rows (A-H) x 12 columns
SAMP = ['LacZ.9_P3', 'ZFP14', 'ZNF583', 'ZNF182', 'ZNF200_CO',
        'GFP.9_P3', 'ZNF600_CO', 'ZNF493', 'ZNF726', 'ZNF735_CO']
NODOX, PLUSDOX = [1, 2, 3], [4, 5, 6]
nd_m = np.nanmean(WELL[NODOX, 1:11], axis=0); nd_s = np.nanstd(WELL[NODOX, 1:11], axis=0, ddof=1)
dx_m = np.nanmean(WELL[PLUSDOX, 1:11], axis=0); dx_s = np.nanstd(WELL[PLUSDOX, 1:11], axis=0, ddof=1)
norm_nd = WELL[NODOX, 1:11] / nd_m
norm_dx = WELL[PLUSDOX, 1:11] / nd_m
nn_m, nn_s = norm_nd.mean(axis=0), norm_nd.std(axis=0, ddof=1)
nx_m, nx_s = norm_dx.mean(axis=0), norm_dx.std(axis=0, ddof=1)

# panel D: genomic distribution of peaks
gl = pd.read_csv(os.path.join(TX, 'GenomicLocation.txt'), sep='\t')
gl['TEs'] = gl['TEs.perc']
gl['TSSp'] = gl['promoter.perc']
gl['Other'] = gl[['upstreamTSS.perc', '5UTR.perc', '3UTR.perc', 'exons.perc', 'introns.perc']].sum(axis=1)
CAND4 = ['ZNF257', 'ZNF498', 'ZNF18', 'ZNF43']
oth = gl[~gl.KZFP.isin(CAND4)]
N_OTHER = len(oth)


def star(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.'


DROWS = []
for lab in ['ZNF257', '__OTHER__', 'ZNF498', 'ZNF18', 'ZNF43']:
    if lab == '__OTHER__':
        DROWS.append(dict(label=f'Every other\nKZFP (n = {N_OTHER})', TEs=oth.TEs.mean(),
                          TSS=oth.TSSp.mean(), Other=oth.Other.mean(), mark='', band=None))
    else:
        r = gl[gl.KZFP == lab].iloc[0]
        if r['TEs.pval'] < r['promoter.pval']:
            mk, band = star(r['TEs.pval']), 'TEs'
        else:
            mk, band = star(r['promoter.pval']), 'TSS'
        DROWS.append(dict(label=lab, TEs=r.TEs, TSS=r.TSSp, Other=r.Other, mark=mk, band=band))
DD = pd.DataFrame(DROWS)

TSS = pd.read_csv('kzfp_TSS_counts.csv')
_tk = {'ZNF257': 'ZNF257', 'ZNF498': 'ZNF498', 'ZNF18': 'ZNF18', 'ZNF43': 'ZNF43'}
oTS = TSS[~TSS.KZNF.isin(_tk.values())]
NTSS = [int(TSS.loc[TSS.KZNF == 'ZNF257', 'n_TSS_genes'].iloc[0]), oTS.n_TSS_genes.mean(),
        int(TSS.loc[TSS.KZNF == 'ZNF498', 'n_TSS_genes'].iloc[0]),
        int(TSS.loc[TSS.KZNF == 'ZNF18', 'n_TSS_genes'].iloc[0]),
        int(TSS.loc[TSS.KZNF == 'ZNF43', 'n_TSS_genes'].iloc[0])]

# panel E: evolutionary age
mg = pd.read_excel(os.path.join(TX, 'OE_ANALYSIS.xlsx'), sheet_name='merged_OE_Tables').drop_duplicates('Row.Labels')
agemap = mg.set_index('Row.Labels')['age_combined'].to_dict()
AGE = s1v.KZFP.map(agemap).dropna()
CAND_AGE = [('ZNF257', 29.1, -11), ('ZNF43', 43.1, 11), ('ZNF498', 105.0, 0), ('ZNF18', 163.7, 0)]

# panel F: per-KZFP TE-family enrichment among toxic KZFPs
def clean(fn):
    fn = re.sub(r'_(pubM|n)?_?vs_293T.*$', '', fn)
    for _ in range(2):
        fn = re.sub(r'_(rep\d?|lowconf|n_lowconf|n|pubM|rep_n|rep_pubM)$', '', fn)
    return fn


fam = pd.read_csv(os.path.join(TX, 'enrich_kzfp_perFam', 'Fam.tsv'), sep='\t')
fam['kzfp2'] = fam.filename.map(clean)
fam['tox2'] = fam.kzfp2.map(tox)
FAM8 = ['Retroposon/SVA', 'LTR/ERV1', 'LTR/ERVK', 'SINE/Alu', 'LTR/ERVL-MaLR', 'LTR/ERVL', 'nonTE', 'LINE/L1']
fw = fam[fam.tox2 == True].copy()
Mx = fw.pivot_table(index='kzfp2', columns='subfam_name', values='padj.hypergeom.reg',
                    aggfunc='min').reindex(columns=FAM8)
Bin = (Mx <= 0.05).astype(int)
Bin = Bin[Bin.sum(axis=1) > 0]
rl = hc.linkage(pdist(Bin.values, metric='jaccard'), method='average')
cl = hc.linkage(pdist(Bin.values.T, metric='jaccard'), method='average')
Bo = Bin.iloc[hc.leaves_list(rl), hc.leaves_list(cl)]

# panel A: plate photographs, carried over from the published supplement
Image.MAX_IMAGE_PIXELS = None
pdf = pdfium.PdfDocument(os.path.join(R, 'Supplementary Figures.pdf'))
PG = pdf[0].render(scale=8.0).to_pil()
PW, PH = PG.size
A_BOX = (0.1195, 0.1955, 0.3765, 0.4140)      # full width: keeps the +/-Dox signs and the right 'Blank'
imA = PG.crop((int(A_BOX[0] * PW), int(A_BOX[1] * PH), int(A_BOX[2] * PW), int(A_BOX[3] * PH)))
# the published panel letters 'A' and 'B' fall inside this box (user annotation 1) -- white them out
aw, ah = imA.size
imA.paste((255, 255, 255), (0, 0, int(0.145 * aw), int(0.090 * ah)))
imA.paste((255, 255, 255), (int(0.875 * aw), 0, aw, int(0.090 * ah)))
imA.save(f'{OUT}/FigS1A_carryover.png')


# ---------------------------------------------------------------- figure
FW, FH = 7.48, 8.70
FS = plt.figure(figsize=(FW, FH))
TAB_CMAP = LinearSegmentedColormap.from_list('gyr', ['#63BE7B', '#FFEB84', '#F8696B'])
GREY, BLACK = '#A6A6A6', '#1A1A1A'

def letter(x, y, s):
    FS.text(x, y, s, fontweight='bold', fontsize=11, va='top', ha='left')

# ---- A : plate photographs -------------------------------------------
axA = FS.add_axes([0.030, 0.700, 0.255, 0.278]); axA.imshow(np.asarray(imA)); axA.axis('off')
letter(0.008, 0.992, 'A')

# ---- B : PrestoBlue quantification -----------------------------------
letter(0.300, 0.992, 'B')
TBX0, TBX1 = 0.362, 0.612
axBt = FS.add_axes([TBX0, 0.880, TBX1 - TBX0, 0.072])
axBt.imshow(WELL, aspect='auto', cmap=TAB_CMAP, vmin=0, vmax=np.nanmax(WELL))
for i in range(8):
    for j in range(12):
        v = WELL[i, j]
        axBt.text(j, i, f'{0.0 if abs(v) < 0.05 else v:.1f}', ha='center', va='center', fontsize=3.2)
axBt.set_xticks(range(12)); axBt.set_xticklabels([''] + SAMP + [''], rotation=90, fontsize=4.2)
axBt.xaxis.set_ticks_position('top')
axBt.set_yticks(range(8)); axBt.set_yticklabels(list('ABCDEFGH'), fontsize=4.2)
axBt.set_xticks(np.arange(-.5, 12, 1), minor=True); axBt.set_yticks(np.arange(-.5, 8, 1), minor=True)
axBt.grid(which='minor', color='white', lw=0.5); axBt.tick_params(which='minor', length=0)
axBt.tick_params(length=1.0, pad=1.0)
for sp in axBt.spines.values(): sp.set_linewidth(0.5)
for r0, r1, lb in [(1, 3, '\u2212Dox'), (4, 6, '+Dox')]:
    yy = 0.880 + 0.072 * (1 - (r0 + r1 + 1) / 16)
    FS.text(0.3465, yy, lb, fontsize=4.4, va='center', ha='center', rotation=90)
FS.text((TBX0 + TBX1) / 2, 0.874, 'Per-well A570/A600 (blank-corrected)', fontsize=5.0, ha='center', va='top')

X = np.arange(10); w = 0.38
ekw = dict(lw=0.5, capsize=1.0, capthick=0.5)
axB1 = FS.add_axes([TBX0, 0.793, TBX1 - TBX0, 0.058])
axB1.bar(X - w/2, nd_m, w, yerr=nd_s, color=GREY, ec='black', lw=0.4, error_kw=ekw)
axB1.bar(X + w/2, dx_m, w, yerr=dx_s, color=BLACK, ec='black', lw=0.4, error_kw=ekw)
axB1.axhline(1, color='black', lw=0.4, ls=':')
axB1.set_ylim(0, 2); axB1.set_yticks([0, 1, 2])
axB1.set_ylabel('A570/A600', fontsize=6.0, labelpad=1.5)
axB1.set_xlim(-0.65, 9.65); axB1.set_xticks(X); axB1.set_xticklabels([])
axB1.tick_params(labelsize=5.4, pad=1.0); axB1.spines[['top', 'right']].set_visible(False)
axB1.legend(handles=[Patch(fc=GREY, ec='black', lw=.4, label='\u2212Dox'),
                     Patch(fc=BLACK, ec='black', lw=.4, label='+Dox')],
            frameon=False, fontsize=5.0, ncol=2, loc='lower right', bbox_to_anchor=(1.01, 0.97),
            handlelength=0.8, handleheight=0.8, handletextpad=0.3, columnspacing=0.8)

axB2 = FS.add_axes([TBX0, 0.718, TBX1 - TBX0, 0.058])
axB2.bar(X - w/2, nn_m, w, yerr=nn_s, color=GREY, ec='black', lw=0.4, error_kw=ekw)
axB2.bar(X + w/2, nx_m, w, yerr=nx_s, color=BLACK, ec='black', lw=0.4, error_kw=ekw)
axB2.axhline(1, color='black', lw=0.4, ls=':')
axB2.axhline(TOX_THR, color='red', lw=0.7, ls='--')
axB2.set_ylim(0, 2); axB2.set_yticks([0, 1, 2])
axB2.set_ylabel('Signal\n(Normalized)', fontsize=6.0, labelpad=1.5)
axB2.set_xlim(-0.65, 9.65); axB2.set_xticks(X)
axB2.set_xticklabels(SAMP, rotation=45, ha='right', fontsize=4.8)
axB2.set_xlabel('KZFPs (1 plate \u2013 10 candidates)', fontsize=6.0, labelpad=0.8)
axB2.tick_params(labelsize=5.4, pad=1.0); axB2.spines[['top', 'right']].set_visible(False)
from matplotlib.lines import Line2D
axB2.legend(handles=[Line2D([0], [0], color='red', lw=0.8, ls='--', label=f'Toxicity cut-off ({TOX_THR})')],
            frameon=False, fontsize=4.6, loc='upper right', bbox_to_anchor=(1.01, 1.04),
            handlelength=1.6, handletextpad=0.35)
FS.patches.append(FancyArrowPatch((0.620, 0.805), (0.620, 0.745), transform=FS.transFigure,
                                  connectionstyle='arc3,rad=0.55', arrowstyle='-|>',
                                  mutation_scale=4.5, lw=0.8, color='#8C8C8C', clip_on=False))

# ---- C : D9 distribution ---------------------------------------------
letter(0.008, 0.672, 'C')
axC = FS.add_axes([0.080, 0.408, 0.205, 0.205])
nTox = int((s1v.D9 <= TOX_THR).sum()); nTot = int(s1v.D9.notna().sum())
axC.hist(s1v.D9.dropna(), bins=np.arange(0, 2.45, 0.05), color='#BFBFBF', ec='black', lw=0.4)
axC.axvline(TOX_THR, color='red', lw=1.1, ls='--')
axC.set_title(f'Toxic (D9 \u2264 {TOX_THR})\n{nTox}/{nTot} ({100*nTox/nTot:.1f}%)',
              color='red', fontsize=7.0, fontweight='bold', pad=2.5)
axC.text(0.97, 0.95, 'D9', transform=axC.transAxes, fontsize=8, va='top', ha='right')
axC.set_xlabel('Signal (Normalized)', fontsize=8, labelpad=1.5)
axC.set_ylabel('KZFPs Count', fontsize=8, labelpad=2)
axC.set_xlim(0, 2.4); axC.tick_params(labelsize=7.5, pad=1.5)
axC.spines[['top', 'right']].set_visible(False)

# ---- D : genomic distribution of peaks -------------------------------
letter(0.300, 0.672, 'D')
TE_C, TSS_C, OT_C = '#4F7F28', '#C45911', '#D9D9D9'
xb = np.arange(5)
axD1 = FS.add_axes([0.400, 0.505, 0.212, 0.108])
axD1.bar(xb, DD.TEs, 0.72, color=TE_C, ec='black', lw=0.4, label='TEs')
axD1.bar(xb, DD.TSS, 0.72, bottom=DD.TEs, color=TSS_C, ec='black', lw=0.4, label='TSS')
axD1.bar(xb, DD.Other, 0.72, bottom=DD.TEs + DD.TSS, color=OT_C, ec='black', lw=0.4, label='Other')
for i, r in DD.iterrows():
    if not r['mark']: continue
    yy = r.TEs / 2 if r.band == 'TEs' else r.TEs + r.TSS / 2
    axD1.text(i, yy, r['mark'], ha='center', va='center', fontsize=7, fontweight='bold')
axD1.set_ylim(0, 100); axD1.set_yticks([0, 25, 50, 75, 100])
axD1.set_ylabel('Genomic\nDistribution (%)', fontsize=8, labelpad=2)
axD1.set_xlim(-0.6, 4.6); axD1.set_xticks(xb); axD1.set_xticklabels([])
axD1.tick_params(labelsize=7.5, pad=1.5); axD1.spines[['top', 'right']].set_visible(False)
axD1.legend(frameon=False, fontsize=7, ncol=3, loc='lower left', bbox_to_anchor=(-0.03, 0.99),
            handlelength=1.0, handleheight=0.9, handletextpad=0.35, columnspacing=0.8)
axD2 = FS.add_axes([0.400, 0.408, 0.212, 0.072])
axD2.bar(xb, NTSS, 0.72, color='#595959', ec='black', lw=0.4)
for i, v in enumerate(NTSS):
    axD2.text(i, v + 130, f'{v:,.0f}', ha='center', va='bottom', fontsize=5.8)
axD2.set_ylim(0, 3450); axD2.set_yticks([0, 1500, 3000])
axD2.set_ylabel('TSS-proximal\ngenes (n)', fontsize=8, labelpad=2)
axD2.set_xlim(-0.6, 4.6); axD2.set_xticks(xb)
axD2.set_xticklabels(DD.label, rotation=45, ha='right', fontsize=6.6)
axD2.tick_params(labelsize=7.5, pad=1.5); axD2.spines[['top', 'right']].set_visible(False)

# ---- E : evolutionary age --------------------------------------------
letter(0.008, 0.335, 'E')
axE = FS.add_axes([0.080, 0.062, 0.455, 0.212])
bins = np.arange(0, AGE.max() + 10, 10)
axE.hist(AGE, bins=bins, orientation='horizontal', density=True, color='#BFBFBF', ec='black', lw=0.4)
axE.invert_yaxis(); axE.set_xlim(0, 0.050)
axE.set_ylabel('Evolutionary age (Myr)', fontsize=8, labelpad=2)
axE.set_xlabel('Fraction of KZFPs (density)', fontsize=8, labelpad=1.5)
axE.tick_params(labelsize=7.5, pad=1.5); axE.spines[['top', 'right']].set_visible(False)
for nm, a_, off in CAND_AGE:
    axE.axhline(a_, color='black', lw=0.6, ls=(0, (4, 3)))
    axE.text(0.0507, a_ + off, nm, fontsize=7.5, style='italic', fontweight='bold', va='center', ha='left')
axE.text(0.97, 0.05, f'n = {len(AGE)} KZFPs', transform=axE.transAxes, fontsize=7, ha='right')

# ---- F : per-KZFP TE-family enrichment -------------------------------
letter(0.628, 0.992, 'F')
FX0, FX1 = 0.700, 0.898
FS.legend(handles=[Patch(fc='white', ec='black', lw=.5, label='pAdj > 0.05'),
                   Patch(fc='#7F7F7F', ec='black', lw=.5, label='pAdj \u2264 0.05')],
          frameon=False, fontsize=7.5, loc='upper left', bbox_to_anchor=(0.665, 0.990),
          ncol=1, handlelength=1.1, handleheight=1.0, handletextpad=0.45, labelspacing=0.30,
          title='Target TE family\n(hypergeometric, BH-adj.)', title_fontsize=7.5, alignment='left')
axFc = FS.add_axes([FX0, 0.862, FX1 - FX0, 0.040])
hc.dendrogram(cl, ax=axFc, color_threshold=0, link_color_func=lambda k: 'black', no_labels=True)
axFc.set_axis_off()
axFr = FS.add_axes([0.655, 0.158, 0.043, 0.700])
hc.dendrogram(rl, ax=axFr, orientation='left', color_threshold=0,
              link_color_func=lambda k: 'black', no_labels=True)
axFr.set_axis_off()
axF = FS.add_axes([FX0, 0.158, FX1 - FX0, 0.700])
axF.imshow(Bo.values, aspect='auto', vmin=0, vmax=1,
           cmap=mpl.colors.ListedColormap(['white', '#7F7F7F']))
axF.set_xticks(range(Bo.shape[1]))
axF.set_xticklabels([c.replace('/', '.') for c in Bo.columns], rotation=45, ha='right', fontsize=7)
axF.set_yticks(range(len(Bo)))
axF.set_yticklabels(['' if n == 'ZNF202' else n for n in Bo.index], fontsize=6.4)
axF.yaxis.tick_right()
axF.set_xticks(np.arange(-.5, Bo.shape[1], 1), minor=True)
axF.set_yticks(np.arange(-.5, len(Bo), 1), minor=True)
axF.grid(which='minor', color='black', lw=0.25)
axF.tick_params(which='minor', length=0); axF.tick_params(length=1.5, pad=1.2)
for sp in axF.spines.values(): sp.set_linewidth(0.5)
FS.text(0.645, 0.508, f'Toxic KZFPs (n = {len(Bo)})', rotation=90, fontsize=8, va='center', ha='center')

FS.savefig(f'{OUT}/Figure_S1_PLOS.png', dpi=600)
FS.savefig(f'{OUT}/Figure_S1_PLOS.svg')
plt.close(FS)
print('written; toxic', nTox, '/', nTot, '| F rows', len(Bo), '| age n', len(AGE))
