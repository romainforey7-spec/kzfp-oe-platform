"""SUPERSEDED - recovery scratch, do not run.

Used once to recover the published S1 panels as images. Sets type down to 3.4 pt.
The delivered S1 Fig is built by build_figS1_fig.py.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image
from scipy import stats as st
import pypdfium2 as pdfium

R = r'C:\Users\forey\Desktop\20260924_For_Revision'
TX = os.path.join(R, '20250110_TE_Fam')

mpl.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'svg.fonttype': 'none'})

ALIAS = {'ZFP69': 'ZNF69', 'ZKSCAN7': 'ZNF167', 'ZKSCAN8': 'ZNF192'}

s1v = pd.read_excel(os.path.join(R, 'Supplementary Table S1.xlsx'))
s1v.columns = ['KZFP', 'D0', 'D4', 'D7', 'D9']
d9 = dict(zip(s1v.KZFP, s1v.D9))

def tox(n):
    v = d9.get(n, d9.get(ALIAS.get(n, n)))
    return None if v is None or pd.isna(v) else bool(v <= 0.85)

TSS = pd.read_csv('C:/Users/forey/.claude-science/orgs/fb5390a7-729b-4e9a-9b05-c86022c81dd9/artifacts/proj_467da4f3d043/f2ed15a8-efce-485b-8dfb-291570589bbe/v87fdd99e_kzfp_TSS_counts.csv')
pkk = set(TSS.KZNF)

gl = pd.read_csv(os.path.join(TX, 'GenomicLocation.txt'), sep='\t')
gl['tox'] = gl.KZFP.map(tox)
gl = gl.dropna(subset=['tox'])
gl['te'] = gl['TEs.pval'] <= 0.05
gl['pr'] = gl['promoter.pval'] <= 0.05
gl['cat'] = pd.Categorical(
    np.select([gl.te & ~gl.pr, ~gl.te & gl.pr, gl.te & gl.pr],
              ['TEs', 'Prom.', 'Both'], default='None'),
    ['TEs', 'Prom.', 'Both', 'None'], ordered=True)
gE = gl[gl.KZFP.isin(pkk)]
ctE2 = pd.crosstab(gE.cat, gE.tox.map({True: 'Toxic', False: 'Not Toxic'})).reindex(['TEs', 'Prom.', 'Both', 'None'])

def clean(fn):
    fn = re.sub(r'_(pubM|n)?_?vs_293T.*$', '', fn)
    for _ in range(2):
        fn = re.sub(r'_(rep\d?|lowconf|n_lowconf|n|pubM|rep_n|rep_pubM)$', '', fn)
    return fn

fam = pd.read_csv(os.path.join(TX, 'enrich_kzfp_perFam', 'Fam.tsv'), sep='\t')
fam['kzfp2'] = fam.filename.map(clean)
fam['tox2'] = fam.kzfp2.map(tox)
fs2 = fam[(fam['padj.hypergeom.reg'] <= 0.05) & fam.tox2.notna()]
FAMORD = ['LINE/L1', 'LTR/ERV1', 'LTR/ERVK', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'nonTE', 'Retroposon/SVA', 'SINE/Alu']
FF = pd.DataFrame([[f,
                    fs2[fs2.subfam_name == f][~fs2[fs2.subfam_name == f].tox2.astype(bool)].kzfp2.nunique(),
                    fs2[fs2.subfam_name == f][fs2[fs2.subfam_name == f].tox2.astype(bool)].kzfp2.nunique()]
                   for f in FAMORD], columns=['Family', 'NotToxic', 'Toxic'])

oe = os.path.join(TX, 'OE_ANALYSIS.xlsx')
mg = pd.read_excel(oe, sheet_name='merged_OE_Tables').drop_duplicates('Row.Labels')
agemap = mg.set_index('Row.Labels')['age_combined'].to_dict()
scanmap = mg.set_index('Row.Labels')['PresenceOfSCAN'].to_dict()
MK = s1v.copy()
MK['toxic'] = MK.D9 <= 0.85
MK['age'] = MK.KZFP.map(agemap)
MK['is_scan'] = MK.KZFP.map(scanmap).fillna('').str.contains('SCAN')

GD = TSS[['KZNF', 'n_TSS_genes']].merge(s1v, left_on='KZNF', right_on='KZFP', how='inner')
GD['logTSS'] = np.log10(GD.n_TSS_genes + 1)

Image.MAX_IMAGE_PIXELS = None
f1 = os.path.join(R, '20260727_Figs', 'Figure 1.TIF')
im1 = Image.open(f1).convert('RGB')
W1, H1 = im1.size
A = im1.crop((int(0.196 * W1), int(0.276 * H1), int(0.330 * W1), int(0.462 * H1)))

os.makedirs('out/figS1', exist_ok=True)

pdf = pdfium.PdfDocument(os.path.join(R, 'Supplementary Figures.pdf'))
pg = pdf[0].render(scale=6.0).to_pil()
PW, PH = pg.size

crops2 = {'A': (0.147, 0.196, 0.375, 0.412), 'B': (0.394, 0.196, 0.588, 0.398), 'D': (0.147, 0.447, 0.352, 0.648)}
for k, (x0, y0, x1, y1) in crops2.items():
    pg.crop((int(x0 * PW), int(y0 * PH), int(x1 * PW), int(y1 * PH))).save(f'out/figS1/FigS1{k}_carryover.png')

# Build FigS1C
NT_C, TX_C = "#7F7FFF", "#DF7F7F"

figC, axC = plt.subplots(figsize=(2.6, 2.1))
axC.hist(s1v.D9.dropna(), bins=np.arange(0, 2.45, 0.05), color='#BFBFBF', edgecolor='black', linewidth=0.4)
axC.axvline(0.85, color='red', lw=1.1, ls='--')
axC.text(0.5, 1.045, "Toxic (signal D9 < 0.85)", transform=axC.transAxes, color='red',
         fontsize=7.5, ha='center', fontweight='bold')
axC.text(0.06, 0.93, "D9", transform=axC.transAxes, fontsize=8, va='top')
axC.set_xlabel("Signal (Normalized)", fontsize=8)
axC.set_ylabel("KZFPs Count", fontsize=8)
axC.set_xlim(0, 2.4)
axC.tick_params(labelsize=7.5)
axC.spines[['top', 'right']].set_visible(False)
figC.savefig('out/figS1/FigS1C.png', dpi=600, bbox_inches='tight')
plt.close(figC)

# Build FigS1E
Ma = MK.dropna(subset=['age'])
figE2, axE2 = plt.subplots(figsize=(2.8, 2.2))
bins = np.arange(0, Ma.age.max() + 10, 10)
axE2.hist(Ma.age, bins=bins, orientation='horizontal', density=True, color='#BFBFBF',
          edgecolor='black', linewidth=0.4)
axE2.invert_yaxis()
axE2.set_xlim(0, 0.05)
axE2.set_ylabel("Age Combined", fontsize=8)
axE2.set_xlabel("Fraction of KZFPs (density)", fontsize=8)
axE2.tick_params(labelsize=7.5)
axE2.spines[['top', 'right']].set_visible(False)

CAND = {'ZNF257': None, 'ZNF43': None, 'ZNF498': None, 'ZNF18': None}
for c in CAND:
    CAND[c] = agemap.get(c, agemap.get({'ZNF498': 'ZSCAN25'}.get(c, c)))

for name, a_, off in [('ZNF257', 29.1, -9), ('ZNF43', 43.1, 9), ('ZNF498', 105.0, 0), ('ZNF18', 163.7, 0)]:
    axE2.text(0.0505, a_ + off, name, fontsize=7, style='italic', fontweight='bold', va='center', ha='left')

figE2.savefig('out/figS1/FigS1E.png', dpi=600, bbox_inches='tight')
plt.close(figE2)

# Build FigS1F
from scipy.cluster import hierarchy as hc
from scipy.spatial.distance import pdist

FAM8 = ['Retroposon/SVA', 'LTR/ERV1', 'LTR/ERVK', 'SINE/Alu', 'LTR/ERVL-MaLR', 'LTR/ERVL', 'nonTE', 'LINE/L1']
fw = fam[fam.tox2 == True].copy()
Mx = (fw.pivot_table(index='kzfp2', columns='subfam_name', values='padj.hypergeom.reg', aggfunc='min')
      .reindex(columns=FAM8))
Bin = (Mx <= 0.05).astype(int)
Bin = Bin[Bin.sum(axis=1) > 0]
rl = hc.linkage(pdist(Bin.values, metric='jaccard'), method='average')
cl = hc.linkage(pdist(Bin.values.T, metric='jaccard'), method='average')
ro = hc.leaves_list(rl)
co = hc.leaves_list(cl)
Bo = Bin.iloc[ro, co]

figF2 = plt.figure(figsize=(4.3, 6.0))
gsF = figF2.add_gridspec(2, 2, width_ratios=[0.13, 1], height_ratios=[0.09, 1],
                         hspace=0.012, wspace=0.012, left=0.04, right=0.80, top=0.925, bottom=0.155)
axd_c = figF2.add_subplot(gsF[0, 1])
hc.dendrogram(cl, ax=axd_c, color_threshold=0, above_threshold_color='black',
              link_color_func=lambda k: 'black', no_labels=True)
axd_c.set_axis_off()
axd_r = figF2.add_subplot(gsF[1, 0])
hc.dendrogram(rl, ax=axd_r, orientation='left', color_threshold=0,
              link_color_func=lambda k: 'black', no_labels=True)
axd_r.set_axis_off()
axh = figF2.add_subplot(gsF[1, 1])
axh.imshow(Bo.values, aspect='auto', cmap=mpl.colors.ListedColormap(['white', '#7F7F7F']), vmin=0, vmax=1)
axh.set_xticks(range(len(co)))
axh.set_xticklabels([Bo.columns[i].replace('/', '.') for i in range(len(co))],
                    rotation=45, ha='right', fontsize=7)
labels = [('' if n == 'ZNF202' else n) for n in Bo.index]
axh.set_yticks(range(len(Bo)))
axh.set_yticklabels(labels, fontsize=3.4)
axh.yaxis.tick_right()
axh.set_ylabel("")
axh.tick_params(length=1.5)
for s in axh.spines.values():
    s.set_linewidth(0.4)
axh.set_xticks(np.arange(-.5, len(co), 1), minor=True)
axh.set_yticks(np.arange(-.5, len(Bo), 1), minor=True)
axh.grid(which='minor', color='black', lw=0.25)
axh.tick_params(which='minor', length=0)
figF2.text(0.015, 0.55, "Toxic KZFPs", rotation=90, fontsize=8, va='center')
lg = [mpl.patches.Patch(fc='white', ec='black', lw=.5, label='pAdj > 0.05'),
      mpl.patches.Patch(fc='#7F7F7F', ec='black', lw=.5, label='pAdj < 0.05')]
figF2.legend(handles=lg, frameon=False, fontsize=6.8, loc='upper left', bbox_to_anchor=(0.02, 1.005),
             ncol=2, handlelength=1.0, columnspacing=0.9,
             title="Target (Enrichment: Fisher exact test)", title_fontsize=6.8)
figF2.savefig('out/figS1/FigS1F.png', dpi=600)
plt.close(figF2)

# Assemble Figure S1
plt.close('all')
FS = plt.figure(figsize=(7.48, 8.2))
gsS = FS.add_gridspec(3, 12, height_ratios=[1.00, 0.66, 0.60], hspace=0.06, wspace=0.06,
                      left=0.015, right=0.985, top=0.985, bottom=0.012)

def imgpanel(gsslot, path, letter, dx=-0.01, dy=1.02):
    ax = FS.add_subplot(gsslot)
    ax.imshow(np.asarray(Image.open(path)))
    ax.axis('off')
    ax.text(dx, dy, letter, transform=ax.transAxes, fontweight='bold', fontsize=10, va='top', ha='left')

imgpanel(gsS[0, 0:5], 'out/figS1/FigS1A_carryover.png', 'A')
imgpanel(gsS[0, 5:9], 'out/figS1/FigS1B_carryover.png', 'B')
imgpanel(gsS[0, 9:12], 'out/figS1/FigS1C.png', 'C')
imgpanel(gsS[1, 0:4], 'out/figS1/FigS1D_carryover.png', 'D')
imgpanel(gsS[2, 0:5], 'out/figS1/FigS1E.png', 'E')
imgpanel(gsS[1:3, 4:12], 'out/figS1/FigS1F.png', 'F')

FS.savefig('out/figS1/Figure_S1_PLOS.png', dpi=600)