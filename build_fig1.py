"""Figure 1 assembly for PLOS Genetics (7.48 in wide, all text >= 8 pt).

NOT a standalone module. This file is the drawing half of Figure 1 and expects its
data objects to exist in the namespace it is executed in - 23 names, among them the
screen table, the age map, the TE-family matrix and the toxicity threshold TOX_THR.
Run it as:

    ns = figure1_inputs()          # the data assembly step, see REPO_GAPS.md
    exec(open("build_fig1.py").read(), ns)

The data-assembly step is the one piece of the figure pipeline not yet in this
repository; the delivered Figure 1 was built with those objects live in the session.
Importing this file on its own raises NameError, by design rather than by accident.
"""

plt.close('all')
from matplotlib.colors import TwoSlopeNorm
from plos_export import PK  # pins output to the figure's declared size
NT_C, TX_C = '#7F7FFF', '#DF7F7F'
F1W, F1H = 7.48, 8.70
G1F = plt.figure(figsize=(F1W, F1H))
def L1(ax, s, dx=-0.235, dy=1.10):
    ax.text(dx, dy, s, transform=ax.transAxes, fontweight='bold', fontsize=11.0, va='bottom', ha='left')

# ---------- A : workflow schematic -----------------------------------
axA = G1F.add_axes([0.030, 0.700, 0.200, 0.272]); axA.axis('off')
axA.set_xlim(0, 1); axA.set_ylim(0, 1)
STEPS = ['Plasmid:\npTRE-KZFP-HA\n(Dox inducible)', 'Lentivector', 'Transduction K562',
         'Puro selection', 'Stable cell line', 'Dox induction',
         'D4 \u2192 D7 \u2192 D9\nPrestoBlue\n(metabolic activity,\nproliferation)']
NL = [t.count('\n') + 1 for t in STEPS]
HH = [0.0235 * n + 0.002 for n in NL]               # half-height of each text block
GAP = 0.034
tot = 2 * sum(HH) + GAP * (len(STEPS) - 1)
top = 0.955
ys = []
cur = top
for k, h in enumerate(HH):
    cur -= h
    ys.append(cur)
    cur -= h + GAP
for k, (txt, y) in enumerate(zip(STEPS, ys)):
    axA.text(0.5, y, txt, ha='center', va='center', fontsize=8.0, linespacing=1.25)
    if k < len(STEPS) - 1:
        axA.plot([0.5, 0.5], [y - HH[k] - 0.003, ys[k+1] + HH[k+1] + 0.003],
                 lw=0.7, color='black', solid_capstyle='butt')
yl = ys[-1] - HH[-1] - 0.006
axA.annotate('', xy=(0.5, max(yl - 0.055, 0.005)), xytext=(0.5, yl),
             arrowprops=dict(arrowstyle='-|>', lw=0.9, color='black', mutation_scale=8,
                             shrinkA=0, shrinkB=0))
L1(axA, 'A', dx=-0.14, dy=1.00)

# ---------- B : per-KZFP signal heatmap ------------------------------
HM = s1v.dropna(subset=['D4', 'D7', 'D9']).sort_values('D9', ascending=False)
axB = G1F.add_axes([0.310, 0.718, 0.120, 0.235])
CMB = LinearSegmentedColormap.from_list('rwb', ['#B40000', '#FF4040', '#FFFFFF', '#6E6EDB', '#00008B'])
NRM = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=2.5)
imh = axB.imshow(HM[['D4', 'D7', 'D9']].values, aspect='auto', cmap=CMB, norm=NRM,
                 interpolation='nearest')
axB.set_xticks([0, 1, 2]); axB.set_xticklabels(['4', '7', '9'], fontsize=8.0)
axB.set_yticks([]); axB.set_ylabel(f'KZFPs (n = {len(HM)})', fontsize=8.0, labelpad=3)
axB.set_xlabel('Days post induction', fontsize=8.0, labelpad=2)
for sp in axB.spines.values(): sp.set_linewidth(0.7)
L1(axB, 'B', dx=-0.36, dy=1.03)

cax = G1F.add_axes([0.452, 0.762, 0.015, 0.160])
cb = G1F.colorbar(imh, cax=cax, ticks=[0, 1, 2])
cb.ax.yaxis.set_ticks_position('left')
cb.ax.tick_params(labelsize=8.0, length=2, pad=1.5)
cb.ax.set_title('Norm.\nA570/A600', fontsize=8.0, pad=4, loc='left')
cb.outline.set_linewidth(0.6)
GFP9 = float(CS.loc[(CS.Timepoint == 'D9') & (CS.Control == 'GFP'), 'mean_ratio'].iloc[0])
LACZ9 = float(CS.loc[(CS.Timepoint == 'D9') & (CS.Control == 'LacZ'), 'mean_ratio'].iloc[0])
for val, lab, col, dyp in [(GFP9, f'GFP {GFP9:.2f}', '#1A1A1A', 2.0),
                           (TOX_THR, f'cut-off {TOX_THR}', 'red', 0.0),
                           (LACZ9, f'LacZ {LACZ9:.2f}', '#1A1A1A', -7.0)]:
    yf = NRM(val)
    cb.ax.plot([0, 1], [yf, yf], transform=cb.ax.transAxes, color=col, lw=1.0,
               ls='--' if col == 'red' else '-', clip_on=False)
    cb.ax.annotate(lab, xy=(1.0, yf), xycoords='axes fraction', xytext=(4.6, dyp),
                   textcoords='offset points', fontsize=8.0, color=col, va='center', ha='left',
                   annotation_clip=False)

# ---------- C : toxic fraction ---------------------------------------
nTox = int((s1v.D9 <= TOX_THR).sum()); nTot = int(s1v.D9.notna().sum()); nNT = nTot - nTox
axC = G1F.add_axes([0.628, 0.715, 0.170, 0.240])
axC.pie([nNT, nTox], colors=[NT_C, TX_C], startangle=90, counterclock=False,
        wedgeprops=dict(ec='black', lw=0.6))
axC.text(1.12, 0.42, f'Toxic: {nTox} ({100*nTox/nTot:.1f}%)', fontsize=8.0, va='center', ha='left')
axC.text(1.12, -0.52, f'Not toxic: {nNT} ({100*nNT/nTot:.1f}%)', fontsize=8.0, va='center', ha='left')
axC.set_xlim(-1.1, 1.1)
L1(axC, 'C', dx=-0.14, dy=0.95)

# ---------- D : evolutionary age by toxicity --------------------------
MK = s1v.copy(); MK['toxic'] = MK.D9 <= TOX_THR
MK['age'] = MK.KZFP.map(agemap); Ma = MK.dropna(subset=['age'])
axD = G1F.add_axes([0.075, 0.400, 0.215, 0.235])
ab = np.arange(0, Ma.age.max() + 10, 10)
for sel, c, lb in [(~Ma.toxic, NT_C, 'Not toxic'), (Ma.toxic, TX_C, 'Toxic')]:
    axD.hist(Ma.age[sel], bins=ab, orientation='horizontal', density=True,
             color=c, ec='black', lw=0.3, alpha=0.62, label=f'{lb} (n = {int(sel.sum())})')
axD.invert_yaxis(); axD.set_xlim(0, 0.050)
axD.set_ylabel('Evolutionary age (Myr)', fontsize=8.0, labelpad=2)
axD.set_xlabel('Fraction of KZFPs (density)', fontsize=8.0, labelpad=2)
axD.tick_params(labelsize=8.0, pad=1.5); axD.spines[['top', 'right']].set_visible(False)
axD.legend(frameon=False, fontsize=8.0, loc='lower right', handlelength=0.9, handleheight=0.9,
           handletextpad=0.35, labelspacing=0.25)
L1(axD, 'D', dx=-0.235, dy=1.02)

# ---------- E / F : stacked composition ------------------------------
def stacked(ax, labels, nt, tx, xlab, letter, dx, rot=45, fs=8.0):
    nt = np.asarray(nt, float); tx = np.asarray(tx, float); tot = nt + tx
    pt = 100 * tx / tot; pn = 100 * nt / tot
    xi = np.arange(len(labels))
    ax.bar(xi, pt, 0.70, color=TX_C, ec='black', lw=0.4, label='Toxic')
    ax.bar(xi, pn, 0.70, bottom=pt, color=NT_C, ec='black', lw=0.4, label='Not toxic')
    for i in range(len(labels)):
        ax.text(i, pt[i] / 2, f'{int(tx[i])}', ha='center', va='center', fontsize=fs)
        ax.text(i, pt[i] + pn[i] / 2, f'{int(nt[i])}', ha='center', va='center', fontsize=fs)
    ax.set_ylim(0, 100); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.75, len(labels) - 0.25)
    ax.set_xticks(xi); ax.set_xticklabels(labels, rotation=rot, ha='right', fontsize=8.0)
    ax.set_ylabel('KZFPs (%)', fontsize=8.0, labelpad=2)
    ax.tick_params(labelsize=8.0, pad=1.5); ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8.0, ncol=2, loc='lower left', bbox_to_anchor=(-0.02, 1.00),
              handlelength=0.9, handleheight=0.9, handletextpad=0.35, columnspacing=0.8)
    L1(ax, letter, dx=dx, dy=1.10)

axE = G1F.add_axes([0.408, 0.400, 0.135, 0.235])
stacked(axE, list(ct2.index), ct2['Not Toxic'].values, ct2['Toxic'].values, None, 'E', -0.40)
axF = G1F.add_axes([0.660, 0.400, 0.320, 0.235])
stacked(axF, [f.replace('/', '.') for f in FF.Family], FF.NotToxic.values, FF.Toxic.values, None, 'F', -0.175)

# ---------- G : signal vs number of bound TSS ------------------------
GD = G2.copy(); GD['lx'] = np.log10(GD.n_TSS_genes + 1)
for k, tp in enumerate(['D4', 'D7', 'D9']):
    ax = G1F.add_axes([0.078 + k * 0.212, 0.068, 0.158, 0.215])
    ax.scatter(GD.lx, GD[tp], s=3.0, color='#333333', alpha=0.55, lw=0, rasterized=True)
    sl, ic, r, p, se = st.linregress(GD.lx, GD[tp])
    xx = np.linspace(GD.lx.min(), GD.lx.max(), 60)
    resid = GD[tp] - (sl * GD.lx + ic); sd = resid.std(ddof=2)
    ax.plot(xx, sl * xx + ic, color='#2E6DA4', lw=1.1)
    ax.fill_between(xx, sl * xx + ic - 1.96 * sd / np.sqrt(len(GD)) * 2.6,
                    sl * xx + ic + 1.96 * sd / np.sqrt(len(GD)) * 2.6,
                    color='#2E6DA4', alpha=0.22, lw=0)
    ax.set_xticks([0, 2, 4]); ax.set_xticklabels(['1', '100', '10,000'], fontsize=8.0)
    ax.set_ylim(-0.1, 2.6); ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax.set_ylabel(f'Signal {tp}', fontsize=8.0, labelpad=2)
    ax.set_xlabel('Genes with a TSS peak', fontsize=8.0, labelpad=1.5)
    ax.text(0.03, 0.97, f'R = {r:.2f}\np = {p:.3f}', transform=ax.transAxes,
            fontsize=8.0, va='top', ha='left')
    ax.tick_params(labelsize=8.0, pad=1.5); ax.spines[['top', 'right']].set_visible(False)
    if k == 0: L1(ax, 'G', dx=-0.305, dy=1.03)
G1F.text(0.290, 0.013, f'log10 scale (+1), n = {len(GD)} KZFPs with ChIP data',
         fontsize=8.0, ha='center')

# ---------- H : SCAN vs no SCAN ---------------------------------------
MK['is_scan'] = MK.KZFP.map(scanmap).fillna('').astype(str).str.contains('SCAN')
sc = MK.loc[MK.is_scan, 'D9'].dropna(); ns = MK.loc[~MK.is_scan, 'D9'].dropna()
U, pw = st.mannwhitneyu(sc, ns, alternative='less')
axH = G1F.add_axes([0.790, 0.068, 0.175, 0.215])
vp = axH.violinplot([sc.values, ns.values], positions=[0, 1], widths=0.78, showextrema=False)
for b in vp['bodies']: b.set_facecolor('#A6A6A6'); b.set_edgecolor('black'); b.set_lw(0.6); b.set_alpha(1.0)
axH.boxplot([sc.values, ns.values], positions=[0, 1], widths=0.16, showfliers=False,
            boxprops=dict(lw=0.6, facecolor='white'), medianprops=dict(lw=0.8, color='black'),
            whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), patch_artist=True)
axH.set_xticks([0, 1]); axH.set_xticklabels([f'SCAN\n(n = {len(sc)})', f'No SCAN\n(n = {len(ns)})'], fontsize=8.0)
axH.set_ylabel('Norm. A570/A600 (D9)', fontsize=8.0, labelpad=2)
axH.set_ylim(-0.1, 3.0); axH.tick_params(labelsize=8.0, pad=1.5)
axH.spines[['top', 'right']].set_visible(False)
yb = 2.62
axH.plot([0, 0, 1, 1], [yb, yb + 0.09, yb + 0.09, yb], lw=0.7, color='black')
axH.text(0.5, yb + 0.13, f'p = {pw:.3f}', ha='center', va='bottom', fontsize=8.0)
L1(axH, 'H', dx=-0.30, dy=1.03)

for ax in G1F.axes:
    for col in ax.collections:
        if isinstance(col, mpl.collections.PathCollection): col.set_rasterized(True)
os.makedirs('out/fig1', exist_ok=True)
G1F.savefig('out/fig1/Figure_1_PLOS.png', dpi=600, **PK(F))
G1F.savefig('out/fig1/Figure_1_PLOS.svg', **PK(F))
plt.close(G1F)
print(f'Fig1: toxic {nTox}/{nTot} | heatmap rows {len(HM)} | G n={len(GD)} | SCAN {len(sc)} vs {len(ns)} p={pw:.4f}')
print(f'GFP D9 {GFP9:.3f}  LacZ D9 {LACZ9:.3f}')

# --- submission TIFF: written from the NATIVE 300-dpi render, never resampled
from PIL import Image as _PILI
_PILI.MAX_IMAGE_PIXELS = None
_im = _PILI.open("out/Figure_1_300.png").convert("RGB")
_im.save("out/Figure_1_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
print("tif", _im.size)
