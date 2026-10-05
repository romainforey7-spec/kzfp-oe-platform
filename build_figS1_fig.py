"""S1 Fig assembly for PLOS Genetics (7.48 x 8.70 in, all text >= 8 pt).

This is the builder of the delivered S1 Fig. Its 264 text elements are all at 8 pt
or above - do not confuse it with build_figS1.py, an earlier draft that sets type as
small as 3.2 pt and writes to the same filenames.

NOT a standalone module. Like build_fig1.py it is the drawing half only, and expects
its data objects in the executing namespace - among them WELL (the per-well plate
table), Bo (the toxic-KZFP by TE-family boolean matrix), NTSS (genes with a TSS peak
per bait), AGE and CAND_AGE (evolutionary ages), SAMPS, DD, oKG, s1v, hc, imA (the
panel A photograph) and TOX_THR. Run it as:

    ns = figureS1_inputs()         # the data assembly step, see REPO_GAPS.md
    exec(open("build_figS1_fig.py").read(), ns)

Importing this file on its own raises NameError, by design rather than by accident.
"""

from plos_export import PK  # pins output to the figure's declared size
plt.close('all')
mpl.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':0.8,
                     'xtick.major.width':0.8,'ytick.major.width':0.8,
                     'xtick.major.size':2.4,'ytick.major.size':2.4,'svg.fonttype':'none'})
FW, FH = 7.48, 8.70
FS = plt.figure(figsize=(FW, FH))
TAB_CMAP = LinearSegmentedColormap.from_list('gyr', ['#63BE7B', '#FFEB84', '#F8696B'])
GREY, BLACK = '#A6A6A6', '#1A1A1A'
def LT(x, y, s): FS.text(x, y, s, fontweight='bold', fontsize=11, va='top', ha='left')

# ---- A : plates ------------------------------------------------------
axA = FS.add_axes([0.020, 0.600, 0.232, 0.385]); axA.imshow(np.asarray(imA)); axA.axis('off')
LT(0.004, 0.992, 'A')

# ---- C : D9 signal distribution --------------------------------------
LT(0.004, 0.556, 'C')
axC = FS.add_axes([0.085, 0.332, 0.195, 0.192])
nTox = int((s1v.D9 <= TOX_THR).sum()); nTot = int(s1v.D9.notna().sum())
axC.hist(s1v.D9.dropna(), bins=np.arange(0, 2.55, 0.05), color='#BFBFBF', ec='black', lw=0.4)
axC.axvline(TOX_THR, color='red', lw=1.1, ls='--')
axC.set_title(f'Toxic (D9 \u2264 {TOX_THR}): {nTox}/{nTot}', color='red', fontsize=8, fontweight='bold', pad=3)
axC.text(0.96, 0.92, 'D9', transform=axC.transAxes, fontsize=8, va='top', ha='right')
axC.set_xlabel('Signal (Normalized)', fontsize=8, labelpad=2)
axC.set_ylabel('KZFPs count', fontsize=8, labelpad=2)
axC.set_xlim(0, 2.5); axC.tick_params(labelsize=8, pad=2)
axC.spines[['top','right']].set_visible(False)

# ---- E : evolutionary age -------------------------------------------
LT(0.645, 0.556, 'E')
axE = FS.add_axes([0.735, 0.332, 0.205, 0.192])
bins = np.arange(0, AGE.max() + 10, 10)
axE.hist(AGE, bins=bins, orientation='horizontal', density=True, color='#BFBFBF', ec='black', lw=0.4)
axE.invert_yaxis(); axE.set_xlim(0, 0.085); axE.set_xticks([0, 0.04, 0.08])
axE.set_ylabel('Evol. age (Myr)', fontsize=8, labelpad=2)
axE.set_xlabel('Fraction of KZFPs', fontsize=8, labelpad=2)
axE.tick_params(labelsize=8, pad=2); axE.spines[['top','right']].set_visible(False)
for nm, a_, off in CAND_AGE:
    axE.plot([0, 0.049], [a_, a_], color='black', lw=0.6, ls=(0, (3.5, 2.5)))
    axE.text(0.0515, a_ + off, nm, fontsize=8, style='italic', fontweight='bold', va='center', ha='left')

# ---- B : PrestoBlue quantification ----------------------------------
LT(0.258, 0.992, 'B')
TX0, TX1 = 0.305, 0.950
axBn = FS.add_axes([TX0, 0.921, TX1 - TX0, 0.001]); axBn.set_axis_off()
axBt = FS.add_axes([TX0, 0.812, TX1 - TX0, 0.110])
axBt.imshow(WELL, aspect='auto', cmap=TAB_CMAP, vmin=0, vmax=np.nanmax(WELL))
for i in range(8):
    for j in range(12):
        v = WELL[i, j]
        axBt.text(j, i, f'{0.0 if abs(v) < 0.05 else v:.1f}', ha='center', va='center', fontsize=8)
axBt.set_xticks(range(12)); axBt.set_xticklabels([''] + SAMPS + [''], rotation=90, fontsize=8)
axBt.xaxis.set_ticks_position('top')
axBt.set_yticks(range(8)); axBt.set_yticklabels(list('ABCDEFGH'), fontsize=8)
axBt.set_xticks(np.arange(-.5, 12, 1), minor=True); axBt.set_yticks(np.arange(-.5, 8, 1), minor=True)
axBt.grid(which='minor', color='white', lw=0.5); axBt.tick_params(which='minor', length=0)
axBt.tick_params(length=1.2, pad=1.5)
for sp in axBt.spines.values(): sp.set_linewidth(0.5)
for r0, r1, lb in [(1, 3, '\u2212Dox'), (4, 6, '+Dox')]:
    FS.text(0.279, 0.812 + 0.110 * (1 - (r0 + r1 + 1) / 16), lb, fontsize=8, va='center',
            ha='center', rotation=90)
X = np.arange(10); w = 0.38
ekw = dict(lw=0.6, capsize=1.2, capthick=0.6)
axB1 = FS.add_axes([TX0, 0.712, TX1 - TX0, 0.066])
axB1.bar(X - w/2, nd_m, w, yerr=nd_s, color=GREY, ec='black', lw=0.5, error_kw=ekw)
axB1.bar(X + w/2, dx_m, w, yerr=dx_s, color=BLACK, ec='black', lw=0.5, error_kw=ekw)
axB1.axhline(1, color='black', lw=0.5, ls=':')
axB1.set_ylim(0, 2); axB1.set_yticks([0, 1, 2])
axB1.set_ylabel('A570/A600', fontsize=8, labelpad=2)
axB1.set_xlim(-0.65, 9.65); axB1.set_xticks(X); axB1.set_xticklabels([])
axB1.tick_params(labelsize=8, pad=2); axB1.spines[['top','right']].set_visible(False)
axB1.legend(handles=[Patch(fc=GREY, ec='black', lw=.5, label='\u2212Dox'),
                     Patch(fc=BLACK, ec='black', lw=.5, label='+Dox')],
            frameon=False, fontsize=8, ncol=2, loc='lower left', bbox_to_anchor=(-0.005, 1.03),
            handlelength=1.1, handleheight=1.0, handletextpad=0.4, columnspacing=1.6)
axB2 = FS.add_axes([TX0, 0.630, TX1 - TX0, 0.066])
axB2.bar(X - w/2, nn_m, w, yerr=nn_s, color=GREY, ec='black', lw=0.5, error_kw=ekw)
axB2.bar(X + w/2, nx_m, w, yerr=nx_s, color=BLACK, ec='black', lw=0.5, error_kw=ekw)
axB2.axhline(1, color='black', lw=0.5, ls=':')
axB2.axhline(TOX_THR, color='red', lw=0.8, ls='--')
axB2.set_ylim(0, 2); axB2.set_yticks([0, 1, 2])
axB2.set_ylabel('Signal\n(Normalized)', fontsize=8, labelpad=2)
axB2.set_xlim(-0.65, 9.65); axB2.set_xticks(X)
axB2.set_xticklabels(SAMPS, rotation=45, ha='right', fontsize=8)
axB2.set_xlabel('KZFPs (1 plate \u2013 10 candidates)', fontsize=8, labelpad=2)
axB2.tick_params(labelsize=8, pad=2); axB2.spines[['top','right']].set_visible(False)
FS.patches.append(FancyArrowPatch((0.967, 0.718), (0.967, 0.674), transform=FS.transFigure,
                                  arrowstyle='-|>', mutation_scale=7, lw=1.3,
                                  color='#8C8C8C', clip_on=False, shrinkA=0, shrinkB=0))

# ---- D : genomic distribution of peaks ------------------------------
LT(0.296, 0.556, 'D')
TE_C, TSS_C, OT_C = '#4F7F28', '#C45911', '#D9D9D9'
xb = np.arange(5)
DLAB = ['ZNF257', f'Other KZFPs\n(n = {len(oKG)})', 'ZNF498', 'ZNF18', 'ZNF43']
axD1 = FS.add_axes([0.400, 0.412, 0.215, 0.112])
axD1.bar(xb, DD.TEs, 0.72, color=TE_C, ec='black', lw=0.5, label='TEs')
axD1.bar(xb, DD.TSS, 0.72, bottom=DD.TEs, color=TSS_C, ec='black', lw=0.5, label='TSS')
axD1.bar(xb, DD.Other, 0.72, bottom=DD.TEs + DD.TSS, color=OT_C, ec='black', lw=0.5, label='Other')
for i, r in DD.iterrows():
    if not r['mark']: continue
    yy = r.TEs / 2 if r.band == 'TEs' else r.TEs + r.TSS / 2
    axD1.text(i, yy, r['mark'], ha='center', va='center', fontsize=8, fontweight='bold')
axD1.set_ylim(0, 100); axD1.set_yticks([0, 50, 100])
axD1.set_ylabel('Genomic\ndistribution (%)', fontsize=8, labelpad=2)
axD1.set_xlim(-0.6, 4.6); axD1.set_xticks(xb); axD1.set_xticklabels([])
axD1.tick_params(labelsize=8, pad=2); axD1.spines[['top','right']].set_visible(False)
axD1.legend(frameon=False, fontsize=8, ncol=3, loc='lower left', bbox_to_anchor=(-0.04, 0.98),
            handlelength=1.0, handleheight=1.0, handletextpad=0.35, columnspacing=0.9)
axD2 = FS.add_axes([0.400, 0.332, 0.215, 0.062])
axD2.bar(xb, NTSS, 0.72, color='#595959', ec='black', lw=0.5)
for i, v in enumerate(NTSS):
    axD2.text(i, v + 260, f'{v:,.0f}', ha='center', va='bottom', fontsize=8)
axD2.set_ylim(0, 7000); axD2.set_yticks([0, 3000, 6000])
axD2.set_ylabel('Genes with\na TSS peak', fontsize=8, labelpad=2)
axD2.set_xlim(-0.6, 4.6); axD2.set_xticks(xb)
axD2.set_xticklabels(DLAB, rotation=45, ha='right', fontsize=8)
axD2.tick_params(labelsize=8, pad=2); axD2.spines[['top','right']].set_visible(False)

# ---- F : TE-family enrichment across toxic KZFPs (transposed) --------
LT(0.004, 0.268, 'F')
BoT = Bo.T                                      # families x KZFPs
FS.legend(handles=[Patch(fc='white', ec='black', lw=.5, label='pAdj > 0.05'),
                   Patch(fc='#7F7F7F', ec='black', lw=.5, label='pAdj \u2264 0.05')],
          frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(0.430, 0.056),
          ncol=2, handlelength=1.0, handleheight=1.0, handletextpad=0.4, columnspacing=1.2,
          title='Target TE family (hypergeometric, BH-adj.)', title_fontsize=8)
axFc = FS.add_axes([0.150, 0.236, 0.835, 0.024])
hc.dendrogram(rl, ax=axFc, color_threshold=0, link_color_func=lambda k: 'black', no_labels=True)
axFc.set_axis_off()
axF = FS.add_axes([0.150, 0.120, 0.835, 0.112])
axF.imshow(BoT.values, aspect='auto', vmin=0, vmax=1,
           cmap=mpl.colors.ListedColormap(['white', '#7F7F7F']))
axF.set_yticks(range(BoT.shape[0]))
axF.set_yticklabels([c.replace('/', '.') for c in BoT.index], fontsize=8)
axF.set_xticks(range(BoT.shape[1]))
axF.set_xticklabels(['' if n == 'ZNF202' else n for n in BoT.columns], rotation=90, fontsize=8)
axF.set_xticks(np.arange(-.5, BoT.shape[1], 1), minor=True)
axF.set_yticks(np.arange(-.5, BoT.shape[0], 1), minor=True)
axF.grid(which='minor', color='black', lw=0.25)
axF.tick_params(which='minor', length=0); axF.tick_params(length=1.5, pad=1.5)
for sp in axF.spines.values(): sp.set_linewidth(0.5)
FS.text(0.160, 0.040, f'Toxic KZFPs (n = {BoT.shape[1]})', fontsize=8, ha='left')

for ax in FS.axes:
    for col in ax.collections: col.set_rasterized(True)
print("sub-8pt:", len(font_audit(FS)), font_audit(FS)[:6])
FS.savefig(f'{OUT}/Figure_S1_PLOS.png', dpi=600, **PK(FS))
FS.savefig(f'{OUT}/Figure_S1_300.png', dpi=300, **PK(FS)); FS.savefig(f'{OUT}/Figure_S1_PLOS.svg', **PK(FS))
plt.close(FS)

# --- submission TIFF: written from the NATIVE 300-dpi render, never resampled
from PIL import Image as _PILI
_PILI.MAX_IMAGE_PIXELS = None
_im = _PILI.open(f"{OUT}/Figure_S1_300.png").convert("RGB")
_im.save(f"{OUT}/Figure_S1_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
print("tif", _im.size)
