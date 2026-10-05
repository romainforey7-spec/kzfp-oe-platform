"""Canonical figure export for the PLOS Genetics submission.

Every submission sheet must come out at an exact physical size: PLOS scales
artwork by width, so a sheet that is even slightly wider than declared has its
type scaled down below the 8 pt floor.

Two matplotlib behaviours make that easy to get wrong, and both are handled here:

1. ``savefig.bbox = "tight"`` -- if this rcParam is set anywhere in the session
   (several style sheets and notebook front-ends set it), matplotlib ignores the
   figure's declared size and crops to the drawn content, silently changing the
   sheet dimensions. We clear it at import and pass an explicit bounding box on
   every write.
2. ``pad_inches`` defaults to 0.1 in, adding 0.2 in to both dimensions. We pass
   0.0 explicitly.

Use ``export_sheet(F, "Figure_3")`` and nothing else. It writes, into ``outdir``:
``<stem>_PLOS.png`` at 600 dpi, ``<stem>_300.png`` at 300 dpi,
``<stem>_PLOS.svg`` with live text, and ``<stem>_PLOS.tif`` (LZW, RGB, 300 dpi)
regenerated from the 300 dpi raster so the TIFF and the PNG cannot disagree.
"""
import os
import matplotlib as mpl
from matplotlib.transforms import Bbox
from PIL import Image

mpl.rcParams["savefig.bbox"] = None
mpl.rcParams["savefig.pad_inches"] = 0.0


def PK(F):
    """savefig keyword arguments that pin the output to the figure's own size."""
    w, h = F.get_size_inches()
    return dict(bbox_inches=Bbox([[0.0, 0.0], [w, h]]), pad_inches=0.0)


def export_sheet(F, stem, outdir=".", dpi_hi=600, dpi_lo=300):
    os.makedirs(outdir, exist_ok=True)
    kw = PK(F)
    p_hi = os.path.join(outdir, f"{stem}_PLOS.png")
    p_lo = os.path.join(outdir, f"{stem}_300.png")
    p_sv = os.path.join(outdir, f"{stem}_PLOS.svg")
    p_tf = os.path.join(outdir, f"{stem}_PLOS.tif")
    F.savefig(p_hi, dpi=dpi_hi, **kw)
    F.savefig(p_lo, dpi=dpi_lo, **kw)
    F.savefig(p_sv, **kw)
    Image.open(p_lo).convert("RGB").save(
        p_tf, format="TIFF", compression="tiff_lzw", dpi=(dpi_lo, dpi_lo))
    return dict(png=p_hi, png300=p_lo, svg=p_sv, tif=p_tf)


PLOS_BOUNDS = dict(w_px_min=789, w_px_max=2250, h_px_max=2625,
                   dpi_min=300, w_in_min=2.63, w_in_max=7.5,
                   h_in_max=8.75, mb_max=10.0)


def check_sheet(path):
    """Measure one submission file against the PLOS bounds. Returns a dict."""
    im = Image.open(path)
    dpi = im.info.get("dpi", (None, None))[0]
    dpi = float(dpi) if dpi else None
    w, h = im.size
    mb = os.path.getsize(path) / 1e6
    r = dict(file=os.path.basename(path), px_w=w, px_h=h, dpi=dpi, mode=im.mode,
             compression=im.info.get("compression"), mb=round(mb, 2))
    if dpi:
        r["in_w"], r["in_h"] = round(w / dpi, 2), round(h / dpi, 2)
    B = PLOS_BOUNDS
    fail = []
    if not (B["w_px_min"] <= w <= B["w_px_max"]): fail.append("width_px")
    if h > B["h_px_max"]: fail.append("height_px")
    if not dpi or dpi < B["dpi_min"]: fail.append("dpi")
    if dpi and not (B["w_in_min"] <= w / dpi <= B["w_in_max"]): fail.append("width_in")
    if dpi and h / dpi > B["h_in_max"]: fail.append("height_in")
    if mb > B["mb_max"]: fail.append("size_mb")
    r["result"] = "PASS" if not fail else "FAIL: " + ",".join(fail)
    return r
