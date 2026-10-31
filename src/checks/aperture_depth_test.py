"""Depth of two phase features (e.g. primary eclipse and a second dip) vs aperture radius, per sector.

A feature from the target keeps a constant depth ratio as the aperture grows; a feature from a
neighbour changes its ratio to the target's eclipse (it grows if the neighbour is outside the small
aperture). Simple photometry on the TESScut cutouts: circular apertures at the target, sky = median
of pixels > 4 px from the target, per-orbit quadratic detrending with both features masked.

Usage: aperture_depth_test.py CUTOUT_GLOB RA DEC P T0 PH1 HW1 PH2 HW2
"""
import glob, re, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
import lightkurve as lk
from astropy.coordinates import SkyCoord

pat, ra, dec, P, T0, ph1, hw1, ph2, hw2 = sys.argv[1], *map(float, sys.argv[2:10])
c = SkyCoord(ra, dec, unit="deg")
print("sector  r(px)  depth1(%)  depth2(ppt)  ratio2/1")
for fn in sorted(glob.glob(pat), key=lambda x: int(re.search(r"s0*(\d+)-", x).group(1))):
    s = int(re.search(r"s0*(\d+)-", fn).group(1))
    pf = lk.read(fn); F = pf.flux.value; t = pf.time.value
    ok = np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1) & (pf.quality == 0)
    F, t = F[ok], t[ok]
    x0, y0 = pf.wcs.world_to_pixel(c); yy, xx = np.indices(F.shape[1:]); rr = np.hypot(xx - x0, yy - y0)
    sky = np.nanmedian(F[:, rr > 4], axis=1)
    ph = ((t - T0) / P) % 1
    d1 = np.abs(((ph - ph1 + 0.5) % 1) - 0.5); d2 = np.abs(((ph - ph2 + 0.5) % 1) - 0.5)
    mask = (d1 > 1.5 * hw1) & (d2 > 1.3 * hw2)
    gaps = np.where(np.diff(t) > 0.5)[0]; segs = np.split(np.arange(len(t)), gaps + 1)
    for r in (1.0, 1.5, 2.0, 3.0):
        f = np.nansum((F - sky[:, None, None])[:, rr < r], axis=1)
        fn_ = np.full_like(f, np.nan)
        for sg in segs:
            m = mask[sg]
            if m.sum() < 20: continue
            p = np.polyfit(t[sg][m] - t[sg].mean(), f[sg][m], 2)
            fn_[sg] = f[sg] / np.polyval(p, t[sg] - t[sg].mean())
        g = np.isfinite(fn_)
        base = np.nanmedian(fn_[mask & g])
        dep1 = base - np.nanmedian(fn_[g & (d1 < 0.3 * hw1)])
        dep2 = base - np.nanmedian(fn_[g & (d2 < 0.4 * hw2)])
        print(f"S{s:<5d} {r:4.1f}  {dep1*100:8.2f}  {dep2*1e3:10.2f}  {dep2/dep1:8.3f}")
