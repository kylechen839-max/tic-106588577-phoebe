"""Difference-image centroid for a feature at a given phase (eclipse or other dip), per sector.

Usage: dip_centroid.py CUTOUT_GLOB RA DEC P T0 PHASE HALFWIDTH [OOT_LO OOT_HI]
In-feature cadences: |phase - PHASE| < HALFWIDTH. Comparison cadences: phase offset from PHASE
between OOT_LO and OOT_HI (default 1.2x and 2.2x HALFWIDTH, both sides). Reports the S/N-weighted
centroid of (out - in) within 2.5 px of the target, its offset from the target and a bootstrap error,
and the peak-pixel S/N.
"""
import glob, re, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
import lightkurve as lk
from astropy.coordinates import SkyCoord

pat, ra, dec, P, T0, ph0, hw = sys.argv[1], *map(float, sys.argv[2:8])
lo = float(sys.argv[8]) if len(sys.argv) > 8 else 1.2 * hw
hi = float(sys.argv[9]) if len(sys.argv) > 9 else 2.2 * hw
c = SkyCoord(ra, dec, unit="deg")
for fn in sorted(glob.glob(pat), key=lambda x: int(re.search(r"s0*(\d+)-", x).group(1))):
    s = int(re.search(r"s0*(\d+)-", fn).group(1))
    pf = lk.read(fn); F = pf.flux.value; t = pf.time.value
    ok = np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1) & (pf.quality == 0)
    F, t = F[ok], t[ok]
    d = ((t - T0) / P - ph0 + 0.5) % 1 - 0.5
    ine = np.abs(d) < hw; oot = (np.abs(d) > lo) & (np.abs(d) < hi)
    if ine.sum() < 5 or oot.sum() < 5:
        print(f"S{s}: too few cadences"); continue
    x0, y0 = pf.wcs.world_to_pixel(c); yy, xx = np.indices(F.shape[1:])
    near = np.hypot(xx - x0, yy - y0) < 2.5
    def diff(i, o): return np.nanmedian(F[o], 0) - np.nanmedian(F[i], 0)
    def cen(dd):
        w = np.clip(dd, 0, None) * near; return (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    dif = diff(ine, oot)
    noise = 1.4826 * np.nanmedian(np.abs(F[oot] - np.nanmedian(F[oot], 0)), 0) / np.sqrt(ine.sum())
    snr = dif / noise
    cx, cy = cen(dif); rng = np.random.default_rng(1); ii, oo = np.where(ine)[0], np.where(oot)[0]
    bs = np.array([cen(np.nanmedian(F[rng.choice(oo, len(oo))], 0) - np.nanmedian(F[rng.choice(ii, len(ii))], 0)) for _ in range(200)])
    j, i = np.unravel_index(np.nanargmax(np.where(near, snr, -np.inf)), snr.shape)
    print(f"S{s}: offset {np.hypot(cx-x0, cy-y0)*21:.0f}\" +/- {np.hypot(*bs.std(0))*21:.0f}\"; peak pixel S/N {snr[j,i]:.1f} at {np.hypot(i-x0,j-y0):.1f} px; n_in {ine.sum()}")
