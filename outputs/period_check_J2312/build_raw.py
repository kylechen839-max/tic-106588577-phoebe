"""Rebuild J2312 _raw.txt from MAST FITS, mirroring src/pipeline/fetch_tess_lc.py (download + normalise)."""
import sys, glob, os
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
dl, outfile = sys.argv[1], sys.argv[2]
LOCALCLIP = "--no-local-clip" not in sys.argv
# lightkurve TessQualityFlags.DEFAULT_BITMASK
DEFAULT = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 512 | 16384 # bits used by lightkurve "default"
def robust_sigma(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))
segs = {}
for fn in sorted(glob.glob(os.path.join(dl, "*.fits"))):
    with fits.open(fn) as h:
        d = h[1].data; sec = int(h[0].header["SECTOR"])
        if "qlp" in fn:
            f, e, author = d["KSPSAP_FLUX"], d["KSPSAP_FLUX_ERR"], "QLP"
        else:
            f, e, author = d["PDCSAP_FLUX"], d["PDCSAP_FLUX_ERR"], "TESS-SPOC"
        t = d["TIME"]; q = d["QUALITY"]
        m = (q & DEFAULT) == 0
        t, f, e = np.asarray(t[m], float), np.asarray(f[m], float), np.asarray(e[m], float)
        good = np.isfinite(t) & np.isfinite(f)
        segs[sec] = (t[good], f[good], e[good], author)
        print(f"sector {sec}: {author} {good.sum()} points")
T, F, E, Sx = [], [], [], []
for sector, (t, f, e, _) in sorted(segs.items()):
    o = np.argsort(t); t, f, e = t[o], f[o], e[o]
    breaks = np.where(np.diff(t) > 0.75)[0] + 1
    for idx in np.split(np.arange(len(t)), breaks):
        if len(idx) < 30: continue
        ff = f[idx]; med = np.nanmedian(ff); ff = ff / med; ee = e[idx] / med
        s = robust_sigma(ff); keep = ff < 1 + 5 * s
        local = median_filter(ff, size=9, mode="nearest"); dev = ff - local
        if LOCALCLIP: keep &= np.abs(dev) < 6 * robust_sigma(dev)
        keep &= ff > 0.05
        T.append(t[idx][keep]); F.append(ff[keep]); E.append(ee[keep]); Sx.append(np.full(keep.sum(), sector))
t, f, e, s = map(np.concatenate, (T, F, E, Sx))
bad = ~np.isfinite(e) | (e <= 0); e[bad] = robust_sigma(f)
np.savetxt(outfile, np.c_[t, f, e, s], header="time flux flux_err sector", fmt="%.8f")
print("total", len(t), "bytes", os.path.getsize(outfile))
