"""J2312 false-alarm test: how often does the nebula alone mimic a ~7 mJy point source?
After removing the 12 um-map nebula model (method B, as in j2312_profile_nebula_removed.py), fit
A*PSF + constant within 25" at random positions on the ridge/nebula (no catalogued W4 or W3 point
source within 25"), and at J2312. Fluxes are scaled to mJy with the same catalogue-like calibration.
Also reports the nebula model's own flux inside the 16.5" aperture at the star.
python src/checks/j2312_false_alarm.py"""
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
sys.stdout.reconfigure(encoding="utf-8")
here = Path(__file__).parent
src = open(here / "j2312_profile_nebula_removed.py", encoding="utf-8").read().split("star_img = s4 - neb")[0]
g = {"__file__": str(here / "j2312_profile_nebula_removed.py")}
exec(compile(src, "nr", "exec"), g)
s4, neb, P4, psf4, R, cx, cy, PIX, H, FB, PHOT, ECAT, cB, c0 = (g[k] for k in
    ("s4", "neb", "P4", "psf4", "R", "cx", "cy", "PIX", "H", "FB", "PHOT", "ECAT", "cB", "c0"))
place = g["place"]
mjy_per_amp = 15.36 / c0[0]          # same scaling as FB = 15.36*cB[0]/c0[0]
res = s4 - neb
yy, xx = np.mgrid[:s4.shape[0], :s4.shape[1]]

# nebula model flux inside the 16.5" aperture at the star, relative to its surroundings
ap = R < 16.5
neb_ap = np.nansum((neb - np.nanmedian(neb[(R > 140) & (R < 160)]))[ap]) * mjy_per_amp / np.nansum(P4[ap])
print(f"nebula model inside the 16.5\" aperture at J2312: {neb_ap:.1f} mJy-equivalent (vs star point source {FB:.1f} mJy)")

# compact sources to avoid: local maxima of the residual > 5 sigma (other stars), plus J2312
sm = ndimage.gaussian_filter(res, 2)
noise = 1.4826 * np.nanmedian(np.abs(res[R > 100] - np.nanmedian(res[R > 100])))
mx = ndimage.maximum_filter(sm, size=11)
pk = np.argwhere((sm == mx) & (sm > 5 * noise))
bad = np.zeros_like(R, bool)
for y, x in pk:
    bad |= np.hypot(xx - x, yy - y) * PIX < 25
bad |= R < 40

def fit_at(x, y):
    P = place(psf4, s4.shape, x, y)
    m = np.hypot(xx - x, yy - y) * PIX < 25
    A, *_ = np.linalg.lstsq(np.c_[P[m], np.ones(m.sum())], res[m], rcond=None)
    return A[0] * mjy_per_amp

rng = np.random.default_rng(1)
nebz = neb - np.nanmedian(neb[(R > 140) & (R < 160)])
star_neb = nebz[int(round(cy)), int(round(cx))]
out = {"ridge (nebula model >= 50% of its level at J2312)": [], "whole field": []}
inside = (yy > H + 1) & (yy < s4.shape[0] - H - 2) & (xx > H + 1) & (xx < s4.shape[1] - H - 2)
cand = np.argwhere(~bad & inside & (R > 40))
for y, x in cand[rng.choice(len(cand), min(3000, len(cand)), replace=False)]:
    f = fit_at(float(x), float(y))
    out["whole field"].append(f)
    if nebz[y, x] >= 0.5 * star_neb:
        out["ridge (nebula model >= 50% of its level at J2312)"].append(f)
fJ = fit_at(cx, cy)
excess = FB - PHOT
print(f"J2312 point-source flux in the nebula-subtracted image (25\" fit): {fJ:.2f} mJy; excess over photosphere {fJ - PHOT:.2f} mJy")
for k, v in out.items():
    v = np.array(v)
    print(f"{k}: n = {len(v)}, fake point flux median {np.median(v):+.2f} mJy, 16-84% [{np.percentile(v,16):+.2f}, {np.percentile(v,84):+.2f}], "
          f"max {v.max():+.2f}; fraction >= {fJ - PHOT:.1f} mJy: {(v >= fJ - PHOT).mean():.4f}; >= 3 mJy: {(v >= 3).mean():.4f}")
    if len(v):
        sd = 1.4826 * np.median(np.abs(v - np.median(v)))
        print(f"   empirical significance of J2312's excess vs this distribution: {(fJ - PHOT - np.median(v)) / sd:.1f} sigma (robust sd {sd:.2f} mJy)")

# The excluded >5-sigma peaks: are they stars (AllWISE W1/W2 counterparts) or nebular knots?
import warnings; warnings.filterwarnings("ignore")
import pandas as pd
cat = pd.read_csv(Path(g["D"]) / "J231201.40+532028.6_allwise_3arcmin.csv")   # AllWISE 3' cone saved 2026-10-26
cat = cat.rename(columns={"ra": "RAJ2000", "dec": "DEJ2000", "w1mpro": "W1mag", "w2mpro": "W2mag", "w3mpro": "W3mag", "w4mpro": "W4mag", "w4snr": "snr4"})
dec0 = np.radians(g["DEC"])
cat["dx"] = -(cat.RAJ2000 - g["RA"]) * np.cos(dec0) * 3600   # arcsec, +x = west on this SIN grid (CDELT1 < 0)
cat["dy"] = (cat.DEJ2000 - g["DEC"]) * 3600
print("\nexcluded peaks (>5 sigma) within 100\" of J2312, other than J2312 itself:")
for y, x in pk:
    dxp, dyp = (x - cx) * PIX, (y - cy) * PIX
    rr_ = np.hypot(dxp, dyp)
    if rr_ < 10 or rr_ > 100 or not (H + 1 < y < s4.shape[0] - H - 2 and H + 1 < x < s4.shape[1] - H - 2):
        continue
    f = fit_at(float(x), float(y))
    d = np.hypot(cat.dx - (-dxp), cat.dy - dyp)   # image x increases west? compare both senses below
    d2 = np.hypot(cat.dx - dxp, cat.dy - dyp)
    j = int(np.argmin(np.minimum(d, d2))); sep = float(min(d[j], d2[j]))
    c_ = cat.iloc[j]
    print(f"  at ({dxp:+.0f}\", {dyp:+.0f}\") r={rr_:.0f}\": point flux {f:.2f} mJy; nearest AllWISE {sep:.1f}\": "
          f"W1 {c_.W1mag:.2f} W2 {c_.W2mag:.2f} W3 {c_.W3mag:.2f} W4 {c_.W4mag:.2f} (W4 S/N {c_.snr4})")
