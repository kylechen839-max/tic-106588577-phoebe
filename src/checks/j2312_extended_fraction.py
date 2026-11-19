"""How much of J2312's 22 um flux can be extended dust? Fit the observed W4 radial profile as
   compact point source (empirical PSF)  +  extended component (PSF blurred by a Gaussian of FWHM g),
for a grid of g, and report the extended fraction and the implied compact (disk-like) excess.
Same tile, stamps and field-star PSF as j2312_w4_profile_plot.py. Total flux at the star is taken from
the nebula-modelled fits (11.8-12.5 mJy); photosphere 4.31 mJy (PHOEBE binary). 2-D least squares on the
background-subtracted stamp within 35".
python src/checks/j2312_extended_fraction.py"""
import importlib.util, sys
from pathlib import Path
import numpy as np
from scipy import ndimage
sys.stdout.reconfigure(encoding="utf-8")
here = Path(__file__).parent
src = open(here / "j2312_w4_profile_plot.py", encoding="utf-8").read().split("edges = np.arange")[0]
pp = {"__file__": str(here / "j2312_w4_profile_plot.py")}; exec(compile(src, "pp", "exec"), pp)  # reuse helpers without re-plotting
load, sin_pix, stamp, PIX, HALF = pp["load"], pp["sin_pix"], pp["stamp"], pp["PIX"], pp["HALF"]
img, h = load("w4"); x0, y0 = sin_pix(h, pp["RA"], pp["DEC"])
st0 = stamp(img, x0, y0)
# empirical 2-D PSF from isolated field stars (same selection as the profile plot)
bg = ndimage.zoom(ndimage.median_filter(img[::4, ::4], size=15), 4, order=1)[:img.shape[0], :img.shape[1]]
res = img - bg; noise = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))
mx = ndimage.maximum_filter(res, size=31); sat = 0.98 * np.nanmax(img); amp = st0[HALF, HALF]
pk = np.argwhere((res == mx) & (res > max(20 * noise, amp / 3)))
stamps = []
for y, x in pk:
    if not (HALF + 5 < y < img.shape[0] - HALF - 5 and HALF + 5 < x < img.shape[1] - HALF - 5): continue
    if np.hypot(x - x0, y - y0) * PIX < 300 or img[y, x] >= sat: continue
    s = stamp(img, float(x), float(y)); yy, xx = np.mgrid[:s.shape[0], :s.shape[1]]
    rr = np.hypot(xx - HALF, yy - HALF) * PIX
    if np.nanmax(ndimage.gaussian_filter(s, 2)[rr > 25]) > 0.12 * s[HALF, HALF]: continue
    stamps.append(s / np.nansum(s[rr < 35]))
psf = np.nanmedian(stamps, axis=0)
yy, xx = np.mgrid[:psf.shape[0], :psf.shape[1]]; rr = np.hypot(xx - HALF, yy - HALF) * PIX
m = rr < 35
psf /= np.nansum(psf[m])
data = st0 / np.nansum(st0[m])          # J2312 normalised to total flux within 35" = 1
print(f"{len(stamps)} field stars in the PSF")
best = None
for g in (4, 6, 8, 10, 13, 16, 20, 25, 30):
    ext = ndimage.gaussian_filter(psf, g / 2.3548 / PIX); ext /= np.nansum(ext[m])
    A = np.c_[psf[m], ext[m]]
    c, *_ = np.linalg.lstsq(A, data[m], rcond=None)
    c = np.clip(c, 0, None)
    chi = np.nansum((data[m] - A @ c) ** 2)
    f_ext = c[1] / c.sum()
    print(f"halo FWHM {g:>2}\" -> extended fraction {f_ext:5.2f}, rms resid {np.sqrt(chi / m.sum()):.2e}")
    if best is None or chi < best[0]: best = (chi, g, f_ext)
chi0 = np.nansum((data[m] - psf[m] * np.nansum(data[m] * psf[m]) / np.nansum(psf[m] ** 2)) ** 2)
print(f"point source alone: rms resid {np.sqrt(chi0 / m.sum()):.2e}")
_, g, f = best
for Ftot in (11.8, 12.3, 12.5):
    Fext = f * Ftot; Fc = Ftot - Fext
    print(f"F_total {Ftot} mJy: extended {Fext:.1f} mJy, compact {Fc:.1f} mJy -> compact excess {Fc - 4.31:.1f} mJy "
          f"({Fc / 4.31:.2f}x photosphere, {(Fc - 4.31) / 0.88:.1f} sigma)")
print(f"best halo FWHM {g}\", extended fraction {f:.2f}")

# Constrained scan: the photosphere (4.31 mJy) must be point-like, so compact fraction fc >= 4.31/F_total.
# Map rms residual over (fc, halo FWHM); "allowed" = within 10% of the best rms.
Ftot = 12.25
yy2, xx2 = np.mgrid[:psf.shape[0], :psf.shape[1]]
far = (rr > 25) & (rr < 35)
sig = 1.4826 * np.nanmedian(np.abs(data[far] - np.nanmedian(data[far])))
grid = []
for g in (8, 10, 13, 16, 20, 25, 30, 40):
    ext = ndimage.gaussian_filter(psf, g / 2.3548 / PIX); ext /= np.nansum(ext[m])
    for fc in np.arange(4.31 / Ftot, 1.0001, 0.01):
        model = fc * psf + (1 - fc) * ext
        grid.append((np.sqrt(np.nanmean((data[m] - model[m]) ** 2)), g, fc))
grid.sort()
rbest = grid[0][0]
ok = [x for x in grid if x[0] <= 1.10 * rbest]
fcs = [x[2] for x in ok]
print(f"\nconstrained (fc >= {4.31 / Ftot:.2f}): best rms {rbest:.2e} at halo {grid[0][1]}\", fc {grid[0][2]:.2f}; pixel noise ~{sig:.1e}")
print(f"within 10% of best rms: fc {min(fcs):.2f}-{max(fcs):.2f}, halo FWHM {min(x[1] for x in ok)}-{max(x[1] for x in ok)}\"")
for lab, fc in (("best", grid[0][2]), ("max allowed", max(fcs))):
    Fc = fc * Ftot
    print(f"  {lab}: compact {Fc:.1f} mJy -> compact excess {Fc - 4.31:.1f} mJy ({Fc / 4.31:.2f}x photosphere), "
          f"extended {Ftot - Fc:.1f} mJy")

# Figure: J2312 radial profile vs best compact+extended model and the max-compact allowed model
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
edges = np.arange(0, 36, PIX); rc = 0.5 * (edges[1:] + edges[:-1])
def prof(a):
    return np.array([np.nanmean(a[(rr >= p) & (rr < q)]) for p, q in zip(edges[:-1], edges[1:])])
def model(g, fc):
    ext = ndimage.gaussian_filter(psf, g / 2.3548 / PIX); ext /= np.nansum(ext[m])
    return fc * psf, (1 - fc) * ext
gb, fb = grid[0][1], grid[0][2]
okmax = max(ok, key=lambda x: x[2]); gm, fm = okmax[1], okmax[2]
pd_ = prof(data); norm = pd_[0]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
for ax, (g, fc, title) in zip(axs, ((gb, fb, f"best fit: compact {fb * Ftot:.1f} mJy (photosphere + {fb * Ftot - 4.31:.1f}), extended {Ftot * (1 - fb):.1f} mJy"),
                                    (gm, fm, f"most compact still allowed: compact {fm * Ftot:.1f} mJy (photosphere + {fm * Ftot - 4.31:.1f})"))):
    c_, e_ = model(g, fc)
    pc, pe = prof(c_) / norm, prof(e_) / norm
    ph = prof(psf * (4.31 / Ftot)) / norm
    ax.fill_between(rc, 0, ph, color="#cde2fb", label="photosphere (4.31 mJy, point source)")
    ax.plot(rc, pc, color="#2a78d6", lw=2, label=f"compact (point source)")
    ax.plot(rc, pe, color="#eb6834", lw=2, label=f"extended (PSF ⊗ {g}″ Gaussian)")
    ax.plot(rc, pc + pe, color="#0b0b0b", lw=1.5, ls="--", label="compact + extended")
    ax.plot(rc, pd_ / norm, "o", color="#0b0b0b", ms=4, label="J2312 W4 data")
    ax.set_title(title, fontsize=9.5, loc="left"); ax.set_xlabel("radius (arcsec)"); ax.set_xlim(0, 35)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axs[0].set_ylabel("brightness / central data value"); axs[0].legend(fontsize=8, frameon=False)
fig.suptitle("J2312 at 22 µm: how much of the light is compact vs extended (total 12.25 mJy at the star)", x=0.01, ha="left", fontsize=12, fontweight="bold")
fig.tight_layout(); fig.savefig(ROOT_OUT := Path(__file__).resolve().parents[2] / "outputs" / "j2312_contamination" / "j2312_compact_vs_extended.png", dpi=150)
print("saved", ROOT_OUT)
