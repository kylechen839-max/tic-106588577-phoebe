#!/usr/bin/env python3
"""TESS FFI light curves from the public AWS cubes, for when MAST is unreachable.

Cuts a small stamp out of the TESScut FFI cubes on s3://stpubdata (anonymous,
byte-range reads, a few seconds per sector), does aperture photometry with a
per-cadence background, drops scattered-light cadences, removes slow trends
per orbit while keeping the periodic signal, then hands the segments to
fetch_tess_lc.py's cleaning, period refinement and sector-consensus binning.

  .venv-tess/bin/python src/pipeline/fetch_tess_ffi_s3.py --designation J013614.42+625507.1 \
      --ra 24.0600938 --dec 62.9186539 --tic 390117626 --period 2.78951

Needs astrocut, s3fs and tess-point in .venv-tess. Outputs as fetch_tess_lc.py,
plus <name>_ffi_stamps.png (median stamps with apertures) and the cutouts in
outputs/candidates/<name>/ffi_cutouts/.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
import lightkurve as lk
from astrocut import cube_cut
from astropy.coordinates import SkyCoord
from tess_stars2px import tess_stars2px_function_entry

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fetch_tess_lc as base  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CUBE = "s3://stpubdata/tess/public/mast/tess-s{:04d}-{}-{}-cube.fits"
SIZE = 15
BG_EXCESS = 0.08   # drop cadences whose background rises > 8% of the star's flux above its quiet level
APER_FRAC = 0.3    # aperture: pixels within 1.6 px of the target above 30% of the peak (background-subtracted)


def cutouts(ra, dec, tic, outdir):
    tic_id = int(tic) if tic and tic != "nan" else 0
    res = tess_stars2px_function_entry(tic_id, ra, dec)
    sectors, cams, ccds = res[3], res[4], res[5]
    coord = SkyCoord(ra, dec, unit="deg")
    outdir.mkdir(parents=True, exist_ok=True)
    files = {}
    for sec, cam, ccd in zip(sectors, cams, ccds):
        sec, cam, ccd = int(sec), int(cam), int(ccd)
        if sec < 1:
            continue
        hits = list(outdir.glob(f"tess-s{sec:04d}-{cam}-{ccd}_*.fits"))
        if hits:
            files[sec] = hits[0]
            continue
        try:
            files[sec] = Path(cube_cut(CUBE.format(sec, cam, ccd), coord, SIZE, output_path=str(outdir), verbose=False))
            print(f"  sector {sec}: cut from camera {cam} CCD {ccd}")
        except Exception as exc:
            print(f"  sector {sec}: no cube ({str(exc).splitlines()[0]})")
    return files


def photometry(fn, ra, dec):
    tpf = lk.read(str(fn), quality_bitmask="default")
    fl = tpf.flux.value
    img = np.nanmedian(fl, axis=0)
    x, y = tpf.wcs.world_to_pixel_values(ra, dec)
    yy, xx = np.mgrid[: img.shape[0], : img.shape[1]]
    r = np.hypot(xx - x, yy - y)
    valid = np.isfinite(img) & (img > 0)
    far = (r > 3) & valid
    bgpix = far & (img < np.nanpercentile(img[far], 50))
    bg = np.nanmedian(fl[:, bgpix], axis=1)
    sub = img - np.nanmedian(img[bgpix])
    near = (r < 1.6) & valid
    ap = near & (sub > APER_FRAC * sub[near].max())
    f = np.nansum(fl[:, ap] - bg[:, None], axis=1)
    t = tpf.time.value
    good = np.isfinite(t) & np.isfinite(f) & np.isfinite(bg)
    return t[good], f[good], bg[good], img, ap, (x, y)


def detrend(t, f, bg, period, gap=0.75, deg=2, nharm=4):
    """Per orbit: drop scattered-light cadences, then fit poly(time) + Fourier(phase) and divide out the poly only."""
    keep = np.ones(len(t), bool)
    trend = np.ones(len(t))
    breaks = np.where(np.diff(t) > gap)[0] + 1
    for idx in np.split(np.arange(len(t)), breaks):
        if len(idx) < 30:
            keep[idx] = False
            continue
        quiet = np.nanpercentile(bg[idx], 20)
        keep[idx] &= (bg[idx] - quiet) < BG_EXCESS * np.nanmedian(f[idx])
    for idx in np.split(np.arange(len(t)), breaks):
        idx = idx[keep[idx]]
        if len(idx) < 30:
            keep[idx] = False
            continue
        tt = (t[idx] - t[idx].mean()) / max(np.ptp(t[idx]), 1e-3)
        ph = 2 * np.pi * (t[idx] / period)
        cols = [tt ** k for k in range(deg + 1)]
        cols += [fn(h * ph) for h in range(1, nharm + 1) for fn in (np.cos, np.sin)]
        A = np.vstack(cols).T
        y = f[idx] / np.nanmedian(f[idx])
        m = np.ones(len(idx), bool)
        for _ in range(5):  # iterate, rejecting eclipse points (low) and outliers
            coef, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
            res = y - A @ coef
            s = base.robust_sigma(res[m])
            m = (res > -2.5 * s) & (res < 4 * s)
        poly = np.vstack(cols[: deg + 1]).T @ coef[: deg + 1]
        trend[idx] = poly / np.median(poly) * np.nanmedian(f[idx])
    return keep, trend


def make_download(ra, dec, tic, period, outdir, name):
    def download(target, coord, max_sectors=12):
        files = cutouts(ra, dec, tic, outdir / "ffi_cutouts")
        out = {}
        n = len(files)
        fig, axes = plt.subplots(1, max(n, 1), figsize=(3.2 * max(n, 1), 3.4), squeeze=False)
        for ax, (sec, fn) in zip(axes[0], sorted(files.items())[:max_sectors]):
            t, f, bg, img, ap, (x, y) = photometry(fn, ra, dec)
            keep, trend = detrend(t, f, bg, period)
            fd = f / trend
            fd = fd[keep]
            print(f"  sector {sec}: {ap.sum()} px aperture, {keep.sum()}/{len(t)} cadences kept")
            out[sec] = (t[keep], fd, np.full(keep.sum(), np.nan), "FFI-S3")
            ax.imshow(np.log10(np.clip(img, 1, None)), origin="lower", cmap="viridis")
            ax.contour(ap, levels=[0.5], colors="r", linewidths=1)
            ax.plot(x, y, "w+", ms=10)
            ax.set_title(f"S{sec}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{name}: median FFI stamps, aperture in red")
        fig.tight_layout(); fig.savefig(outdir / f"{name}_ffi_stamps.png", dpi=120); plt.close(fig)
        return out
    return download


if __name__ == "__main__":
    argv = sys.argv[1:]
    get = lambda k: argv[argv.index(k) + 1]
    name, ra, dec, period = get("--designation"), float(get("--ra")), float(get("--dec")), float(get("--period"))
    tic = get("--tic") if "--tic" in argv else None
    outdir = ROOT / "outputs" / "candidates" / name
    base.download = make_download(ra, dec, tic, period, outdir, name)
    base.main()
