#!/usr/bin/env python3
"""Kyle's "LightkurveBasicPipelineForTESS" (Colab, Google Drive) run on TESScut cubes from AWS.

The Colab notebook gets its cutouts with lightkurve.search_tesscut, which needs MAST (blocked in the
cloud). Here the same 10x10 px cutouts come from the public TESScut cubes on s3://stpubdata via
astrocut, and the notebook's steps are applied unchanged to every sector:
  - aperture: pf.create_threshold_mask(threshold=2, reference_pixel=(5, 5))
  - background regressors: pixels [8:10, 5:10]; DesignMatrix(...).pca(5).append_constant()
  - RegressionCorrector; corrected = raw - model + median(model)
  - no time cuts, no flatten
Each sector is then divided by its median so sectors can be combined. Writes <outdir>/<des>_S<nn>_lk.txt
(time, flux, flux_err, normalised) and <outdir>/<des>_lk_all.txt.

  .venv-tess/bin/python src/pipeline/lightkurve_pipeline_s3.py --designation J013614.42+625507.1 \
      --ra 24.0600938 --dec 62.9186539 --tic 390117626
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
import lightkurve as lk
from astrocut import cube_cut
from astropy.coordinates import SkyCoord
from tess_stars2px import tess_stars2px_function_entry

ROOT = Path(__file__).resolve().parents[2]
CUBE = "s3://stpubdata/tess/public/mast/tess-s{:04d}-{}-{}-cube.fits"

ap = argparse.ArgumentParser()
ap.add_argument("--designation", required=True)
ap.add_argument("--ra", type=float, required=True)
ap.add_argument("--dec", type=float, required=True)
ap.add_argument("--tic", type=int, default=0)
ap.add_argument("--size", type=int, default=10)
ap.add_argument("--threshold", type=float, default=2)
ap.add_argument("--components", type=int, default=5)
ap.add_argument("--outdir", default=None)
args = ap.parse_args()
name = args.designation
out = Path(args.outdir) if args.outdir else ROOT / "outputs" / "candidates" / name / "lk_pipeline"
cdir = out / "cutouts"
cdir.mkdir(parents=True, exist_ok=True)
coord = SkyCoord(args.ra, args.dec, unit="deg")
res = tess_stars2px_function_entry(args.tic, args.ra, args.dec)
allrows = []
for sec, cam, ccd in zip(res[3], res[4], res[5]):
    sec, cam, ccd = int(sec), int(cam), int(ccd)
    hits = list(cdir.glob(f"tess-s{sec:04d}-{cam}-{ccd}_*.fits"))
    try:
        fn = hits[0] if hits else Path(cube_cut(CUBE.format(sec, cam, ccd), coord, args.size, output_path=str(cdir), verbose=False))
    except Exception as exc:
        print(f"S{sec}: no cube ({str(exc).splitlines()[0]})")
        continue
    pf = lk.read(str(fn))
    t_mask = pf.create_threshold_mask(threshold=args.threshold, reference_pixel=(args.size // 2, args.size // 2))
    bg_mask = np.zeros(pf.shape[1:], dtype=bool)
    bg_mask[8:10, 5:10] = True
    lc_raw = pf.to_lightcurve(aperture_mask=t_mask)
    good = np.isfinite(lc_raw.flux.value) & np.isfinite(lc_raw.time.value)
    lc_raw, pfg = lc_raw[good], pf[good]
    regs = np.nan_to_num(pfg.flux[:, bg_mask].value)
    dm = lk.DesignMatrix(regs, name="background_regressors").pca(args.components).append_constant()
    rc = lk.RegressionCorrector(lc_raw)
    rc.correct(dm)
    corr = lc_raw - rc.model_lc + np.median(rc.model_lc.flux)
    t, f, e = corr.time.value, corr.flux.value, corr.flux_err.value
    ok = np.isfinite(f)
    med = np.median(f[ok])
    np.savetxt(out / f"{name}_S{sec:02d}_lk.txt", np.c_[t[ok], f[ok], e[ok], f[ok] / med],
               header="time flux flux_err flux_norm (lightkurve RegressionCorrector, Kyle's Colab settings)", fmt="%.6f")
    allrows.append(np.c_[t[ok], f[ok] / med, e[ok] / med, np.full(ok.sum(), sec)])
    print(f"S{sec}: {t_mask.sum()} px aperture, {ok.sum()} cadences, median flux {med:.1f} e/s")
np.savetxt(out / f"{name}_lk_all.txt", np.vstack(allrows), header="time flux flux_err sector (each sector / its median)", fmt="%.6f")
