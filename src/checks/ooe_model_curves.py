#!/usr/bin/env python3
"""Out-of-eclipse model light curves at a per-sector fit's best parameters, sphere vs Roche.

The per-sector fitter uses spherical stars for P > 3 d, which have no ellipsoidal variation.
This evaluates the same parameters with both distortion methods on a fine phase grid, so the
data's out-of-eclipse shape can be compared with what tidal distortion and reflection predict.

  .venv-phoebe/bin/python src/checks/ooe_model_curves.py --designation J231201.40+532028.6 \
      --persector outputs/candidates/J231201.40+532028.6/J231201.40+532028.6_persector_persectorffi2.json \
      --teff 9932 --no-irrad
Writes <cand>/<des>_ooe_models.npz (phase, one flux array per variant).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
import run_phoebe_candidate as rpc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--designation", required=True)
ap.add_argument("--persector", required=True)
ap.add_argument("--teff", type=float, required=True)
ap.add_argument("--q", type=float, default=0.8)
ap.add_argument("--ntri", type=int, default=1500)
ap.add_argument("--no-irrad", action="store_true")
args = ap.parse_args()

name = args.designation
cdir, eph, *_ = rpc.load_inputs(name)
period, t0 = eph["period"], eph["t0"]
sf = json.loads(Path(args.persector).read_text())
rpc.PARAMS[:] = list(sf["params"])
m1, _ = rpc.ms_mass_radius(args.teff)
phase = np.linspace(-0.5, 0.5, 401)
times = t0 + phase * period
out = {"phase": phase}
variants = [("sphere", "sphere", not args.no_irrad), ("roche", "roche", not args.no_irrad)]
if args.no_irrad:
    variants.append(("roche_irrad", "roche", True))
for label, dist, irrad in variants:
    m = rpc.Model(times, np.ones_like(times), np.full_like(times, 1e-3), period, t0, args.teff, m1, args.q,
                  dist, args.ntri, irrad=irrad)
    m.b.set_value("pblum_mode@lc01", "component-coupled")
    m.set(sf["v"])
    m.b.run_compute(compute="phoebe01", model="ooemodel", overwrite=True, progressbar=False)
    f = np.asarray(m.b.get_value("fluxes@lc01@ooemodel@model"), float)
    out[label] = f / np.median(f[np.abs(phase) > 0.15])
    print(label, "min", out[label].min().round(4), flush=True)
np.savez(cdir / f"{name}_ooe_models.npz", **out)
