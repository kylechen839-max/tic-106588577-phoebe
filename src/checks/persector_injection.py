"""Injection test for src/pipeline/fit_per_sector_dilution.py.

Builds a synthetic per-point light curve on J1122's real time stamps and sectors from a PHOEBE
model with known parameters, applies a known per-sector dilution (flux = a_s + b_s * model) and
white noise at each sector's real scatter, and writes it in the --raw format. Running the fitter
on it from perturbed parameters should recover the parameters and the depth scales b_s/mean(b).
Writes outputs/checks/persector_injection_raw.txt and _truth.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import run_phoebe_candidate as rpc  # noqa: E402

DES = "J112238.89-592027.5"
RAW = ROOT / "outputs" / "j1122_tess_multisector_s10_s11_s37" / f"{DES}_LightCurve_multisector.txt"
TRUE = {"incl": 89.5, "dt0": 0.0026, "rsum": 0.134, "k": 0.32, "teffratio": 0.59}
SCALE = {10: 0.70, 11: 1.40, 37: 0.90}   # injected depth scales (like the J1122 fit)
TEFF = 7705.0

cdir, eph, times, flux, sig = rpc.load_inputs(DES)
raw = np.genfromtxt(RAW, comments="#")
raw = raw[np.isfinite(raw).all(axis=1)]
t, s = raw[:, 0], raw[:, 3].astype(int)
rpc.PARAMS[:] = list(TRUE)
m1, _ = rpc.ms_mass_radius(TEFF)
# model on a fine phase grid, interpolated to the real time stamps
P, t0 = eph["period"], eph["t0"]
grid = t0 + np.linspace(-0.2, 0.8, 2001) * P
m = rpc.Model(grid, np.ones_like(grid), np.full_like(grid, 1e-3), P, t0, TEFF, m1, 0.8, "sphere", 1500, l3=0.0)
m.set([TRUE[k] for k in TRUE])
m.b.set_value("pblum_mode@lc01", "component-coupled")
m.b.run_compute(compute="phoebe01", model="truth", overwrite=True, progressbar=False)
mod = np.asarray(m.b.get_value("fluxes@lc01@truth@model"), float)
mod = mod / np.median(mod)
gph = rpc.phase_values(grid, t0, P)
o = np.argsort(gph)
ph = rpc.phase_values(t, t0, P)
mt = np.interp(ph, gph[o], mod[o])
rng = np.random.default_rng(7)
f = np.empty_like(mt)
for sec, bscale in SCALE.items():
    k = s == sec
    noise = 1.4826 * np.median(np.abs(np.diff(raw[k, 1]))) / np.sqrt(2)
    f[k] = (1 - bscale) + bscale * mt[k] + rng.normal(0, noise, k.sum())
out = ROOT / "outputs" / "checks"
out.mkdir(parents=True, exist_ok=True)
np.savetxt(out / "persector_injection_raw.txt", np.c_[t, f, raw[:, 2], s], header="time flux flux_err sector", fmt="%.8f")
(out / "persector_injection_truth.json").write_text(json.dumps(
    {"params": list(TRUE), "v": list(TRUE.values()), "depth_scale": {str(k): v / np.mean(list(SCALE.values())) for k, v in SCALE.items()},
     "note": "start the fitter from perturbed params (persector_injection_start.json)"}, indent=1))
start = {"params": list(TRUE), "v": [89.0, 0.0020, 0.128, 0.36, 0.55]}
(out / "persector_injection_start.json").write_text(json.dumps(start))
print("wrote injection files; true depth scales", {k: round(v / np.mean(list(SCALE.values())), 3) for k, v in SCALE.items()})
