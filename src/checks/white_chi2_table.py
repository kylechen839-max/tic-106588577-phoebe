"""Recompute the photometry check, including the white-noise chi2 added on 2026-10-20, for
every stored PHOEBE best-fit bundle (outputs/candidates/*/*_phoebe_best*.phoebe) that has a
matching _phoebe_result*.json. Uses the stored 'final' model; no new PHOEBE computation.
Writes outputs/checks/white_chi2_table.csv. Offline (PHOEBE venv).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import phoebe

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
from run_phoebe_candidate import photometry_check  # noqa: E402

rows = []
for bf in sorted((ROOT / "outputs" / "candidates").glob("*/*_phoebe_best*.phoebe")):
    des = bf.parent.name
    tag = bf.name[len(f"{des}_phoebe_best"):-len(".phoebe")]
    rj = bf.parent / f"{des}_phoebe_result{tag}.json"
    if not rj.exists():
        continue
    res = json.loads(rj.read_text())
    eph = json.loads((bf.parent / f"{des}_ephemeris.json").read_text())
    try:
        b = phoebe.load(str(bf))
        pred = np.asarray(b.get_value(qualifier="fluxes", context="model", model="final"), float)
        obs = np.asarray(b.get_value(qualifier="fluxes", context="dataset", dataset="lc01"), float)
        sig = np.asarray(b.get_value(qualifier="sigmas", context="dataset", dataset="lc01"), float)
        times = np.asarray(b.get_value(qualifier="times", context="dataset", dataset="lc01"), float)
    except Exception as e:  # old bundles without a stored model
        print(f"{des}{tag}: skipped ({e.__class__.__name__})", flush=True)
        continue
    phase2 = eph.get("eclipses", {}).get("phase2", 0.5)
    nparams = len(res.get("params", [])) or 5
    ph = photometry_check(times, obs, pred, sig, res["t0"], res["period"], nparams, phase2)
    row = {"designation": des, "tag": tag or "(round 1)", "chi2_red": round(ph["chi2_red"], 3),
           "chi2_red_stored": round(res.get("chi2_final_hires", np.nan) / ph["dof"], 3),
           "white_noise_ppt": round(1e3 * ph.get("white_noise", np.nan), 3),
           "chi2_red_white": round(ph.get("chi2_red_white", np.nan), 1),
           "chi2_red_white_eclipse": round(ph.get("chi2_red_white_eclipse", np.nan), 1),
           "bias_prim": round(ph.get("primary_bias_sigma", np.nan), 2),
           "bias_sec": round(ph.get("secondary_bias_sigma", np.nan), 2),
           "passes_all": res.get("passes_all")}
    rows.append(row)
    print(row, flush=True)
out = ROOT / "outputs" / "checks" / "white_chi2_table.csv"
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(f"wrote {out} ({len(rows)} fits)")
