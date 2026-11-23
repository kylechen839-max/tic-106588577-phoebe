"""Binary-photosphere SED / W4 check with the fractional radii of a per-sector fit.

Absolute radii use the Kepler sma for M1 = main-sequence mass at T1 (assumed, as in run_phoebe_candidate) and the
fitted q. Photosphere = R1^2 B(T1) + R2^2 B(T2), scaled to J,H,Ks,W1 by run_phoebe_candidate.blackbody_check.
Usage: .venv-phoebe/bin/python src/checks/persector_sed.py DESIGNATION PERSECTOR_JSON TEFF1 BASE_RESULT_JSON
(BASE_RESULT_JSON: any _phoebe_result*.json of the star, for the period and Gaia distance.)
"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import run_phoebe_candidate as rp

name, fit_path, T1, base_path = sys.argv[1], Path(sys.argv[2]), float(sys.argv[3]), Path(sys.argv[4])
fit, base = json.load(open(fit_path)), json.load(open(base_path))
p = dict(zip(fit["params"], fit["v"]))
m1, _ = rp.ms_mass_radius(T1)
sma = rp.kepler_sma(base["period"], m1 * (1 + p["q"]))
R1, R2, T2 = p["rsum"] / (1 + p["k"]) * sma, p["rsum"] * p["k"] / (1 + p["k"]) * sma, T1 * p["teffratio"]
bc = rp.blackbody_check(name, R1, R2, T1, T2, base["blackbody_check"]["gaia_distance_pc"])
keep = {k: v for k, v in bc.items() if not isinstance(v, list)}
row = {"fit": fit_path.name, "T1": T1, "m1": m1, "q": p["q"], "sma_rsun": sma, "R1": R1, "R2": R2, "T2": T2, **keep}
print(json.dumps(row, indent=1, default=float))
out = fit_path.with_name(fit_path.stem + "_sed.json")
json.dump(row, open(out, "w"), indent=1, default=float); print(out)
