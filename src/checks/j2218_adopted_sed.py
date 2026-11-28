"""J2218: binary-photosphere SED / W4 check redone with the adopted light-curve model (_persectorffi_albs_final).

Fractional radii (rsum, k) and T2/T1 come from the per-sector fit. Absolute radii use the Kepler sma for
M1 = main-sequence mass at T1 = 9500 K (assumed, as in run_phoebe_candidate) and the fitted q
(the fit itself ran with sma from --q 0.8; only the fractional radii are constrained by the light curve).
Photosphere = R1^2 B(T1) + R2^2 B(T2), scaled to J,H,Ks,W1 with run_phoebe_candidate.blackbody_check.
Usage: .venv-phoebe/bin/python src/checks/j2218_adopted_sed.py
"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import run_phoebe_candidate as rp

N = "J221843.61+544715.0"
C = ROOT / f"outputs/candidates/{N}"
old = json.load(open(C / f"{N}_phoebe_result_t9500.json"))
fit = json.load(open(C / f"{N}_persector_persectorffi_albs_final.json"))
p = dict(zip(fit["params"], fit["v"]))
eph = json.load(open(C / "ffi" / f"{N}_ephemeris.json"))
T1 = 9500.0
m1, _ = rp.ms_mass_radius(T1)
rows = []
for label, q in (("adopted q", p["q"]), ("fit's sma (q 0.8)", 0.8)):
    sma = rp.kepler_sma(eph["period"], m1 * (1 + q))
    R1 = p["rsum"] / (1 + p["k"]) * sma
    R2 = p["rsum"] * p["k"] / (1 + p["k"]) * sma
    T2 = T1 * p["teffratio"]
    bc = rp.blackbody_check(N, R1, R2, T1, T2, old["blackbody_check"]["gaia_distance_pc"])
    keep = {k: v for k, v in bc.items() if not isinstance(v, list)}
    rows.append({"label": label, "q": q, "m1": m1, "sma_rsun": sma, "R1": R1, "R2": R2, "T2": T2, **keep})
    print(f"{label}: q {q:.3f} sma {sma:.2f} R1 {R1:.2f} R2 {R2:.2f} T2 {T2:.0f}  ->  " +
          ", ".join(f"{k} {v:.3g}" if isinstance(v, float) else f"{k} {v}" for k, v in keep.items()))
print("old (_t9500):", {k: v for k, v in old["blackbody_check"].items() if not isinstance(v, list)})
out = C / "timing" / f"{N}_adopted_sed.json"
json.dump(rows, open(out, "w"), indent=1, default=float); print(out)
