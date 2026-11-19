"""J2218: does a main-sequence tertiary (the 25-yr light-travel-time companion, m3 sin i >~ 1.5 Msun) change the
binary-photosphere SED check, the W4 excess or the photometric distance?

Photosphere = R1^2 B(T1) + R2^2 B(T2) [+ R3^2 B(T3)], scaled to J,H,Ks,W1 exactly as run_phoebe_candidate.blackbody_check.
Binary from _phoebe_result_t9500.json (T1 9500 K, R1 3.22, T2 6703 K, R2 1.10 Rsun). Tertiary: rough main-sequence
T and R per mass (assumed, +-10%). Also reports the tertiary's share of the TESS-band light (l3).
Usage: .venv-phoebe/bin/python src/checks/j2218_tertiary_sed.py
"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import run_phoebe_candidate as rp

N = "J221843.61+544715.0"
j = json.load(open(ROOT / f"outputs/candidates/{N}/{N}_phoebe_result_t9500.json"))
R1, R2, T1, T2 = j["R1_rsun"], j["R2_rsun"], j["teff1"], j["teff2"]
d_gaia = j["blackbody_check"]["gaia_distance_pc"]
base = rp.blackbody_check(N, R1, R2, T1, T2, d_gaia)
lam, f, s = map(np.array, (base["lam"], base["flux"], base["sigma"]))
MS = {0.0: (0, 0), 1.0: (5800, 1.0), 1.2: (6300, 1.2), 1.5: (7000, 1.5), 1.8: (7900, 1.75), 2.0: (8600, 1.9)}
tess = np.linspace(0.6, 1.0, 41)
rows = []
for m3, (T3, R3) in MS.items():
    shape = R1**2 * rp.planck(T1, lam) + R2**2 * rp.planck(T2, lam) + (R3**2 * rp.planck(T3, lam) if m3 else 0)
    ip = np.arange(4); w = s[ip] ** -2
    A = np.sum(w * shape[ip] * f[ip]) / np.sum(w * shape[ip] ** 2); model = A * shape
    chi2 = float(np.sum(((f[ip] - model[ip]) / s[ip]) ** 2))
    d = rp.RSUN_M * np.sqrt(np.pi * 1e26 / A) / rp.PC_M
    l3 = 0.0 if not m3 else float(np.trapezoid(R3**2 * rp.planck(T3, tess), tess) / np.trapezoid(
        R1**2 * rp.planck(T1, tess) + R2**2 * rp.planck(T2, tess) + R3**2 * rp.planck(T3, tess), tess))
    r = dict(m3=m3, T3=T3, R3=R3, l3_tess=l3, chi2_JHKW1=chi2, d_phot=d, d_ratio=d / d_gaia,
             **{f"{b}_ratio": float(f[i] / model[i]) for i, b in ((4, "W2"), (5, "W3"), (6, "W4"))},
             **{f"{b}_sigma": float((f[i] - model[i]) / s[i]) for i, b in ((4, "W2"), (5, "W3"), (6, "W4"))})
    rows.append(r)
    print(f"m3 {m3:.1f} (T {T3:5d} K, R {R3:.2f}): l3(TESS) {l3:.3f}  JHKW1 chi2 {chi2:6.1f}  d_phot/d_Gaia {d/d_gaia:.3f}  "
          f"W2 {r['W2_ratio']:.3f}x ({r['W2_sigma']:+.1f}s)  W3 {r['W3_ratio']:.2f}x ({r['W3_sigma']:+.1f}s)  "
          f"W4 {r['W4_ratio']:.1f}x ({r['W4_sigma']:+.1f}s)")
out = ROOT / f"outputs/candidates/{N}/timing/{N}_tertiary_sed.json"
json.dump(rows, open(out, "w"), indent=1); print(out)
