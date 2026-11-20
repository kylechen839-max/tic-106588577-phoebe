"""Tier A: binary photosphere, distance check and dust location vs. binary orbit, from existing PHOEBE fits.

For each tier-A object take the adopted PHOEBE fit (phoebe_fit column of candidates_all_in_one.csv),
read its blackbody_check (two-star photosphere R1^2 B(T1) + R2^2 B(T2) scaled to J,H,Ks,W1) and
geometry, and place the single-temperature dust at its blackbody equilibrium radius
    r_dust = 0.5 * sqrt(R1^2 T1^4 + R2^2 T2^4) / T_dust^2
(grey grains; small grains run hotter, so this is a lower limit on the distance). Compare with the
separation a: circumbinary orbits are stable beyond about 2-4 a (Holman & Wiegert 1999).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "tierA_geometry"
RSUN_AU = 1 / 215.032

# J2312's adopted light-curve model is the per-sector Roche fit; its SED check lives in this shared-fit file
OVERRIDE = {"J231201.40+532028.6": "_gteffnoirreccesm003"}


def main():
    cands = pd.read_csv(ROOT / "outputs" / "crossmatch" / "candidates_all_in_one.csv")
    rows = []
    for _, c in cands[cands.tier == "A"].iterrows():
        n = c.designation
        tag = OVERRIDE.get(n, str(c.phoebe_fit))
        tag = "" if tag == "(base)" else tag
        j = json.load(open(ROOT / "outputs" / "candidates" / n / f"{n}_phoebe_result{tag}.json"))
        bb = j.get("blackbody_check") or {}
        ex = bb.get("excess_vs_binary_photosphere", {})
        R1, R2, T1, T2 = j["R1_rsun"], j["R2_rsun"], j["teff1"], j["teff2"]
        td = float(c.get("dust_temp", np.nan))
        a = j["sma_rsun"]
        r_d = 0.5 * np.sqrt(R1 ** 2 * T1 ** 4 + R2 ** 2 * T2 ** 4) / td ** 2 if np.isfinite(td) else np.nan
        rows.append({
            "designation": n, "fit": tag or "(base)", "lc_chi2_red": c.phoebe_chi2_red,
            "P_d": j["period"], "incl": j["incl"], "ecc": j.get("ecc", 0.0), "q_assumed": j.get("q", j.get("q_assumed")),
            "T1": T1, "T2": T2, "R1": R1, "R2": R2, "a_rsun": a, "a_au": a * RSUN_AU,
            "binphot_chi2_jhkw1": bb.get("binary_photosphere_chi2_JHKW1"),
            "d_ratio": bb.get("distance_ratio_phot_over_gaia"),
            "W3_sigma_bin": ex.get("W3", {}).get("sigma"), "W4_sigma_bin": ex.get("W4", {}).get("sigma"),
            "W4_ratio_bin": ex.get("W4", {}).get("ratio"),
            "excess_survives": bb.get("pass_excess_survives"),
            "T_dust": td, "dust_at_grid_floor": bool(td <= 40.5) if np.isfinite(td) else None,
            "r_dust_au": r_d * RSUN_AU, "r_dust_over_a": r_d / a,
        })
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "tierA_geometry.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
