"""Offline YSO / dipper grading from cached 2MASS + AllWISE photometry and Gaia values.

Reports the Koenig & Leisawitz (2014) WISE colour class, the 2-22 um spectral index alpha (K to W4),
the absolute G magnitude with a rough extinction from (BP-RP), and for a given dip period the
Keplerian corotation radius and the blackbody dust temperature there.

Usage: yso_colour_grade.py DESIGNATION TEFF RSTAR MSTAR PERIOD_D
(photometry from outputs/crossmatch/dd1_support_cache/allwise.csv, Gaia from dd1_candidates_scored.csv)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LAM = {"K": 2.159, "W1": 3.353, "W2": 4.603, "W3": 11.56, "W4": 22.09}  # um
F0 = {"K": 666.7, "W1": 309.54, "W2": 171.787, "W3": 31.674, "W4": 8.363}  # Jy (2MASS, Jarrett+2011)
RSUN_AU = 0.00465047


def main():
    des, teff, rstar, mstar, period = sys.argv[1], *map(float, sys.argv[2:6])
    ph = next(r for r in csv.DictReader(open(ROOT / "outputs/crossmatch/dd1_support_cache/allwise.csv")) if r["designation"] == des)
    sc = next(r for r in csv.DictReader(open(ROOT / "outputs/crossmatch/dd1_candidates_scored.csv")) if r["designation"] == des)
    m = {b: float(ph[f"{b}mag"]) for b in ("J", "H", "K", "W1", "W2", "W3", "W4")}
    e = {b: float(ph[f"e_{b}mag"]) for b in ("J", "H", "K", "W1", "W2", "W3", "W4")}
    w12, w23 = m["W1"] - m["W2"], m["W2"] - m["W3"]
    s12, s23 = np.hypot(e["W1"], e["W2"]), np.hypot(e["W2"], e["W3"])
    # Koenig & Leisawitz 2014, sec. 3 (after their contaminant cuts)
    if w12 > 1.0 and w23 > 2.0:
        kl = "Class I"
    elif (w12 - s12) > 0.25 and (w23 - s23) > 1.0:
        kl = "Class II"
    else:
        kl = "not YSO by WISE colours"
    # spectral index from K to W4 in lambda F_lambda
    lam = np.array([LAM[b] for b in ("K", "W1", "W2", "W3", "W4")])
    fnu = np.array([F0[b] * 10 ** (-0.4 * m[b]) for b in ("K", "W1", "W2", "W3", "W4")])
    lfl = fnu / lam  # proportional to lambda F_lambda
    alpha = float(np.polyfit(np.log10(lam), np.log10(lfl), 1)[0])
    cls_alpha = ("Class I" if alpha > 0.3 else "flat" if alpha > -0.3 else "Class II" if alpha > -1.6
                 else "Class III / transitional" if alpha > -2.5 else "photosphere")
    G, bprp, dist = float(sc["Gmag"]), float(sc["bp_rp"]), float(sc["gaia_dist"])
    # intrinsic BP-RP vs Teff (Pecaut & Mamajek dwarf table, coarse)
    tt = np.array([3500, 3800, 4100, 4400, 4600, 5000, 5300, 5800, 6200])
    bb = np.array([2.26, 1.95, 1.62, 1.40, 1.25, 1.04, 0.92, 0.78, 0.66])
    bprp0 = float(np.interp(teff, tt, bb))
    ebr = max(bprp - bprp0, 0.0)
    av = 2.0 * ebr  # A_V ~ 2 E(BP-RP) (rough)
    ag = 0.83 * av
    MG = G - 5 * np.log10(dist / 10) - ag
    # dwarf M_G vs Teff (Pecaut & Mamajek, coarse)
    mg = np.array([10.3, 9.0, 7.7, 6.9, 6.4, 5.8, 5.3, 4.5, 3.9])
    MG_ms = float(np.interp(teff, tt, mg))
    a_au = (mstar * (period / 365.25) ** 2) ** (1 / 3)
    a_rstar = a_au / (rstar * RSUN_AU)
    tdust = teff * np.sqrt(1 / (2 * a_rstar))
    out = {"designation": des, "colours": {"W1-W2": round(w12, 3), "W2-W3": round(w23, 3), "W3-W4": round(m["W3"] - m["W4"], 3),
                                          "H-K": round(m["H"] - m["K"], 3), "K-W1": round(m["K"] - m["W1"], 3)},
           "koenig_leisawitz_2014": kl, "alpha_K_W4": round(alpha, 2), "alpha_class": cls_alpha,
           "bp_rp": bprp, "bp_rp0_assumed": bprp0, "AV_rough": round(av, 2), "M_G": round(float(MG), 2),
           "M_G_main_sequence": MG_ms, "dM_G_above_MS": round(MG_ms - float(MG), 2),
           "corotation_au": round(a_au, 4), "corotation_rstar": round(a_rstar, 2), "T_dust_corotation_K": round(float(tdust)),
           "f_ir_cached": float(sc["f_ir"]), "fixed_f_ir_cached": float(sc["fixed_f_ir"]),
           "simbad_otype_cached": sc.get("simbad_otype"), "tic_lumclass": sc.get("tic_lumclass")}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
