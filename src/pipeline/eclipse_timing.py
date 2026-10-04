#!/usr/bin/env python3
"""Per-sector primary-eclipse timings (O-C) for a candidate with a fetched raw TESS light curve.

Fits a Gaussian dip to each sector's eclipses folded on a reference ephemeris,
then fits a linear ephemeris to the sector mid-times and compares with an
optional older ephemeris.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[2]


def dip(x, a, mu, sig, c):
    return c - a * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designation", required=True)
    ap.add_argument("--old-period", type=float, default=None)
    ap.add_argument("--old-t0", type=float, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--secondary", action="store_true",
                    help="time the secondary eclipse instead (uses eclipses.phase2/halfwidth2 from the ephemeris "
                         "JSON); a primary/secondary period difference measures apsidal motion")
    a = ap.parse_args()
    d = ROOT / "outputs" / "candidates" / a.designation
    t, f, e, s = np.loadtxt(d / f"{a.designation}_raw.txt").T
    eph = json.loads((d / f"{a.designation}_ephemeris.json").read_text())
    P, t0 = eph["period"], eph["t0"]
    ecl = eph.get("eclipses", {})
    hw = ecl.get("halfwidth1", 0.03)
    depth0 = eph.get("primary_depth", 0.05)
    if a.secondary:  # shift the reference epoch to the secondary and time it like a primary
        t0 = t0 + ecl["phase2"] * P
        hw, depth0 = ecl.get("halfwidth2", hw), ecl.get("depth2", depth0)
    rows = []
    for sec in np.unique(s):
        m = s == sec
        ph = ((t[m] - t0) / P + 0.5) % 1 - 0.5
        w = np.abs(ph) < 2 * hw
        if w.sum() < 20:
            continue
        x = ph[w] * P * 24
        try:
            popt, pcov = curve_fit(dip, x, f[m][w], p0=[depth0, 0, hw * P * 24 / 2, 1.0])
        except RuntimeError:
            continue
        n = np.round((np.median(t[m]) - t0) / P)
        rows.append((sec, n, t0 + n * P + popt[1] / 24, np.sqrt(pcov[1, 1]) / 24))
    rows = np.array(rows)
    sec, n, tm, err = rows.T
    A = np.vstack([np.ones_like(n), n]).T
    W = 1 / err**2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    T0, Pfit = cov @ (A.T @ (W * tm))
    res = {"designation": a.designation, "T0": T0, "P": Pfit, "e_T0": np.sqrt(cov[0, 0]), "e_P": np.sqrt(cov[1, 1]),
           "sectors": sec.tolist(), "cycles": n.tolist(), "mid_times": tm.tolist(), "errors_d": err.tolist(),
           "oc_new_h": ((tm - (T0 + n * Pfit)) * 24).tolist()}
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(tm, (tm - (T0 + n * Pfit)) * 24, yerr=err * 24, fmt="o", label=f"new P={Pfit:.7f} d")
    if a.old_period:
        nold = np.round((tm - a.old_t0) / a.old_period)
        oc_old = (tm - (a.old_t0 + nold * a.old_period)) * 24
        res["oc_old_h"] = oc_old.tolist()
        ax.errorbar(tm, oc_old, yerr=err * 24, fmt="s", label=f"old P={a.old_period:.7f} d")
    for x, y, lab in zip(tm, (tm - (T0 + n * Pfit)) * 24, sec):
        ax.annotate(f"S{int(lab)}", (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7)
    ax.axhline(0, color="0.5")
    ax.set_xlabel("BTJD"); ax.set_ylabel("O-C (hours)"); ax.legend(fontsize=8)
    which = "secondary" if a.secondary else "primary"
    ax.set_title(f"{a.designation} {which}-eclipse timings")
    fig.tight_layout()
    sfx = "_sec" if a.secondary else ""
    fig.savefig(a.out or d / f"{a.designation}_oc{sfx}.png", dpi=140)
    (d / f"{a.designation}_timing{sfx}.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps({k: res[k] for k in ("P", "e_P", "T0", "e_T0")}, indent=2))
    if "oc_old_h" in res:
        print("O-C vs old ephemeris (h):", dict(zip([int(x) for x in sec], np.round(res["oc_old_h"], 2))))


if __name__ == "__main__":
    main()
