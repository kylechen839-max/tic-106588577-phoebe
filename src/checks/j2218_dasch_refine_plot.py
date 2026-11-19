"""Plot of the refined J2218 O-C (series zero points, decade windows) with the best circular LTTE orbits at
fixed P3 = 25, 36 and 62 yr. Reads timing/J221843.61+544715.0_dasch_refined.json and timing/..._oc_points.json
(--root: a checkout of claude/project-thread-5mcxzc for the oc_points file). Writes timing/..._oc_refined.png."""
import argparse, json
import numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("--root", default="."); ap.add_argument("--out", default=".")
a = ap.parse_args()
TD = "outputs/candidates/J221843.61+544715.0/timing/"
R = json.load(open(f"{a.out}/{TD}J221843.61+544715.0_dasch_refined.json"))
pts = json.load(open(f"{a.root}/{TD}J221843.61+544715.0_oc_points.json"))["points"]
yrf = lambda bt: 2000 + (np.asarray(bt) + 2457000 - 2451544.5) / 365.25
fig, ax = plt.subplots(figsize=(9.5, 4.8))
w = R["zp_decade"]["windows"]
T, Y, E, lab = [], [], [], []
for r in w:
    ok = r["dchi2"] >= 9
    ax.errorbar(r["yr"], r["oc"], yerr=[[r["oc"] - r["lo16"]], [r["hi84"] - r["oc"]]], fmt="o" if ok else "o",
                mfc="C0" if ok else "none", color="C0", alpha=1 if ok else 0.35, capsize=2)
    if ok: T.append(r["btjd"]); Y.append(r["oc"]); E.append(max(0.5 * (r["hi84"] - r["lo16"]), 10))
T.append(581); Y.append(50.5); E.append(5.8)
for p in pts: T.append(p["btjd"]); Y.append(p["oc_min"]); E.append(p["err_min"])
T, Y, E = map(np.array, (T, Y, E)); ty = yrf(T)
ax.errorbar(yrf(581), 50.5, 5.8, fmt="s", color="C2", label="ASAS-SN")
for p in pts: ax.errorbar(yrf(p["btjd"]), p["oc_min"], p["err_min"], fmt="D", ms=4, color="C3")
ax.plot([], [], "D", color="C3", label="Gaia, TESS")
xx = np.linspace(1890, 2040, 800)
for P3, c in ((25, "k"), (36, "C1"), (62, "C4")):
    best = None
    for phs in np.linspace(0, P3, 12, endpoint=False):
        f = lambda p, t: p[0] + p[1] * (t - 2000) + p[2] * np.sin(2 * np.pi * (t - p[3]) / P3)
        r = least_squares(lambda p: (f(p, ty) - Y) / E, [0, 0, 40, 1950 + phs], bounds=([-500, -50, 0, 1800], [500, 50, 500, 2200]))
        if best is None or r.cost < best[0]: best = (r.cost, r.x, f)
    ax.plot(xx, best[2](best[1], xx), "-", color=c, lw=1.2,
            label=f"P₃ = {P3} yr, amplitude {best[1][2]:.0f} min, χ² {2 * best[0]:.1f}")
ax.axhline(0, color="0.7", lw=0.6); ax.set_ylim(-260, 260); ax.set_xlim(1890, 2040)
ax.set_xlabel("year"); ax.set_ylabel("O − C (minutes)")
ax.set_title("J2218 eclipse timing, 1890–2024 (plates: series zero points removed; open = too weak to time)", fontsize=10)
ax.legend(fontsize=8, loc="lower left"); fig.tight_layout()
fig.savefig(f"{a.out}/{TD}J221843.61+544715.0_oc_refined.png", dpi=140); print("saved")
