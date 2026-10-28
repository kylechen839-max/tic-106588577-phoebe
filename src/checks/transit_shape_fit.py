"""Primary-eclipse shape fit with a dark-companion (Mandel-Agol / batman) model, for EBs with no
detected secondary. Fits k = R2/R1, a/R1, b, a time offset and a 2-harmonic out-of-eclipse wave;
reports a chi2 profile in b (grazing vs central), and turns a secondary-depth upper limit into an
upper limit on the surface-brightness ratio S2/S1. For a circular orbit the overlap area is the same
at both conjunctions, so S2/S1 ~ secondary depth / primary depth (edge-intensity weighted; exact for
a grazing eclipse up to limb darkening, approximate for a central one).

Usage: transit_shape_fit.py TIMESERIES P T0 OUTPREFIX [sec_limit_ppt] [u1 u2]
Depths are as measured in the given light curve: unknown dilution makes k a lower limit.
"""
import sys
import numpy as np
import batman
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fn, P, T0, outp = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
seclim = float(sys.argv[5]) if len(sys.argv) > 5 else 2.0
u = [float(sys.argv[6]), float(sys.argv[7])] if len(sys.argv) > 7 else [0.30, 0.25]
d = np.loadtxt(fn)
t, f = d[:, 0], d[:, 1]
ph = ((t - T0) / P + 0.5) % 1 - 0.5
pars = batman.TransitParams()
pars.per, pars.t0, pars.ecc, pars.w, pars.limb_dark, pars.u = P, 0.0, 0.0, 90.0, "quadratic", u
tt = ph * P
pars.rp, pars.a, pars.inc = 0.15, 6.0, 85.0
m = batman.TransitModel(pars, tt)

def model(v, b_fixed=None):
    k, a, b, dt, c0, c1, s1, c2, s2 = v
    if b_fixed is not None:
        b = b_fixed
    pars.rp, pars.a = k, a
    pars.inc = np.degrees(np.arccos(np.clip(b / a, 0, 1)))
    pars.t0 = dt
    w = c0 + c1 * np.cos(2 * np.pi * ph) + s1 * np.sin(2 * np.pi * ph) + c2 * np.cos(4 * np.pi * ph) + s2 * np.sin(4 * np.pi * ph)
    return m.light_curve(pars) * w

oot = np.abs(ph) > 0.06
sig = 1.4826 * np.median(np.abs(np.diff(f[oot]))) / np.sqrt(2)

def resid(v, b_fixed=None):
    return (f - model(v, b_fixed)) / sig

v0 = [0.15, 6.0, 0.5, 0.0, 1.0, 0, 0, 0, 0]
lo = [0.01, 1.5, 0.0, -0.02, 0.9, -0.05, -0.05, -0.05, -0.05]
hi = [1.0, 30.0, 1.5, 0.02, 1.1, 0.05, 0.05, 0.05, 0.05]
best = None
for b0 in (0.1, 0.5, 0.8, 1.0, 1.2):
    v0[2] = b0
    r = least_squares(resid, v0, bounds=(lo, hi))
    if best is None or r.cost < best.cost:
        best = r
v = best.x
chi2 = 2 * best.cost
n = len(f)
res = f - model(v)
# red-noise beta from 30-min bins of residuals
nb = max(1, int(round(30 * 60 / np.median(np.diff(t) * 86400))))
rb = [res[i:i + nb].mean() for i in range(0, n - nb, nb)]
beta = max(1.0, np.std(rb) / (np.std(res) / np.sqrt(nb)))
J = best.jac
cov = np.linalg.pinv(J.T @ J) * beta ** 2
err = np.sqrt(np.diag(cov))
names = ["k", "a/R1", "b", "dt", "c0", "c1", "s1", "c2", "s2"]
lines = [f"# {fn}  P={P} T0={T0}  LD u={u}  N={n}  sigma_pt={sig*1e3:.2f} ppt  beta={beta:.2f}",
         f"chi2 {chi2:.1f}  chi2_red {chi2/(n-9):.3f}"]
for nm, x, e in zip(names, v, err):
    lines.append(f"{nm} {x:.5f} +/- {e:.5f}")
k, a, b = v[:3]
inc = np.degrees(np.arccos(b / a))
r1, r2 = 1 / a, k / a
lines.append(f"incl {inc:.2f} deg; r1=R1/a {r1:.4f}; r2=R2/a {r2:.4f}; r1+r2 {r1+r2:.4f}")
lines.append(f"stellar density (R1, circular, M2<<M1) rho1 = {3*np.pi*a**3/(6.674e-8*(P*86400)**2):.3f} g/cm3")
# b profile
prof = []
for bf in np.arange(0.0, 1.31, 0.1):
    vv = v.copy(); vv[2] = bf
    rr = least_squares(lambda x: resid(np.r_[x[:2], bf, x[2:]]), np.r_[vv[:2], vv[3:]],
                       bounds=(np.r_[lo[:2], lo[3:]], np.r_[hi[:2], hi[3:]]))
    prof.append((bf, 2 * rr.cost, rr.x[0], rr.x[1]))
lines.append("b profile (b, dchi2/beta^2, k, a/R1):")
c0 = min(p[1] for p in prof)
for bf, c, kk, aa in prof:
    lines.append(f"  {bf:.1f} {(c - c0)/beta**2:8.1f} {kk:.4f} {aa:.2f}")
ok = [p for p in prof if (p[1] - c0) / beta ** 2 < 4]
kmin, kmax = min(p[2] for p in ok), max(p[2] for p in ok)
lines.append(f"2-sigma range: b {min(p[0] for p in ok):.1f}-{max(p[0] for p in ok):.1f}, k {kmin:.3f}-{kmax:.3f}")
sl = seclim / 1e3
pdepth = 1 - model(v)[np.abs(ph) < 0.002].min() / v[4]
lines.append(f"model primary depth {pdepth*100:.2f}%; secondary < {seclim} ppt -> S2/S1 < {sl/pdepth:.3f} (TESS band)")
open(outp + ".txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
o = np.argsort(ph)
sel = np.abs(ph[o]) < 0.12
ax[0].plot(ph[o][sel], f[o][sel], ".", ms=1.5, color="0.6")
edges = np.linspace(-0.12, 0.12, 97); idx = np.digitize(ph, edges) - 1
cb = [(edges[i] + edges[i + 1]) / 2 for i in range(96)]
ax[0].plot(cb, [np.median(f[idx == i]) for i in range(96)], "ko", ms=3)
ax[0].plot(ph[o][sel], model(v)[o][sel], "C3-", lw=1.2, label=f"k={k:.3f}, b={b:.2f}, a/R1={a:.2f}")
ax[0].legend(); ax[0].set_ylabel("flux")
ax[1].plot(cb, [np.median(res[idx == i]) * 1e3 for i in range(96)], "ko", ms=3); ax[1].axhline(0, color="C3")
ax[1].set_ylabel("resid (ppt)"); ax[1].set_xlabel("phase")
fig.tight_layout(); fig.savefig(outp + ".png", dpi=120)
