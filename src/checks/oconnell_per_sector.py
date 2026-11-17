"""O'Connell effect and out-of-eclipse asymmetry per TESS sector (model-independent).

For each sector, fold on the ephemeris and measure
  dmax = <flux, phase 0.20-0.30> - <flux, phase 0.70-0.80>   (Max I minus Max II, quadrature)
and the odd (sin) Fourier terms of the out-of-eclipse light (|phase - eclipse| > 0.1),
  f = c0 + a1 cos + b1 sin + a2 cos2 + b2 sin2.
A hot spot from accretion is fixed in the binary frame (stable b1); cool starspots migrate
and change in size, so b1 varies from sector to sector and year to year.

Usage: oconnell_per_sector.py RAW.txt EPHEMERIS.json OUT.json
"""
import json, sys
import numpy as np

raw, ephf, out = sys.argv[1:4]
d = np.loadtxt(raw); d = d[np.isfinite(d).all(axis=1)]
eph = json.load(open(ephf)); P, t0 = eph["period"], eph["t0"]
res = {}
for s in np.unique(d[:, 3]):
    m = d[:, 3] == s
    ph = ((d[m, 0] - t0) / P) % 1; f = d[m, 1] / np.median(d[m, 1]); e = d[m, 2]
    q1 = (ph > 0.2) & (ph < 0.3); q2 = (ph > 0.7) & (ph < 0.8)
    dmax = f[q1].mean() - f[q2].mean(); derr = np.hypot(f[q1].std() / np.sqrt(q1.sum()), f[q2].std() / np.sqrt(q2.sum()))
    oot = (np.abs(ph - 0.0) > 0.1) & (np.abs(ph - 1.0) > 0.1) & (np.abs(ph - 0.5) > 0.1)
    x = 2 * np.pi * ph[oot]
    A = np.vstack([np.ones_like(x), np.cos(x), np.sin(x), np.cos(2 * x), np.sin(2 * x)]).T
    c, *_ = np.linalg.lstsq(A, f[oot], rcond=None)
    r = f[oot] - A @ c; cov = np.linalg.inv(A.T @ A) * r.var()
    err = np.sqrt(np.diag(cov))
    res[str(int(s))] = dict(t_mid=float(np.median(d[m, 0])), dmax_ppt=1e3 * dmax, dmax_err_ppt=1e3 * derr,
                            b1_ppt=1e3 * c[2], b1_err=1e3 * err[2], b2_ppt=1e3 * c[4], b2_err=1e3 * err[4],
                            a1_ppt=1e3 * c[1], a2_ppt=1e3 * c[3])
    v = res[str(int(s))]
    print(f"S{int(s):3d}: MaxI-MaxII {v['dmax_ppt']:+6.2f}±{v['dmax_err_ppt']:.2f} ppt   sin1 {v['b1_ppt']:+6.2f}±{v['b1_err']:.2f}  "
          f"sin2 {v['b2_ppt']:+5.2f}±{v['b2_err']:.2f}  cos1 {v['a1_ppt']:+6.2f}  cos2 {v['a2_ppt']:+6.2f} ppt")
json.dump(res, open(out, "w"), indent=1)
