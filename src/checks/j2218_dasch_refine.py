"""J2218: refined DASCH timing (per-plate-series zero points, finer windows) and LTTE refit.

Builds on j2218_dasch_timing.py / j2218_ltte_fit.py (branch claude/project-thread-5mcxzc). Changes:
  1. each plate series' weighted-median magnitude is subtracted (series with < 15 points dropped), removing
     zero-point offsets between telescopes;
  2. windows are tried at two widths (decade edges as before, and ~5 yr), with 300 bootstrap draws;
  3. the circular LTTE grid is refit to DASCH + ASAS-SN + Gaia + TESS for each window set.
Needs, relative to --root (default: repo root): ffi/ raw + ephemeris, timing/dasch csv, timing/oc_points json.
Output: timing/J221843.61+544715.0_dasch_refined.json
"""
import argparse, json
import numpy as np, pandas as pd
from scipy.optimize import least_squares

ap = argparse.ArgumentParser(); ap.add_argument("--root", default="."); a = ap.parse_args()
C = f"{a.root}/outputs/candidates/J221843.61+544715.0/"
d = np.loadtxt(C + "ffi/J221843.61+544715.0_raw.txt"); eph = json.load(open(C + "ffi/J221843.61+544715.0_ephemeris.json"))
P0, T0 = eph["period"], eph["t0"]
ref = np.isin(d[:, 3], [56, 57, 76, 77, 83, 84]); ph = ((d[ref, 0] - T0) / P0) % 1
e = np.arange(0, 1.0001, 0.005); x = 0.5 * (e[1:] + e[:-1]); idx = np.digitize(ph, e) - 1
y = np.array([np.median(d[ref, 1][idx == i]) for i in range(len(x))]); y /= np.median(y)
tm = -2.5 * np.log10(y); TX = np.r_[x - 1, x, x + 1]; TY = np.r_[tm, tm, tm]

g = pd.read_csv(C + "timing/dasch_dr7_apass_422184361544715.csv")
g = g[np.isfinite(g.magcal_magdep) & (g.magcal_local_rms < 0.25) & (g.time_accuracy_days < 0.01)
      & (g.limiting_mag_local > 13.0) & (g.reject_flag == 0)].copy()
med = g.magcal_magdep.median(); mad = 1.4826 * np.median(np.abs(g.magcal_magdep - med))
g = g[np.abs(g.magcal_magdep - med) < 4 * mad].copy()
g["s"] = np.maximum(g.magcal_local_rms, 0.08)


def prepare(series_zp):
    h = g.copy()
    if series_zp:
        cnt = h.series.value_counts(); h = h[h.series.isin(cnt[cnt >= 15].index)].copy()
        # out-of-eclipse points define each series' zero point (eclipses are ~15% of phase)
        oot = np.interp(((h.date_jd - 2457000 - T0) / P0) % 1, TX, TY) < 0.02
        zp = h[oot].groupby("series").magcal_magdep.median()
        h["m"] = h.magcal_magdep - h.series.map(zp)
        h = h[np.isfinite(h.m)]
    else:
        h["m"] = h.magcal_magdep
    t = h.date_jd.values - 2457000
    return h.m.values, h.s.values, t, 2000 + (h.date_jd.values - 2451544.5) / 365.25


shs = np.arange(-0.15, 0.15, 0.0005)


def prof(mm_, ww, ph):
    mm = mm_ - np.sum(ww * mm_) / ww.sum(); out = []
    for sh in shs:
        T = np.interp((ph - sh) % 1, TX, TY); Tc = T - np.sum(ww * T) / ww.sum()
        b = np.sum(ww * Tc * mm) / np.sum(ww * Tc * Tc); out.append(b * np.sum(ww * Tc * mm) if b > 0 else 0)
    return np.array(out)


def windows(m, s, t, yr, edges, nboot=300):
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (yr >= lo) & (yr < hi); n = int(sel.sum())
        if n < 40: continue
        mm, ss, tt = m[sel], s[sel], t[sel]; ww = 1 / ss**2; ph = ((tt - T0) / P0) % 1
        p = prof(mm, ww, ph); i = p.argmax()
        rng = np.random.default_rng(1); bs = []
        for _ in range(nboot):
            j = rng.integers(0, n, n); bs.append(shs[prof(mm[j], ww[j], ph[j]).argmax()])
        bs = np.array(bs) * P0 * 1440; l16, h84 = np.percentile(bs, [16, 84])
        rows.append(dict(lo=lo, hi=hi, yr=float(yr[sel].mean()), btjd=float(tt.mean()), n=n, dchi2=float(p[i]),
                         oc=float(shs[i] * P0 * 1440), lo16=float(l16), hi84=float(h84)))
    return rows


pts = json.load(open(C + "timing/J221843.61+544715.0_oc_points.json"))["points"]
yrf = lambda bt: 2000 + (np.asarray(bt) + 2457000 - 2451544.5) / 365.25


def ltte(rows, label):
    T, Y, E = [], [], []
    for r in rows:
        if r["dchi2"] >= 9:
            T.append(r["btjd"]); Y.append(r["oc"]); E.append(max(0.5 * (r["hi84"] - r["lo16"]), 10))
    nd = len(T)
    T.append(581); Y.append(50.5); E.append(5.8)                       # ASAS-SN 2016.6
    for p in pts: T.append(p["btjd"]); Y.append(p["oc_min"]); E.append(p["err_min"])
    T, Y, E = map(np.array, (T, Y, E)); ty = yrf(T)
    lin = lambda p: p[0] + p[1] * (ty - 2000)
    circ = lambda p: lin(p) + p[2] * np.sin(2 * np.pi * (ty - p[4]) / p[3])
    r = least_squares(lambda p: (lin(p) - Y) / E, [0, 0]); out = dict(n=len(Y), n_dasch=nd, linear_chi2=float(np.sum(r.fun**2)))
    scan = []
    for P3 in np.arange(12, 120, 1.0):
        best = None
        for phs in np.linspace(0, P3, 12, endpoint=False):
            rr = least_squares(lambda p: (circ(p) - Y) / E, [0, 0, 40, P3, 1950 + phs],
                               bounds=([-500, -50, 0, P3 - 0.5, 1800], [500, 50, 500, P3 + 0.5, 2200]))
            c = float(np.sum(rr.fun**2))
            if best is None or c < best[0]: best = (c, rr.x.tolist())
        scan.append((float(P3), *best))
    scan.sort(key=lambda z: z[1])
    out["circ_best"] = dict(P3=scan[0][2][3], chi2=scan[0][1], A_min=scan[0][2][2], params=scan[0][2])
    sc = sorted(scan); out["profile"] = [(p, c) for p, c, _ in sc]
    # local minima of the P3 profile
    cs = np.array([c for _, c, _ in sc]); ps = np.array([p for p, _, _ in sc])
    mins = [i for i in range(1, len(cs) - 1) if cs[i] <= cs[i - 1] and cs[i] <= cs[i + 1]]
    out["minima"] = sorted([(float(ps[i]), float(cs[i])) for i in mins], key=lambda z: z[1])[:4]
    print(f"{label}: {len(Y)} pts ({nd} DASCH); linear chi2 {out['linear_chi2']:.1f}; best circ P3 {out['circ_best']['P3']:.1f} yr "
          f"A {out['circ_best']['A_min']:.0f} min chi2 {out['circ_best']['chi2']:.1f} (dof {len(Y)-5}); minima {out['minima']}")
    return out


res = {}
for zp in (False, True):
    m, s, t, yr = prepare(zp)
    print(f"series zero points {'ON ' if zp else 'OFF'}: N {len(m)}, rms {np.std(m):.3f} mag")
    for lab, edges in (("decade", [1885, 1905, 1915, 1925, 1930, 1935, 1940, 1945, 1950, 1960, 1975, 1991]),
                       ("5yr", list(range(1885, 1995, 5)))):
        rows = windows(m, s, t, yr, edges)
        for r in rows:
            flag = "" if r["dchi2"] >= 9 else "  (dropped, dchi2<9)"
            print(f"   {lab:6s} {r['lo']}-{r['hi']} n={r['n']:4d} <yr>={r['yr']:.1f} dchi2={r['dchi2']:5.1f} "
                  f"O-C {r['oc']:+5.0f} [{r['lo16']:+.0f},{r['hi84']:+.0f}]{flag}")
        key = f"{'zp' if zp else 'raw'}_{lab}"
        res[key] = dict(rms=float(np.std(m)), N=len(m), windows=rows, fit=ltte(rows, key))
json.dump(res, open(C + "timing/J221843.61+544715.0_dasch_refined.json", "w"), indent=1)
