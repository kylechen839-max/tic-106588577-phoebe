"""Model-independent primary vs secondary eclipse timing per TESS sector.

Builds an empirical eclipse template (primary and secondary separately) by folding the
short-cadence sectors, then for every sector shifts the template in time and fits
flux = a + b * template(t - dt). The template is boxcar-smoothed to each sector's
exposure time, so 30-min FFI sectors are compared with a 30-min-smeared template.

Why: apsidal motion moves the primary and secondary in opposite directions;
a light-travel-time effect (third body) or a period change moves both the same way.

Usage: eclipse_split_timing.py RAW.txt EPHEMERIS.json OUT_PREFIX [--ref-sectors 56,57,...]
Writes OUT_PREFIX_split_timing.json and OUT_PREFIX_split_timing.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def bin_profile(ph, f, centre, half, step):
    edges = np.arange(centre - half, centre + half + step / 2, step)
    idx = np.digitize(ph, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    x = 0.5 * (edges[1:] + edges[:-1])
    y = np.array([np.median(f[ok][idx[ok] == i]) if np.any(idx[ok] == i) else np.nan for i in range(len(x))])
    g = np.isfinite(y)
    return x[g], y[g]


def smear(x, y, width):
    """Boxcar-average the template over an exposure of `width` (phase units)."""
    if width <= 0:
        return y
    offs = np.linspace(-width / 2, width / 2, 15)
    return np.mean([np.interp(x + o, x, y) for o in offs], axis=0)


def scan(tph, fl, sg, tx, ty, P, shifts_min):
    chi = []
    for dt in shifts_min:
        m = np.interp(tph - dt / 1440.0 / P, tx, ty, left=np.nan, right=np.nan)
        ok = np.isfinite(m)
        A = np.vstack([np.ones(ok.sum()), m[ok]]).T / sg[ok, None]
        coef, *_ = np.linalg.lstsq(A, fl[ok] / sg[ok], rcond=None)
        r = (fl[ok] - np.vstack([np.ones(ok.sum()), m[ok]]).T @ coef) / sg[ok]
        chi.append((r @ r, ok.sum()))
    return np.array([c[0] for c in chi]), np.array([c[1] for c in chi])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw")
    ap.add_argument("ephem")
    ap.add_argument("out")
    ap.add_argument("--ref-sectors", default="56,57,76,77,83,84")
    ap.add_argument("--half", type=float, default=0.075, help="half window in phase around each eclipse")
    args = ap.parse_args()

    d = np.loadtxt(args.raw)
    d = d[np.isfinite(d).all(axis=1)]
    eph = json.loads(Path(args.ephem).read_text())
    P, t0 = eph["period"], eph["t0"]
    ph2 = eph["eclipses"]["phase2"]
    t, f, s, sec = d.T
    ph = ((t - t0) / P + 0.5) % 1 - 0.5  # primary at 0
    ref = np.isin(sec, [float(x) for x in args.ref_sectors.split(",")])

    templates = {}
    for name, c in (("primary", 0.0), ("secondary", ph2)):
        p = ph if c == 0 else (((t - t0) / P - c + 0.5) % 1 - 0.5) + c
        templates[name] = bin_profile(p[ref], f[ref], c, args.half + 0.02, 0.001)

    shifts = np.arange(-25, 25.01, 0.1)
    out = {"period": P, "t0": t0, "phase2": ph2, "ref_sectors": args.ref_sectors, "sectors": {}}
    for sv in np.unique(sec):
        m = sec == sv
        cad = np.median(np.diff(t[m]))
        res = {"t_mid": float(np.median(t[m])), "cadence_min": float(cad * 1440)}
        for name, c in (("primary", 0.0), ("secondary", ph2)):
            p = ph[m] if c == 0 else (((t[m] - t0) / P - c + 0.5) % 1 - 0.5) + c
            w = np.abs(p - c) < args.half
            tx, ty = templates[name]
            ty_s = smear(tx, ty, cad / P)
            chi, n = scan(p[w], f[m][w], s[m][w], tx, ty_s, P, shifts)
            i = int(np.argmin(chi))
            red = chi[i] / max(n[i] - 3, 1)
            scale = max(red, 1.0)
            within = shifts[(chi - chi[i]) / scale <= 1.0]
            res[name] = {
                "dt_min": float(shifts[i]),
                "err_min": float(0.5 * (within.max() - within.min())),
                "chi2_red": float(red),
                "n_points": int(n[i]),
                "n_eclipses": int(len(np.unique(np.round((t[m][w] - t0) / P - c)))),
            }
        res["prim_minus_sec_min"] = res["primary"]["dt_min"] - res["secondary"]["dt_min"]
        res["prim_minus_sec_err"] = float(np.hypot(res["primary"]["err_min"], res["secondary"]["err_min"]))
        out["sectors"][str(int(sv))] = res
        print(f"S{int(sv):3d} cad {cad*1440:5.1f} min  prim {res['primary']['dt_min']:+6.2f}±{res['primary']['err_min']:.2f}"
              f"  sec {res['secondary']['dt_min']:+6.2f}±{res['secondary']['err_min']:.2f}"
              f"  prim-sec {res['prim_minus_sec_min']:+6.2f}±{res['prim_minus_sec_err']:.2f} min")

    Path(args.out + "_split_timing.json").write_text(json.dumps(out, indent=1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    for name, col, off in (("primary", "C0", -8), ("secondary", "C3", 8)):
        tm = np.array([v["t_mid"] for v in out["sectors"].values()])
        y = np.array([v[name]["dt_min"] for v in out["sectors"].values()])
        e = np.array([v[name]["err_min"] for v in out["sectors"].values()])
        ax.errorbar(tm + off, y, e, fmt="o", color=col, label=name)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("BTJD")
    ax.set_ylabel("O - C (min), vs template from S" + args.ref_sectors.replace(",", "/S"))
    ax.legend()
    ax.set_title(Path(args.out).name + ": primary vs secondary timing")
    fig.tight_layout()
    fig.savefig(args.out + "_split_timing.png", dpi=130)


if __name__ == "__main__":
    main()
