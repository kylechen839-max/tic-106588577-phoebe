#!/usr/bin/env python3
"""BLS period search on TESS FFI light curves (AWS cubes) for objects with no catalogue period.

For each row of a CSV (designation,ra,dec,tic[,...]): cut up to MAX_SECTORS FFI stamps, do the
fetch_tess_ffi_s3 photometry and scattered-light cut, divide out a per-orbit quadratic (no phase
model, so nothing periodic is removed except on > orbit timescales), drop isolated outliers, then run
BLS over 0.35-15 d with white-noise errors. A light curve already in <outdir>/<des>/ is reused.
Writes <outdir>/<des>/ (cutouts, raw LC, BLS plot) and <outdir>/bls_results.csv.

  .venv-tess/bin/python src/checks/bls_ffi.py objects.csv outputs/candidates/tierB_tess/noperiod [MAX_SECTORS]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
import fetch_tess_ffi_s3 as ffi  # noqa: E402

PERIODS = np.exp(np.linspace(np.log(0.35), np.log(15), 30000))
DURATIONS = [0.04, 0.08, 0.15, 0.3]


def cutout_lc(r, d, maxsec):
    files = ffi.cutouts(r.ra, r.dec, str(int(r.tic)) if pd.notna(r.tic) else None, d / "ffi_cutouts", maxsec)
    T, F, S = [], [], []
    for sec, fn in sorted(files.items())[:maxsec]:
        try:
            t, f, bg, *_ = ffi.photometry(fn, r.ra, r.dec)
            keep, trend = ffi.detrend(t, f, bg, 1.0, nharm=0)
        except Exception as exc:
            print(r.designation, sec, "photometry failed", exc)
            continue
        T.append(t[keep]); F.append((f / trend)[keep]); S.append(np.full(keep.sum(), sec))
    return T, F, S


def search(t, f):
    rms = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
    rm = np.array([np.median(f[max(0, i - 3):i + 4]) for i in range(len(f))])
    ok = np.abs(f - rm) < 5 * rms   # isolated outliers only; eclipses span many points
    bls = BoxLeastSquares(t[ok], f[ok], dy=np.full(ok.sum(), rms))
    pw = bls.power(PERIODS[PERIODS < 0.5 * np.ptp(t)], DURATIONS, objective="likelihood")
    i = int(np.argmax(pw.power))
    P, t0, dep, dur = pw.period[i], pw.transit_time[i], pw.depth[i], pw.duration[i]
    st = bls.compute_stats(P, dur, t0)
    odd, even = st["depth_odd"], st["depth_even"]
    return ok, rms, dict(P=round(P, 5), t0=round(t0, 4), depth=round(dep, 5), dur=dur,
                         snr=round(dep / pw.depth_err[i], 1),
                         odd_even_sig=round(abs(odd[0] - even[0]) / np.hypot(odd[1], even[1]), 1),
                         depth_phase05=round(st["depth_half"][0], 5))


def plot(name, d, t, f, s, res, rms):
    P, t0 = res["P"], res["t0"]
    fig, ax = plt.subplots(2, 1, figsize=(9, 6))
    for sec in np.unique(s):
        m = s == sec
        ax[0].plot(t[m], f[m], ".", ms=1.5, label=f"S{int(sec)}")
        ax[1].plot(((t[m] - t0) / P + 0.25) % 1 - 0.25, f[m], ".", ms=1.5)
    ax[0].legend(fontsize=7)
    ax[1].set_title(f"BLS P={P:.5f} d, depth {res['depth']*100:.2f}%, S/N {res['snr']:.0f}, rms {rms*100:.2f}%", fontsize=9)
    ax[1].set_xlabel("phase")
    fig.suptitle(name)
    fig.tight_layout()
    fig.savefig(d / f"{name}_bls.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    rows_in = pd.read_csv(sys.argv[1])
    outdir = Path(sys.argv[2])
    maxsec = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    rows = []
    for r in rows_in.itertuples():
        d = outdir / r.designation
        rawf = d / f"{r.designation}_raw.txt"
        if rawf.exists():
            a = np.genfromtxt(rawf, comments="#")
            T, F, S = [a[:, 0]], [a[:, 1]], [a[:, 3]]
        else:
            T, F, S = cutout_lc(r, d, maxsec)
        if not T:
            rows.append(dict(designation=r.designation, n=0))
            continue
        t, f, s = map(np.concatenate, (T, F, S))
        g = np.isfinite(f)
        t, f, s = t[g], f[g], s[g]
        np.savetxt(rawf, np.c_[t, f, np.full_like(f, np.nan), s], header="time flux flux_err sector")
        ok, rms, res = search(t, f)
        rows.append(dict(designation=r.designation, n=int(ok.sum()),
                         sectors=",".join(map(str, sorted(set(s.astype(int))))), rms=round(rms, 5), **res))
        plot(r.designation, d, t[ok], f[ok], s[ok], res, rms)
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(outdir / "bls_results.csv", index=False)
