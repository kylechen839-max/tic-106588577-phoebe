#!/usr/bin/env python3
"""BLS period search on TESS FFI light curves (AWS cubes) for objects with no catalogue period.

For each row of a CSV (designation,ra,dec,tic[,...]): cut up to --max-sectors FFI stamps, do the
fetch_tess_ffi_s3 photometry and scattered-light cut, divide out a per-orbit quadratic (no phase
model, so nothing periodic is removed except on > orbit timescales), then run BLS over 0.3-15 d.
Writes <outdir>/<des>/ (cutouts, raw LC, fold plot) and <outdir>/bls_results.csv.

  .venv-tess/bin/python src/checks/bls_ffi.py objects.csv outputs/candidates/tierB_tess/noperiod
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

rows_in = pd.read_csv(sys.argv[1])
outdir = Path(sys.argv[2])
maxsec = int(sys.argv[3]) if len(sys.argv) > 3 else 2
outdir.mkdir(parents=True, exist_ok=True)
periods = np.exp(np.linspace(np.log(0.3), np.log(15), 30000))
durations = [0.04, 0.08, 0.15, 0.3]
res_rows = []
for r in rows_in.itertuples():
    d = outdir / r.designation
    try:
        files = ffi.cutouts(r.ra, r.dec, str(int(r.tic)) if pd.notna(r.tic) else None, d / "ffi_cutouts", maxsec)
    except Exception as exc:
        print(r.designation, "cutout failed", exc)
        continue
    T, F, S = [], [], []
    for sec, fn in sorted(files.items())[:maxsec]:
        try:
            t, f, bg, *_ = ffi.photometry(fn, r.ra, r.dec)
            keep, trend = ffi.detrend(t, f, bg, 1.0, nharm=0)
        except Exception as exc:
            print(r.designation, sec, "photometry failed", exc)
            continue
        T.append(t[keep]); F.append((f / trend)[keep]); S.append(np.full(keep.sum(), sec))
    if not T:
        res_rows.append(dict(designation=r.designation, n=0))
        continue
    t, f, s = map(np.concatenate, (T, F, S))
    g = np.isfinite(f)
    t, f, s = t[g], f[g], s[g]
    np.savetxt(d / f"{r.designation}_raw.txt", np.c_[t, f, np.full_like(f, np.nan), s], header="time flux flux_err sector")
    rms = 1.4826 * np.median(np.abs(f - np.median(f)))
    bls = BoxLeastSquares(t, f)
    pw = bls.power(periods[periods < 0.5 * np.ptp(t)], durations, objective="snr")
    i = int(np.argmax(pw.power))
    P, t0, dep, dur, snr = pw.period[i], pw.transit_time[i], pw.depth[i], pw.duration[i], pw.depth_snr[i]
    st = bls.compute_stats(P, dur, t0)
    odd, even = st["depth_odd"], st["depth_even"]
    # secondary at phase 0.5 vs half period: if odd and even depths agree, 2P is also possible
    res_rows.append(dict(designation=r.designation, n=len(t), sectors=",".join(map(str, sorted(set(s.astype(int))))),
                         rms=round(rms, 5), P=round(P, 5), t0=round(t0, 4), depth=round(dep, 5), dur=dur,
                         snr=round(snr, 1), odd=round(odd[0], 5), even=round(even[0], 5),
                         odd_even_sig=round(abs(odd[0] - even[0]) / np.hypot(odd[1], even[1]), 1)))
    fig, ax = plt.subplots(2, 1, figsize=(9, 6))
    for sec in np.unique(s):
        m = s == sec
        ax[0].plot(t[m], f[m], ".", ms=1.5, label=f"S{int(sec)}")
        ax[1].plot(((t[m] - t0) / P + 0.25) % 1 - 0.25, f[m], ".", ms=1.5)
    ax[0].legend(fontsize=7)
    ax[1].set_title(f"BLS P={P:.5f} d, depth {dep*100:.2f}%, SNR {snr:.1f}, rms {rms*100:.2f}%", fontsize=9)
    ax[1].set_xlabel("phase")
    fig.suptitle(r.designation)
    fig.tight_layout()
    fig.savefig(d / f"{r.designation}_bls.png", dpi=100)
    plt.close(fig)
    print(res_rows[-1], flush=True)
pd.DataFrame(res_rows).to_csv(outdir / "bls_results.csv", index=False)
