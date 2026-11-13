"""Cleaned multi-sector TESS light curve for J013155.91-004225.4 = TIC 248986951.

Follows the project's convention (docs/PROJECT_CONTEXT.md): per sector SPOC, then TESS-SPOC,
then QLP; per-orbit median normalisation; outlier rejection that keeps the eclipses (residuals
from a phase-binned running median, clipped at 5 sigma).
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import lightkurve as lk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
NAME, TIC, P_VSX = "J013155.91-004225.4", 248986951, 0.351176
PREF = ["SPOC", "TESS-SPOC", "QLP"]

res = lk.search_lightcurve(f"TIC {TIC}")
tab = res.table.to_pandas()
tab["sector"] = tab["mission"].str.extract(r"(\d+)").astype(int)
rows, log = [], []
for sec in sorted(tab.sector.unique()):
    for auth in PREF:
        idx = tab.index[(tab.sector == sec) & (tab.author == auth)]
        if len(idx):
            break
    else:
        continue
    i = idx[np.argmin(tab.loc[idx, "exptime"].values)]  # shortest cadence for that author
    lc = res[int(i)].download(quality_bitmask="default")
    meta = lc.meta
    col = "pdcsap_flux" if auth != "QLP" else ("kspsap_flux" if "kspsap_flux" in lc.colnames else "sap_flux")
    lc = lc.select_flux(col).remove_nans()
    t, f, e = lc.time.value, lc.flux.value.astype(float), lc.flux_err.value.astype(float)
    # per-orbit normalisation: split at the largest gap
    gap = np.argmax(np.diff(t)) + 1 if len(t) > 1 and np.max(np.diff(t)) > 0.5 else len(t)
    orbit = np.where(np.arange(len(t)) < gap, 1, 2)
    for o in (1, 2):
        m = orbit == o
        if m.any():
            med = np.nanmedian(f[m]); f[m] /= med; e[m] /= med
    rows.append(pd.DataFrame(dict(time_btjd=t, flux=f, flux_err=e, sector=sec, orbit=orbit, author=auth,
                                  exptime_s=float(tab.loc[i, "exptime"]))))
    log.append(dict(sector=sec, author=auth, exptime_s=float(tab.loc[i, "exptime"]), flux_col=col, n_raw=len(t),
                    crowdsap=meta.get("CROWDSAP"), flfrcsap=meta.get("FLFRCSAP")))
    print(f"S{sec:03d} {auth:9s} {tab.loc[i, 'exptime']:5.0f}s {col:12s} n={len(t)} CROWDSAP={meta.get('CROWDSAP')}")

d = pd.concat(rows, ignore_index=True)

# period check on the 2-min / short-cadence data
ls = LombScargle(d.time_btjd, d.flux)
fr, pw = ls.autopower(minimum_frequency=0.5, maximum_frequency=20, samples_per_peak=20)
p_ls = 1 / fr[np.argmax(pw)]
print(f"Lomb-Scargle peak period {p_ls:.6f} d (x2 = {2 * p_ls:.6f}); VSX {P_VSX}")

# refine P near VSX by phase-dispersion minimisation over the full baseline
grid = np.linspace(P_VSX * (1 - 2e-4), P_VSX * (1 + 2e-4), 4001)
def disp(P):
    ph = (d.time_btjd.values / P) % 1
    b = np.floor(ph * 100).astype(int)
    s = pd.Series(d.flux.values).groupby(b).agg(["var", "count"])
    return np.nansum(s["var"] * (s["count"] - 1)) / (len(d) - len(s))
dd = np.array([disp(P) for P in grid]); P = grid[np.argmin(dd)]
print(f"PDM-refined period {P:.7f} d")

# eclipse-preserving outlier filter: residual from phase-binned median, 5 sigma MAD, iterated
keep = np.ones(len(d), bool)
for _ in range(3):
    ph = (d.time_btjd.values / P) % 1
    b = np.floor(ph * 200).astype(int)
    med = pd.Series(np.where(keep, d.flux, np.nan)).groupby(b).transform("median").values
    r = d.flux.values - med
    out = np.zeros(len(d), bool)
    for s in d.sector.unique():
        m = (d.sector == s).values
        sig = 1.4826 * np.nanmedian(np.abs(r[m & keep]))
        out |= m & (np.abs(r) > 5 * sig)
    keep = ~out
d["clipped"] = ~keep
clean = d[keep].copy()
t0 = clean.time_btjd.values[np.argmin(pd.Series(clean.flux.values).rolling(15, center=True, min_periods=5).median().values)]
clean["phase"] = (((clean.time_btjd - t0) / P + 0.25) % 1) - 0.25
d.to_csv(OUT / f"{NAME}_tess_all_with_flags.csv", index=False)
clean.to_csv(OUT / f"{NAME}_tess_clean.csv", index=False)
lg = pd.DataFrame(log)
lg["n_clean"] = [int(((clean.sector == s)).sum()) for s in lg.sector]
lg["rms_ppt"] = [1e3 * 1.4826 * np.median(np.abs(r[(d.sector == s).values & keep])) for s in lg.sector]
lg.to_csv(OUT / f"{NAME}_sector_log.csv", index=False)
print(lg.to_string(index=False))
print(f"t0 (min) = {t0:.5f} BTJD, P = {P:.7f} d, total clean points {len(clean)}, clipped {int((~keep).sum())}")

# figures
secs = sorted(clean.sector.unique())
fig, axs = plt.subplots(len(secs), 1, figsize=(11, 1.9 * len(secs)), squeeze=False)
for ax, s in zip(axs[:, 0], secs):
    a = d[d.sector == s]
    ax.plot(a.time_btjd[~a.clipped], a.flux[~a.clipped], ".", ms=1, c="k")
    ax.plot(a.time_btjd[a.clipped], a.flux[a.clipped], "x", ms=3, c="r")
    ax.set_ylabel(f"S{s}\n{a.author.iloc[0]}", fontsize=7); ax.tick_params(labelsize=7)
axs[-1, 0].set_xlabel("BTJD")
fig.suptitle(f"{NAME} (TIC {TIC}) per-orbit normalised; red = clipped", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_tess_sectors.png", dpi=110); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
cols = plt.cm.viridis(np.linspace(0, 0.9, len(secs)))
for c, s in zip(cols, secs):
    a = clean[clean.sector == s]
    axs[0].plot(a.phase, a.flux, ".", ms=1, color=c, alpha=0.4, label=f"S{s}")
bins = np.linspace(-0.25, 0.75, 101)
mid = 0.5 * (bins[1:] + bins[:-1])
bm = pd.Series(clean.flux.values).groupby(np.digitize(clean.phase, bins)).median()
for ax in axs:
    ax.plot(mid[bm.index.values[(bm.index > 0) & (bm.index <= 100)] - 1], bm[(bm.index > 0) & (bm.index <= 100)], "r-", lw=1)
axs[0].legend(markerscale=8, fontsize=7); axs[0].set_xlabel("phase"); axs[0].set_ylabel("normalised flux")
axs[1].set_xlabel("phase"); axs[1].set_title("100-bin median", fontsize=8)
fig.suptitle(f"{NAME} folded, P = {P:.7f} d, t0 = {t0:.4f} BTJD", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_tess_folded.png", dpi=110); plt.close(fig)
