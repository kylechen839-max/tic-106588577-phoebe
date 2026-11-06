"""Cleaned multi-sector QLP light curve for J112238.89-592027.5 = TIC 451126962.

QLP SAP flux, cadences with QLP quality == 0 only (QLP's scattered-light bits 4096 and 2^29-2^30 are
not in lightkurve's default mask). QLP's own detrended fluxes (KSPSAP / DET_FLUX) are not used:
their splines partly fit the 8-h primary and leave +0.5% shoulders beside it. Instead each segment
(split at gaps > 0.5 d) is divided by a running median of the out-of-eclipse flux (window 1.0 d,
eclipses masked at |phase| < 0.04 and |phase - 0.5| < 0.04, linear interpolation across them).
Errors: QLP has no SAP error column, so flux_err is each segment's out-of-eclipse point-to-point
scatter. Then 5-sigma (MAD) clipping of residuals from a phase-binned median, per sector.
Ephemeris from docs/J1122_EPHEMERIS_CORRECTION.md: P 4.7917045 d, T0 BTJD 3762.8427.
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

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
NAME, TIC, P, T0 = "J112238.89-592027.5", 451126962, 4.7917045, 3762.8427
HW = 0.03  # eclipse half-width in phase for the out-of-eclipse mask (generous)


def phase(t):
    return ((t - T0) / P + 0.25) % 1 - 0.25


res = lk.search_lightcurve(f"TIC {TIC}", author="QLP")
rows, log = [], []
for r in res:
    lc = r.download(quality_bitmask="none")
    sec = int(lc.meta["SECTOR"])
    col = "sap_flux"
    n_all = len(lc)
    lc = lc[lc.quality.value == 0]
    lc = lc.select_flux(col).remove_nans()
    t, f, e = lc.time.value, lc.flux.value.astype(float), lc.flux_err.value.astype(float)
    ph = phase(t)
    ooe = (np.abs(ph) > 0.04) & (np.abs(ph - 0.5) > 0.04)
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
    orbit = seg + 1
    trend = np.full(len(t), np.nan)
    for g in np.unique(seg):
        m = seg == g
        tm, fm = t[m & ooe], f[m & ooe]
        if len(tm) < 10:
            continue
        rm = np.array([np.median(fm[np.abs(tm - x) < 0.5]) for x in tm])
        trend[m] = np.interp(t[m], tm, rm)
    ok = np.isfinite(trend)
    t, f, e, orbit = t[ok], f[ok] / trend[ok], e[ok] / trend[ok], orbit[ok]
    # QLP has no SAP error column: use each segment's out-of-eclipse point-to-point scatter
    ooe = ooe[ok]
    for g in np.unique(orbit):
        m = orbit == g
        if not np.isfinite(e[m]).any():
            e[m] = 1.4826 * np.median(np.abs(np.diff(f[m & ooe]))) / np.sqrt(2)
    exp = float(r.exptime.value[0]) if hasattr(r.exptime, "value") else float(r.exptime[0])
    rows.append(pd.DataFrame(dict(time_btjd=t, flux=f, flux_err=e, sector=sec, orbit=orbit, exptime_s=exp)))
    log.append(dict(sector=sec, exptime_s=exp, flux_col=col, n_all=n_all, n_quality0=len(t)))
    print(f"S{sec:03d} QLP {exp:5.0f}s {col} n={len(t)}")

d = pd.concat(rows, ignore_index=True).sort_values("time_btjd", ignore_index=True)
d["phase"] = phase(d.time_btjd.values)

keep = np.ones(len(d), bool)
for _ in range(3):
    b = np.floor((d.phase.values + 0.25) * 500).astype(int)
    out = np.zeros(len(d), bool)
    for s in d.sector.unique():
        m = (d.sector == s).values
        med = pd.Series(np.where(keep & m, d.flux, np.nan)).groupby(b).transform("median").values
        r = d.flux.values - med
        sig = 1.4826 * np.nanmedian(np.abs(r[m & keep]))
        out |= m & (np.abs(r) > 5 * sig)
    keep = ~out
d["clipped"] = ~keep
d.to_csv(OUT / f"{NAME}_qlp_all_with_flags.csv", index=False)
c = d[keep].drop(columns="clipped")
c.to_csv(OUT / f"{NAME}_qlp_clean.csv", index=False)


def depth(a, centre, w=0.005):
    ph = a.phase.values
    inn = np.abs(ph - centre) < w
    near = (np.abs(ph - centre) > 0.04) & (np.abs(ph - centre) < 0.08)
    return 1 - np.median(a.flux[inn]) / np.median(a.flux[near]), int(inn.sum())


lg = pd.DataFrame(log).sort_values("sector")
stats = []
for s in lg.sector:
    a = c[c.sector == s]
    ooe = (np.abs(a.phase) > HW) & (np.abs(a.phase - 0.5) > HW)
    d1, n1 = depth(a, 0.0); d2, n2 = depth(a, 0.5)
    stats.append(dict(n_clean=len(a), n_clipped=int(d[(d.sector == s)].clipped.sum()),
                      rms_ppt=1e3 * 1.4826 * np.median(np.abs(np.diff(a.flux[ooe]))) / np.sqrt(2),
                      primary_depth_pct=100 * d1, n_in_primary=n1, secondary_depth_pct=100 * d2, n_in_secondary=n2))
lg = pd.concat([lg.reset_index(drop=True), pd.DataFrame(stats)], axis=1)
lg.to_csv(OUT / f"{NAME}_qlp_sector_log.csv", index=False)
print(lg.round(3).to_string(index=False))
d1, _ = depth(c, 0.0); d2, _ = depth(c, 0.5)
print(f"all sectors: primary {100 * d1:.2f}%, secondary {100 * d2:.2f}%, clean points {len(c)}, clipped {int(d.clipped.sum())}")

bins = np.linspace(-0.25, 0.75, 401); mid = 0.5 * (bins[1:] + bins[:-1])
g = c.groupby(np.digitize(c.phase, bins))
bm = pd.DataFrame(dict(phase=mid[g.flux.median().index.values - 1], flux=g.flux.median().values,
                       err=g.flux.agg(lambda x: 1.4826 * np.median(np.abs(x - np.median(x))) / np.sqrt(len(x))).values,
                       n=g.size().values))
bm.to_csv(OUT / f"{NAME}_qlp_folded_binned400.csv", index=False)

secs = sorted(c.sector.unique())
fig, axs = plt.subplots(len(secs), 1, figsize=(11, 1.7 * len(secs)), squeeze=False)
for ax, s in zip(axs[:, 0], secs):
    a = d[d.sector == s]
    ax.plot(a.time_btjd[~a.clipped], a.flux[~a.clipped], ".", ms=1, c="k")
    ax.plot(a.time_btjd[a.clipped], a.flux[a.clipped], "x", ms=3, c="r")
    ax.set_ylabel(f"S{s}", fontsize=7); ax.tick_params(labelsize=7)
axs[-1, 0].set_xlabel("BTJD")
fig.suptitle(f"{NAME} (TIC {TIC}) QLP, per-orbit normalised; red = clipped", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_qlp_sectors.png", dpi=110); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.3), gridspec_kw=dict(width_ratios=[2, 1, 1]))
cols = plt.cm.viridis(np.linspace(0, 0.9, len(secs)))
for col_, s in zip(cols, secs):
    a = c[c.sector == s]
    for ax in axs:
        ax.plot(a.phase, a.flux, ".", ms=1, color=col_, alpha=0.35, label=f"S{s}" if ax is axs[0] else None)
for ax, (lo, hi, t) in zip(axs, ((-0.25, 0.75, "full orbit"), (-0.04, 0.04, "primary"), (0.46, 0.54, "secondary"))):
    ax.plot(bm.phase, bm.flux, "r-", lw=1); ax.set_xlim(lo, hi); ax.set_title(t, fontsize=8); ax.set_xlabel("phase")
lo_ = np.nanpercentile(c.flux, 0.3); axs[0].set_ylim(lo_ - 0.005, np.nanpercentile(c.flux, 99.9) + 0.005)
for ax in axs[1:]:
    ax.set_ylim(lo_ - 0.005, np.nanpercentile(c.flux, 99.9) + 0.005)
axs[0].legend(markerscale=8, fontsize=7, ncol=2); axs[0].set_ylabel("normalised flux")
fig.suptitle(f"{NAME} QLP 7 sectors folded, P = {P} d, T0 = {T0} BTJD (red = 400-bin median)", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_qlp_folded.png", dpi=110); plt.close(fig)
