#!/usr/bin/env python3
"""Fetch, clean, period-refine and phase-bin TESS photometry for a DD x EB candidate.

Generic version of the J0719/J1122 extraction + robust sector-consensus binning.
Uses pipeline light curves (SPOC 2-min, TESS-SPOC FFI, QLP FFI) through
lightkurve; falls back to a simple TESScut aperture extraction.

Outputs in outputs/candidates/<designation>/:
  <name>_raw.txt       time flux flux_err sector
  <name>_binned.txt    Time Flux Flux_Err Phase N NSector  (PHOEBE input)
  <name>_ephemeris.json
  <name>_lc_diagnostic.png
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
import lightkurve as lk
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[2]
FLUX_COLUMN = {"SPOC": "pdcsap_flux", "TESS-SPOC": "pdcsap_flux", "QLP": "kspsap_flux"}


def robust_sigma(x):
    x = np.asarray(x, float)
    mad = np.nanmedian(np.abs(x - np.nanmedian(x)))
    s = 1.4826 * mad
    return float(s) if np.isfinite(s) and s > 0 else float(np.nanstd(x))


def download(target, coord, max_sectors=12):
    """Return list of (sector, time, flux, err) from the best pipeline products."""
    res = lk.search_lightcurve(target if target else coord, radius=5 * u.arcsec if not target else None)
    if len(res) == 0 and target:
        res = lk.search_lightcurve(coord, radius=5 * u.arcsec)
    out = {}
    if len(res):
        tab = res.table
        for author in ["SPOC", "TESS-SPOC", "QLP"]:
            sel = np.where(np.array([str(a) for a in tab["author"]]) == author)[0]
            for i in sel:
                sector = int(str(res[i].mission[0]).split()[-1])
                if sector in out:
                    continue
                if len(out) >= max_sectors:
                    break
                lc = None
                for col in (FLUX_COLUMN[author], "sap_flux", "flux"):
                    try:
                        lc = res[i].download(flux_column=col, quality_bitmask="default")
                        break
                    except Exception as exc:
                        err = exc
                if lc is None:
                    print(f"  sector {sector} {author} download failed: {str(err).splitlines()[0]}")
                    continue
                if lc is None:
                    continue
                t = np.asarray(lc.time.value, float)
                f = np.asarray(lc.flux.value, float)
                e = np.asarray(lc.flux_err.value, float) if lc.flux_err is not None else np.full_like(f, np.nan)
                good = np.isfinite(t) & np.isfinite(f)
                if good.sum() < 100:
                    continue
                out[sector] = (t[good], f[good], e[good], author)
                print(f"  sector {sector}: {author} {good.sum()} points")
    if not out:
        out = tesscut_fallback(coord, max_sectors)
    return out


def tesscut_fallback(coord, max_sectors):
    print("  no pipeline light curves; using TESScut 3x3 aperture")
    out = {}
    res = lk.search_tesscut(coord)
    for i in range(min(len(res), max_sectors)):
        try:
            tpf = res[i].download(cutout_size=15)
        except Exception as exc:
            print(f"  tesscut {i} failed: {exc}")
            continue
        sector = int(tpf.sector)
        ap = np.zeros(tpf.shape[1:], bool)
        cy, cx = ap.shape[0] // 2, ap.shape[1] // 2
        ap[cy - 1:cy + 2, cx - 1:cx + 2] = True
        lc = tpf.to_lightcurve(aperture_mask=ap)
        # background: median of faint pixels per cadence
        flat = tpf.flux.value.reshape(len(tpf.time), -1)
        med_img = np.nanmedian(flat, axis=0)
        bg_pix = med_img < np.nanpercentile(med_img, 30)
        bg = np.nanmedian(flat[:, bg_pix], axis=1) * ap.sum()
        f = lc.flux.value - bg
        t = lc.time.value
        good = np.isfinite(t) & np.isfinite(f) & (tpf.quality == 0)
        if good.sum() < 100:
            continue
        out[sector] = (t[good], f[good], np.full(good.sum(), np.nan), "TESScut")
    return out


def normalise(segments, gap=0.75):
    """Median-normalise each continuous segment (TESS orbit) and clip upward outliers only."""
    T, F, E, S = [], [], [], []
    for sector, (t, f, e, _) in sorted(segments.items()):
        o = np.argsort(t)
        t, f, e = t[o], f[o], e[o]
        breaks = np.where(np.diff(t) > gap)[0] + 1
        for idx in np.split(np.arange(len(t)), breaks):
            if len(idx) < 30:
                continue
            ff = f[idx]
            med = np.nanmedian(ff)
            ff = ff / med
            ee = e[idx] / med
            s = robust_sigma(ff)
            keep = ff < 1 + 5 * s  # keep eclipses; drop flares / cosmic rays
            # isolated low/high points (bad cadences): compare with a running median of neighbours.
            # Eclipses are smooth over many consecutive cadences, so they survive this filter.
            from scipy.ndimage import median_filter
            local = median_filter(ff, size=9, mode="nearest")
            dev = ff - local
            keep &= np.abs(dev) < 6 * robust_sigma(dev)
            keep &= ff > 0.05
            T.append(t[idx][keep]); F.append(ff[keep]); E.append(ee[keep]); S.append(np.full(keep.sum(), sector))
    t, f, e, s = map(np.concatenate, (T, F, E, S))
    bad = ~np.isfinite(e) | (e <= 0)
    e[bad] = robust_sigma(f)
    return t, f, e, s


def binned_dispersion(phase, flux, nbins=200):
    idx = np.floor(phase * nbins).astype(int) % nbins
    total, n = 0.0, 0
    for k in range(nbins):
        m = idx == k
        if m.sum() > 1:
            total += np.sum((flux[m] - np.median(flux[m])) ** 2)
            n += m.sum()
    return total / max(n, 1)


def refine_period(t, f, p0, frac=0.003, n=301):
    trial = np.linspace(p0 * (1 - frac), p0 * (1 + frac), n)
    disp = [binned_dispersion(((t - t[0]) / p) % 1, f) for p in trial]
    best = trial[int(np.argmin(disp))]
    # second, finer pass
    trial2 = np.linspace(best * (1 - frac / 20), best * (1 + frac / 20), 201)
    disp2 = [binned_dispersion(((t - t[0]) / p) % 1, f) for p in trial2]
    return float(trial2[int(np.argmin(disp2))])


def find_t0(t, f, p, nbins=200):
    ph = ((t - t[0]) / p) % 1
    idx = np.floor(ph * nbins).astype(int)
    prof = np.array([np.median(f[idx == k]) if np.any(idx == k) else np.nan for k in range(nbins)])
    k = int(np.nanargmin(prof))
    ks = [(k + d) % nbins for d in (-1, 0, 1)]
    y = prof[ks]
    denom = y[0] - 2 * y[1] + y[2]
    off = 0.5 * (y[0] - y[2]) / denom if np.isfinite(denom) and denom != 0 else 0.0
    phase_min = (k + 0.5 + np.clip(off, -1, 1)) / nbins
    t0 = t[0] + phase_min * p
    # move t0 near the middle of the data
    ncyc = np.round((np.median(t) - t0) / p)
    return float(t0 + ncyc * p), prof


def odd_even_test(t, f, p, t0, core_frac=None):
    """Per-cycle primary-eclipse depths, split into odd and even cycles.

    Returns (n_sigma difference, depth_odd, depth_even, quarter_dip) where the error
    comes from the eclipse-to-eclipse scatter, so sector systematics are included.
    quarter_dip is the depth at phase 0.5 of P relative to the primary depth;
    if the folded-at-P light curve already has a clear secondary at 0.5, the
    catalogue period is right even if odd/even scatter is large.
    """
    ph = ((t - t0) / p + 0.5) % 1 - 0.5
    cyc = np.round((t - t0) / p).astype(int)
    # core half-width from the folded profile: points below half depth
    oot_level = np.median(f[np.abs(ph) > 0.15])
    prim = np.median(f[np.abs(ph) < 0.004]) if np.any(np.abs(ph) < 0.004) else np.min(f)
    depth = oot_level - prim
    if core_frac is None:
        inside = np.abs(ph[(f < oot_level - 0.5 * depth) & (np.abs(ph) < 0.15)])
        core_frac = max(np.percentile(inside, 50) if len(inside) > 5 else 0.005, 0.002)
    depths = {}
    for c in np.unique(cyc):
        m = cyc == c
        core = m & (np.abs(ph) < core_frac)
        ring = m & (np.abs(ph) > 3 * core_frac) & (np.abs(ph) < 6 * core_frac + 0.02)
        if core.sum() >= 3 and ring.sum() >= 5:
            depths[c] = np.median(f[ring]) - np.median(f[core])
    if len(depths) < 4:
        return np.nan, np.nan, np.nan, np.nan
    c = np.array(list(depths)); d = np.array(list(depths.values()))
    odd, even = d[c % 2 == 1], d[c % 2 == 0]
    if len(odd) < 2 or len(even) < 2:
        return np.nan, np.nan, np.nan, np.nan
    err = np.sqrt(np.var(odd, ddof=1) / len(odd) + np.var(even, ddof=1) / len(even))
    ph2 = ((t - t0) / p) % 1
    sec = np.median(f[np.abs(ph2 - 0.5) < core_frac]) if np.any(np.abs(ph2 - 0.5) < core_frac) else oot_level
    quarter = (oot_level - sec) / depth if depth > 0 else np.nan
    return float(abs(np.mean(odd) - np.mean(even)) / err), float(np.mean(odd)), float(np.mean(even)), float(quarter)


def phase_values(t, t0, p):
    ph = ((t - t0) / p) % 1.0
    return np.where(ph >= 0.8, ph - 1.0, ph)


def locate_eclipses(t, f, t0, p, nfine=1000):
    """Find primary/secondary phase, depth and half-width from a fine folded profile.

    Out-of-eclipse variation (ellipsoidal, reflection, spots) is removed first
    with a 2nd-order Fourier series fitted outside the eclipses, iterating the
    eclipse mask, so widths are not inflated by the continuum curvature.
    """
    ph = ((t - t0) / p) % 1.0
    idx = np.floor(ph * nfine).astype(int)
    prof = np.array([np.median(f[idx == k]) if np.any(idx == k) else np.nan for k in range(nfine)])
    x = (np.arange(nfine) + 0.5) / nfine
    good = np.isfinite(prof)
    X = np.vstack([np.ones(nfine)] + [fn(2 * np.pi * h * x) for h in (1, 2) for fn in (np.cos, np.sin)]).T
    dphase = lambda a, b: np.abs(((a - b + 0.5) % 1) - 0.5)
    mask = good & (dphase(x, 0.0) > 0.1) & (dphase(x, 0.5) > 0.1)
    res = prof.copy()
    hw1_prev = 0.05
    for _ in range(4):
        coef, *_ = np.linalg.lstsq(X[mask], prof[mask], rcond=None)
        cont = X @ coef
        res = prof - cont
        noise = robust_sigma(res[mask])
        k1 = int(np.nanargmin(np.where(dphase(x, 0.0) < 0.05, res, np.nan)))
        d1 = -res[k1]
        ecl = good & (res < -max(3 * noise, 0.05 * d1))
        far1 = dphase(x, x[k1]) > max(0.1, 1.5 * hw1_prev)
        cand = np.where(far1 & good, res, np.nan)
        k2 = int(np.nanargmin(cand))
        d2 = -res[k2]
        # grow the mask around both eclipses until the residual returns to the continuum
        new_mask = good.copy()
        widths = []
        for k, d in ((k1, d1), (k2, d2)):
            thr = -max(3 * noise, 0.05 * d)
            w = 0
            for step in range(1, nfine // 4):
                a, b = res[(k - step) % nfine], res[(k + step) % nfine]
                if (np.isfinite(a) and a < thr) or (np.isfinite(b) and b < thr):
                    w = step
                elif step > w + 4:
                    break
            hw = (w + 1) / nfine
            widths.append(hw)
            new_mask &= dphase(x, x[k]) > 1.3 * hw + 0.005
        mask = new_mask
        hw1_prev = widths[0]
    # centroid of each dip (more precise than the minimum bin)
    def centroid(k, hw):
        sel = good & (dphase(x, x[k]) < hw) & (res < 0)
        if sel.sum() < 3:
            return x[k]
        off = ((x[sel] - x[k] + 0.5) % 1) - 0.5
        return (x[k] + np.sum(off * -res[sel]) / np.sum(-res[sel])) % 1
    c1, c2 = centroid(k1, widths[0]), centroid(k2, widths[1])
    phase2 = float((c2 - c1) % 1.0)
    detected = bool(d2 > 4 * noise)
    return {"phase1": float(((c1 + 0.5) % 1) - 0.5), "phase2": phase2, "depth1": float(d1), "depth2": float(d2),
            "halfwidth1": float(widths[0]), "halfwidth2": float(widths[1] if detected else widths[0]),
            "secondary_detected": detected, "continuum_amplitude": float(np.ptp(cont)), "profile_noise": float(noise)}


def adaptive_edges(ecl, bins, min_width):
    """Uniform bins out of eclipse, ~30 bins across each eclipse."""
    edges = list(np.linspace(-0.2, 0.8, bins + 1))
    for c, w in ((0.0, ecl["halfwidth1"]), (ecl["phase2"], ecl["halfwidth2"])):
        c = c if c < 0.8 else c - 1.0
        lo, hi = c - 1.3 * w, c + 1.3 * w
        fine_w = max(2.6 * w / 30, min_width)
        edges = [e for e in edges if not (lo < e < hi)] + list(np.arange(lo, hi + 1e-12, fine_w))
    edges = np.array(sorted(e for e in edges if -0.2 <= e <= 0.8))
    keep = np.concatenate([[True], np.diff(edges) > 1e-9])
    return edges[keep]


def sector_consensus_bin(t, f, s, t0, p, bins=200, min_points=3, edges=None):
    ph = phase_values(t, t0, p)
    sectors = np.unique(s)
    need = min(2, len(sectors))
    if edges is None:
        edges = np.linspace(-0.2, 0.8, bins + 1)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (ph >= lo) & (ph < hi)
        meds, ns = [], 0
        for sec in sectors:
            sm = m & (s == sec)
            if sm.sum() >= min_points:
                meds.append(np.median(f[sm])); ns += sm.sum()
        if len(meds) < need:
            # sparse phase coverage: fall back to all points if enough
            if m.sum() >= min_points:
                meds = [np.median(f[m])]; ns = int(m.sum())
            else:
                continue
        mid = float(np.median(ph[m]))
        flux = float(np.median(meds))
        if len(meds) >= 2:
            sig = robust_sigma(meds) / np.sqrt(len(meds))
        else:
            sig = robust_sigma(f[m]) / np.sqrt(max(m.sum(), 1))
        rows.append((t0 + mid * p, flux, sig, mid, ns, len(meds)))
    rows = np.array(rows)
    floor = max(np.median(rows[:, 2]), 2e-4)
    rows[:, 2] = np.maximum(rows[:, 2], floor)
    return rows


def eclipse_depths(rows, phase2=0.5):
    ph, fl = rows[:, 3], rows[:, 1]
    p2 = phase2 if phase2 < 0.8 else phase2 - 1
    far = (np.abs(ph) > 0.1) & (np.abs(ph - p2) > 0.1)
    oot = np.median(fl[far]) if np.any(far) else np.nan
    prim = oot - np.min(fl[np.abs(ph) < 0.08]) if np.any(np.abs(ph) < 0.08) else np.nan
    sec = oot - np.min(fl[np.abs(ph - p2) < 0.05]) if np.any(np.abs(ph - p2) < 0.05) else np.nan
    return float(oot), float(prim), float(sec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designation", required=True)
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--tic", default=None)
    ap.add_argument("--period", type=float, required=True, help="catalogue period (d)")
    ap.add_argument("--bins", type=int, default=200)
    ap.add_argument("--max-sectors", type=int, default=12)
    ap.add_argument("--no-double", action="store_true", help="do not test 2P")
    ap.add_argument("--fixed-period", action="store_true",
                    help="use --period as given (e.g. from an eclipse-timing O-C fit) instead of refining")
    args = ap.parse_args()

    out = ROOT / "outputs" / "candidates" / args.designation
    out.mkdir(parents=True, exist_ok=True)
    coord = SkyCoord(args.ra, args.dec, unit="deg")
    target = f"TIC {args.tic}" if args.tic and args.tic != "nan" else None
    print(f"{args.designation}: downloading TESS data ({target or coord.to_string()})")
    segs = download(target, coord, args.max_sectors)
    if not segs:
        raise SystemExit("no TESS data")
    t, f, e, s = normalise(segs)
    np.savetxt(out / f"{args.designation}_raw.txt", np.c_[t, f, e, s], header="time flux flux_err sector", fmt="%.8f")

    p = args.period if args.fixed_period else refine_period(t, f, args.period)
    t0, _ = find_t0(t, f, p)
    oe = (np.nan,) * 4
    if not args.no_double:
        oe = odd_even_test(t, f, p, t0)
        print(f"  odd/even depth test at P: {oe[0]:.1f} sigma (odd {oe[1]:.4f}, even {oe[2]:.4f}), "
              f"phase-0.5 dip / primary = {oe[3]:.2f}")
        # Only double when odd/even primaries differ significantly AND there is no secondary at phase 0.5
        if np.isfinite(oe[0]) and oe[0] > 4 and (not np.isfinite(oe[3]) or oe[3] < 0.15):
            print("  odd and even eclipses differ: adopting 2P")
            p = refine_period(t, f, 2 * p, frac=0.0005)
            t0, _ = find_t0(t, f, p)
    p2_chi2 = oe[0]
    ecl = locate_eclipses(t, f, t0, p, nfine=int(np.clip(len(t) / 8, 200, 1000)))
    cadence = np.median(np.diff(np.sort(t)))
    edges = adaptive_edges(ecl, args.bins, min_width=1.5 * cadence / p)
    rows = sector_consensus_bin(t, f, s, t0, p, edges=edges)
    print(f"  eclipses: {json.dumps(ecl)}; {len(rows)} adaptive bins")
    np.savetxt(out / f"{args.designation}_binned.txt", rows, header="Time Flux Flux_Err Phase N NSector", fmt="%.8e")
    oot, d1, d2 = eclipse_depths(rows, ecl["phase2"])
    eph = {"designation": args.designation, "catalog_period": args.period, "period": p, "t0": t0,
           "odd_even_sigma": p2_chi2, "odd_depth": oe[1], "even_depth": oe[2],
           "phase05_dip_over_primary": oe[3], "sectors": sorted(int(x) for x in segs),
           "authors": sorted(set(v[3] for v in segs.values())), "n_points": int(len(t)),
           "n_bins": int(len(rows)), "oot_flux": oot, "primary_depth": d1, "secondary_depth": d2,
           "raw_rms": robust_sigma(f), "eclipses": ecl,
           "ecosw_estimate": float(np.pi / 2 * (ecl["phase2"] - 0.5))}
    (out / f"{args.designation}_ephemeris.json").write_text(json.dumps(eph, indent=2) + "\n")

    ph = phase_values(t, t0, p)
    fig, ax = plt.subplots(2, 1, figsize=(9, 7))
    for sec in np.unique(s):
        m = s == sec
        ax[0].plot(t[m], f[m], ".", ms=1, label=f"S{int(sec)}")
        ax[1].plot(ph[m], f[m], ".", ms=1, alpha=0.3)
    ax[1].plot(rows[:, 3], rows[:, 1], "k.", ms=4, label="sector-consensus bins")
    ax[0].set_xlabel("BTJD"); ax[0].set_ylabel("norm flux"); ax[0].legend(fontsize=6, ncol=6, markerscale=5)
    ax[1].set_xlabel("phase"); ax[1].set_ylabel("norm flux"); ax[1].legend(fontsize=8)
    ax[1].set_title(f"P={p:.6f} d  t0={t0:.5f}  depths {d1:.4f}/{d2:.4f}", fontsize=9)
    fig.suptitle(args.designation)
    fig.tight_layout(); fig.savefig(out / f"{args.designation}_lc_diagnostic.png", dpi=140); plt.close(fig)
    print(json.dumps(eph, indent=2))


if __name__ == "__main__":
    main()
