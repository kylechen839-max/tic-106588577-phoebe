#!/usr/bin/env python3
"""Generic PHOEBE pipeline for a DD x EB candidate.

Same workflow that was validated on J1122 and on V501 Mon / W UMa
(docs/ECB_TESTDATA_LITERATURE_COMPARISON.md):

  sector-consensus binned TESS LC -> coarse grid -> bounded Powell -> short emcee
  -> high-resolution recompute -> post-fit checks

Post-fit checks (written to <name>_phoebe_result.json):
  photometry  : reduced chi2, RMS, mean residual inside each eclipse (bias)
  blackbody   : the fitted binary (R1, R2, T1, T2 in absolute units, from the
                assumed masses + Kepler's law) as a two-blackbody photosphere is
                scaled to J,H,Ks,W1; the implied photometric distance is compared
                with Gaia, and the W3/W4 excess is recomputed against the binary
                photosphere instead of a single star.

Run with the PHOEBE environment, e.g.
  ~/phoebe-env/bin/python src/pipeline/run_phoebe_candidate.py --designation J045912.77+165543.8
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.optimize import minimize

import phoebe
from phoebe import u as pu

ROOT = Path(__file__).resolve().parents[2]
SCORED = ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv"
ALLWISE = ROOT / "outputs" / "crossmatch" / "dd1_support_cache" / "allwise.csv"
EXTRA_SCORED = ROOT / "outputs" / "crossmatch" / "extra_scored.csv"
EXTRA_ALLWISE = ROOT / "outputs" / "crossmatch" / "extra_support_allwise.csv"

ZP = {"J": 1594.0, "H": 1024.0, "K": 666.7, "W1": 309.540, "W2": 171.787, "W3": 31.674, "W4": 8.363}
LAM = {"J": 1.235, "H": 1.662, "K": 2.159, "W1": 3.35, "W2": 4.60, "W3": 11.56, "W4": 22.09}
RSUN_M = 6.957e8
PC_M = 3.0857e16

PARAMS = ["incl", "dt0", "rsum", "k", "teffratio"]


# ----------------------------------------------------------------------------- stellar priors
def ms_mass_radius(teff):
    """Rough main-sequence mass and radius from Teff (Eker et al. 2018 style relations)."""
    logt = math.log10(teff)
    # piecewise log M - log Teff relation for MS stars
    if teff < 3900:
        m = 0.45 * (teff / 3800) ** 2
    elif teff < 6000:
        m = 10 ** (2.0 * (logt - math.log10(5800)))
    elif teff < 10000:
        m = 10 ** (1.7 * (logt - math.log10(5800))) * 1.0
    else:
        m = 10 ** (2.3 * (logt - math.log10(10000))) * 2.5
    m = float(np.clip(m, 0.2, 40))
    r = m ** 0.8 if m < 1.4 else m ** 0.65 * 1.07
    return m, float(r)


def kepler_sma(period_d, mtot):
    return (mtot * (period_d / 365.25) ** 2) ** (1 / 3) * 215.032  # Rsun


# ----------------------------------------------------------------------------- data
def load_inputs(name):
    cdir = ROOT / "outputs" / "candidates" / name
    eph = json.loads((cdir / f"{name}_ephemeris.json").read_text())
    data = np.genfromtxt(cdir / f"{name}_binned.txt", comments="#")
    times, flux, sig = data[:, 0], data[:, 1], data[:, 2]
    sig = np.maximum(sig, float(os.environ.get("SIGMA_FLOOR", "0.0003")))
    return cdir, eph, times, flux, sig


def phase_values(times, t0, period):
    ph = ((times - t0) / period) % 1.0
    return np.where(ph >= 0.8, ph - 1.0, ph)


# ----------------------------------------------------------------------------- bundle
class Model:
    def __init__(self, times, flux, sig, period, t0, teff1, mass1, q, distortion, ntri, l3=0.0):
        self.period, self.t0, self.teff1 = period, t0, teff1
        self.mtot = mass1 * (1 + q)
        self.sma = kepler_sma(period, self.mtot)
        b = phoebe.default_binary()
        b.add_dataset("lc", times=times, fluxes=flux, sigmas=sig, dataset="lc01", passband="TESS:T")
        b.set_value_all("atm", "blackbody")
        b.set_value_all("ld_mode_bol", "manual")
        b.set_value_all("ld_func_bol", "linear")
        b.set_value_all("ld_coeffs_bol", [0.5])
        b.set_value_all("ld_mode", "manual")
        b.set_value_all("ld_func", "linear")
        b.set_value_all("ld_coeffs", [0.5 if teff1 < 9000 else 0.35])
        b.set_value("period@binary", period)
        b.set_value("t0_supconj@binary", t0)
        b.set_value("ecc@binary", 0.0)
        b.set_value("q@binary", q)
        b.set_value("sma@binary", self.sma)
        b.set_value("teff@primary", teff1)
        gb = 1.0 if teff1 > 7500 else 0.32
        alb = 1.0 if teff1 > 7500 else 0.6
        for comp in ("primary", "secondary"):
            b.set_value(f"gravb_bol@{comp}", gb)
            b.set_value(f"irrad_frac_refl_bol@{comp}", alb)
            b.set_value(f"distortion_method@{comp}@phoebe01", distortion)
        b.set_value("pblum_mode@lc01", "dataset-scaled")
        b.set_value_all("ntriangles", ntri)
        b.set_value("l3_mode@lc01", "fraction")
        b.set_value("l3_frac@lc01", l3)
        self.b = b

    def set(self, v):
        d = dict(zip(PARAMS, v))
        incl, dt0, rsum, k, tr = (d[n] for n in PARAMS[:5])
        if "l3" in d:
            self.b.set_value("l3_frac@lc01", d["l3"])
        if "q" in d:
            self.b.set_value("q@binary", d["q"])
        if "ecosw" in d:
            ecc = math.hypot(d["ecosw"], d["esinw"])
            if ecc >= 0.9:
                raise ValueError("ecc too large")
            self.b.set_value("ecc@binary", ecc)
            self.b.set_value("per0@binary", math.degrees(math.atan2(d["esinw"], d["ecosw"])) % 360)
        r1 = rsum / (1 + k) * self.sma
        r2 = rsum * k / (1 + k) * self.sma
        b = self.b
        b.set_value("incl@binary", incl)
        b.set_value("t0_supconj@binary", self.t0 + dt0)
        b.set_value("requiv@primary", r1)
        b.set_value("requiv@secondary", r2)
        b.set_value("teff@secondary", self.teff1 * tr)

    def evaluate(self, v, label="candmodel", keep=False):
        try:
            self.set(v)
            self.b.run_compute(compute="phoebe01", model=label, overwrite=True, progressbar=False)
            pred = np.asarray(self.b.get_value(f"fluxes@lc01@{label}@model"), float)
            obs = np.asarray(self.b.get_value("fluxes@lc01@dataset"), float)
            sig = np.asarray(self.b.get_value("sigmas@lc01@dataset"), float)
            chi2 = float(np.sum(((obs - pred) / sig) ** 2))
            if not np.isfinite(chi2):
                return np.inf, None
            if not keep:
                self.b.remove_model(label)
            return chi2, pred
        except Exception:
            try:
                self.b.remove_model(label)
            except Exception:
                pass
            return np.inf, None


# ----------------------------------------------------------------------------- fitting
def initial_guesses(eph, flux, times):
    """Seed values from the eclipse shape: width -> rsum, depth ratio -> teffratio."""
    ph = phase_values(times, eph["t0"], eph["period"])
    oot = eph.get("oot_flux", 1.0)
    d1 = max(eph.get("primary_depth", 0.05), 1e-3)
    d2 = max(eph.get("secondary_depth", 0.0), 0.0)
    ecl = eph.get("eclipses")
    if ecl and ecl.get("halfwidth1"):
        # half-width to first contact (i ~ 90 deg, circular): r1 + r2 ~ sin(2 pi hw)
        width = 2 * ecl["halfwidth1"]
        # the 5%-depth threshold reaches slightly beyond first contact, hence 0.85
        rsum = float(np.clip(0.85 * math.sin(2 * math.pi * min(ecl["halfwidth1"], 0.2)), 0.04, 0.7))
    else:
        below = ph[(np.abs(ph) < 0.2) & (flux < oot - 0.15 * d1)]
        width = (below.max() - below.min()) if len(below) > 2 else 0.06
        rsum = float(np.clip(np.pi * width / 1.9, 0.06, 0.7))
    tr = float(np.clip((d2 / d1) ** 0.25 if d2 > 0 else 0.6, 0.35, 1.05))
    return rsum, tr, width


def bounds_for(rsum0, distortion):
    return {
        "incl": (60.0, 90.0),
        "dt0": (-0.01, 0.01),
        "rsum": (max(0.03, 0.4 * rsum0), min(0.78 if distortion == "roche" else 0.9, 2.2 * rsum0)),
        "k": (0.2, 1.4),
        "teffratio": (0.3, 1.25),
        "l3": (0.0, 0.8),
        "ecosw": (-0.6, 0.6),
        "esinw": (-0.7, 0.7),
        "q": (0.15, 1.5),
    }


def in_bounds(v, B):
    return all(B[n][0] <= x <= B[n][1] for n, x in zip(PARAMS, v))


def grid_search(m, rsum0, tr0, B, log, ntop=3):
    """Coarse grid; returns the best (chi2, v) and the top-ntop distinct grid points."""
    results = []
    incls = [74, 80, 84, 87, 89.5]
    rsums = [rsum0 * f for f in (0.7, 0.85, 1.0, 1.2, 1.4)]
    ks = [0.5, 0.75, 1.0, 1.3]
    trs = sorted(set([float(np.clip(tr0 * f, 0.3, 1.2)) for f in (0.8, 1.0, 1.2)]))
    best = np.inf
    for incl in incls:
        for rs in rsums:
            for k in ks:
                for tr in trs:
                    extra = {"l3": m.l3_start, "ecosw": m.ecosw_start, "esinw": m.esinw_start, "q": m.q_start}
                    v = np.array([incl, m.dt0_start, rs, k, tr] + [extra[n] for n in PARAMS[5:]])
                    if not in_bounds(v, B):
                        continue
                    chi2, _ = m.evaluate(v)
                    if np.isfinite(chi2):
                        results.append((chi2, v))
                        if chi2 < best:
                            best = chi2
                            log(f"  grid {len(results)}: chi2={chi2:.1f} v={np.round(v, 4).tolist()}")
    results.sort(key=lambda r: r[0])
    top = []
    for c, v in results:
        # distinct = differs from every kept point in incl by >2 deg or in k or rsum
        if all(abs(v[0] - t[1][0]) > 2 or abs(v[3] - t[1][3]) > 0.05 or abs(v[2] - t[1][2]) > 0.01 for t in top):
            top.append((c, v))
        if len(top) >= ntop:
            break
    return (top[0] if top else (np.inf, None)), top


def powell(m, v0, B, maxiter, log, maxfev=None):
    """Local optimisation from the grid point (name kept for the pipeline stage).

    Nelder-Mead in normalised [0, 1] coordinates with a small initial simplex
    (2% of each parameter range), bounds enforced by a penalty. scipy's bounded
    Powell was replaced because its line searches span the whole feasible range
    and miss the narrow chi2 basin on a rugged PHOEBE surface. Returns the best
    vector actually evaluated.
    """
    lo = np.array([B[p][0] for p in PARAMS])
    hi = np.array([B[p][1] for p in PARAMS])
    span = hi - lo
    best = {"chi2": np.inf, "v": np.array(v0, float)}
    n = {"evals": 0}

    def f(x):
        if np.any(x < 0) or np.any(x > 1):
            return 1e15
        v = lo + x * span
        c, _ = m.evaluate(v)
        n["evals"] += 1
        if np.isfinite(c) and c < best["chi2"]:
            best.update(chi2=c, v=v.copy())
        if n["evals"] % 25 == 0:
            log(f"  local-opt eval {n['evals']}: best chi2 {best['chi2']:.2f}")
        return c if np.isfinite(c) else 1e15

    x0 = np.clip((np.asarray(v0, float) - lo) / span, 1e-4, 1 - 1e-4)
    ndim = len(x0)
    for restart in range(max(1, maxiter // 3)):
        simplex = [x0]
        for i in range(ndim):
            xi = x0.copy()
            step = 0.02 if x0[i] < 0.98 else -0.02
            xi[i] = np.clip(xi[i] + step, 0, 1)
            simplex.append(xi)
        minimize(f, x0, method="Nelder-Mead",
                 options={"initial_simplex": np.array(simplex), "maxfev": maxfev or 60 * ndim,
                          "xatol": 1e-4, "fatol": 1e-3 * max(best["chi2"], 1)})
        x_new = (best["v"] - lo) / span
        if np.allclose(x_new, x0, atol=1e-5):
            break
        x0 = x_new  # restart NM from the improved point with a fresh simplex
    return best["chi2"], best["v"], n["evals"]


def run_emcee(m, v0, B, nwalkers, niter, scale_chi2, log):
    import emcee
    rng = np.random.default_rng(1)
    steps = {"incl": 0.3, "dt0": 0.0004, "rsum": 0.01 * v0[2], "k": 0.03, "teffratio": 0.02, "l3": 0.02,
             "ecosw": 0.002, "esinw": 0.01, "q": 0.03}
    step = np.array([steps[n] for n in PARAMS])
    p0 = []
    while len(p0) < nwalkers:
        v = v0 + step * rng.normal(size=len(v0))
        if in_bounds(v, B):
            p0.append(v)

    def lnp(v):
        if not in_bounds(v, B):
            return -np.inf
        c, _ = m.evaluate(v)
        return -0.5 * c / scale_chi2 if np.isfinite(c) else -np.inf

    s = emcee.EnsembleSampler(nwalkers, len(v0), lnp)
    t = time.time()
    s.run_mcmc(np.array(p0), niter, progress=False)
    log(f"  emcee {nwalkers}x{niter} in {time.time() - t:.0f}s, acc={np.mean(s.acceptance_fraction):.2f}")
    chain, lp = s.get_chain(), s.get_log_prob()
    i, j = np.unravel_index(np.nanargmax(lp), lp.shape)
    flat = s.get_chain(discard=niter // 2, flat=True)
    return chain[i, j], flat, chain, lp


# ----------------------------------------------------------------------------- post-fit checks
def planck(T, lam_um):
    h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
    nu = c / (np.asarray(lam_um) * 1e-6)
    return 2 * h * nu**3 / c**2 / np.expm1(np.clip(h * nu / (k * T), 1e-8, 700))


def read_rows(path, name):
    if not Path(path).exists():  # e.g. cloud container without the AllWISE support cache
        return []
    with open(path, newline="") as fh:
        return [r for r in csv.DictReader(fh) if r["designation"] == name]


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def blackbody_check(name, R1, R2, T1, T2, gaia_dist_pc):
    rows = read_rows(ALLWISE, name) or (read_rows(EXTRA_ALLWISE, name) if EXTRA_ALLWISE.exists() else [])
    rows = sorted(rows, key=lambda r: float(r["angDist"]))
    if not rows:
        return {"status": "no AllWISE row"}
    a = rows[0]
    bands = ["J", "H", "K", "W1", "W2", "W3", "W4"]
    cols = ["Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag"]
    f, s, lam = [], [], []
    for b_, c in zip(bands, cols):
        mag, err = fnum(a[c]), fnum(a["e_" + c])
        err = max(err if err is not None else 0.1, 0.03)
        f.append(ZP[b_] * 10 ** (-0.4 * mag)); s.append(f[-1] * 0.4 * math.log(10) * err); lam.append(LAM[b_])
    f, s, lam = map(np.array, (f, s, lam))
    # binary photosphere in Jy at distance D: pi (R/D)^2 B_nu * 1e26
    shape = R1**2 * planck(T1, lam) + R2**2 * planck(T2, lam)
    ip = np.arange(4)
    w = s[ip] ** -2
    A = np.sum(w * shape[ip] * f[ip]) / np.sum(w * shape[ip] ** 2)
    model = A * shape
    chi2_phot = float(np.sum(((f[ip] - model[ip]) / s[ip]) ** 2))
    # A = pi * Rsun^2 / D^2 * 1e26  -> D
    d_phot = RSUN_M * math.sqrt(math.pi * 1e26 / A) / PC_M if A > 0 else np.nan
    excess = {b_: {"ratio": float(f[i] / model[i]), "sigma": float((f[i] - model[i]) / s[i])}
              for i, b_ in enumerate(bands) if b_ in ("W2", "W3", "W4")}
    dist_ratio = d_phot / gaia_dist_pc if gaia_dist_pc and np.isfinite(gaia_dist_pc) else np.nan
    return {
        "binary_photosphere_chi2_JHKW1": chi2_phot,
        "photometric_distance_pc": d_phot,
        "gaia_distance_pc": gaia_dist_pc,
        "distance_ratio_phot_over_gaia": dist_ratio,
        "excess_vs_binary_photosphere": excess,
        "pass_photosphere": bool(chi2_phot < 60),
        # extinction makes stars look fainter -> photometric distance too large; allow 0.6-1.8
        "pass_distance": bool(np.isfinite(dist_ratio) and 0.6 < dist_ratio < 1.8) if np.isfinite(dist_ratio) else None,
        "pass_excess_survives": bool(max(excess["W3"]["sigma"], excess["W4"]["sigma"]) >= 5),
        "lam": lam.tolist(), "flux": f.tolist(), "sigma": s.tolist(), "model": model.tolist(),
    }


def photometry_check(times, obs, pred, sig, t0, period, nparams, phase2=0.5):
    ph = phase_values(times, t0, period)
    r = obs - pred
    chi2 = float(np.sum((r / sig) ** 2))
    dof = max(len(obs) - nparams - 1, 1)
    out = {"chi2": chi2, "dof": dof, "chi2_red": chi2 / dof, "rms": float(np.sqrt(np.mean(r**2)))}
    for lab, centre in (("primary", 0.0), ("secondary", phase2 if phase2 < 0.8 else phase2 - 1)):
        m = np.abs(ph - centre) < 0.02
        if m.sum():
            mean = float(np.mean(r[m]))
            err = float(np.sqrt(np.mean(sig[m] ** 2)) / np.sqrt(m.sum()))
            out[f"{lab}_mean_resid"] = mean
            out[f"{lab}_bias_sigma"] = mean / err if err > 0 else np.nan
    out["pass_fit"] = bool(out["chi2_red"] < 5)
    out["pass_no_eclipse_bias"] = bool(all(abs(out.get(f"{l}_bias_sigma", 0)) < 3 for l in ("primary", "secondary")))
    return out


def plot_fit(path, name, times, obs, sig, pred, t0, period, title):
    ph = phase_values(times, t0, period)
    o = np.argsort(ph)
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.errorbar(ph, obs, yerr=sig, fmt=".", color="0.3", ms=3, lw=0.5, label="TESS sector-consensus bins")
    ax.plot(ph[o], pred[o], color="#dc2626", lw=1.8, label="PHOEBE")
    ax.set_ylabel("norm flux"); ax.legend(fontsize=8); ax.set_title(title, fontsize=9)
    rx.axhline(0, color="0.5"); rx.plot(ph, obs - pred, ".", color="0.3", ms=3)
    rx.set_xlabel("phase"); rx.set_ylabel("O-C")
    fig.suptitle(name); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_sed(path, name, bb):
    lam, f, s, mdl = (np.array(bb[k]) for k in ("lam", "flux", "sigma", "model"))
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.errorbar(lam, f, yerr=s, fmt="o", color="k", label="2MASS + AllWISE")
    ax.plot(lam, mdl, "s--", color="#4363d8", label="PHOEBE binary blackbody photosphere")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("micron"); ax.set_ylabel("Jy")
    ax.set_title(f"{name}: d_phot/d_Gaia={bb['distance_ratio_phot_over_gaia']:.2f}", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designation", required=True)
    ap.add_argument("--teff", type=float, default=None)
    ap.add_argument("--q", type=float, default=None)
    ap.add_argument("--distortion", default=None, choices=[None, "roche", "sphere"])
    ap.add_argument("--ntri", type=int, default=400)
    ap.add_argument("--ntri-final", type=int, default=1500)
    ap.add_argument("--powell-maxiter", type=int, default=6)
    ap.add_argument("--walkers", type=int, default=12)
    ap.add_argument("--iters", type=int, default=30)
    ap.add_argument("--l3", type=float, default=0.0, help="third light (fraction of total); start value if --fit-l3")
    ap.add_argument("--fit-l3", action="store_true", help="fit third light as a free parameter")
    ap.add_argument("--fit-q", choices=["auto", "yes", "no"], default="auto",
                    help="fit mass ratio; auto = for Roche-distortion (P < 3 d) systems")
    ap.add_argument("--eccentric", choices=["auto", "yes", "no"], default="auto",
                    help="fit ecosw/esinw; auto = when the secondary eclipse is >0.01 in phase from 0.5")
    ap.add_argument("--period", type=float, default=None, help="override ephemeris period")
    ap.add_argument("--t0", type=float, default=None, help="override ephemeris t0")
    ap.add_argument("--tag", default="", help="suffix for output file names")
    args = ap.parse_args()
    name = args.designation
    cdir, eph, times, flux, sig = load_inputs(name)
    logf = open(cdir / f"{name}_phoebe{args.tag}.log", "a")

    def log(msg):
        print(msg, flush=True); logf.write(msg + "\n"); logf.flush()

    rows = read_rows(SCORED, name) or (read_rows(EXTRA_SCORED, name) if EXTRA_SCORED.exists() else [])
    row = rows[0] if rows else {}
    teff1 = args.teff or fnum(row.get("tic_teff")) or fnum(row.get("gaia_teff")) or fnum(row.get("bb_teff")) or 6000.0
    plx = fnum(row.get("plx"))
    gaia_dist = fnum(row.get("gaia_dist")) or (1000.0 / plx if plx and plx > 0 else np.nan)
    m1, r1_ms = ms_mass_radius(teff1)
    q = args.q if args.q else 0.8
    period, t0 = args.period or eph["period"], args.t0 or eph["t0"]
    if args.fit_l3:
        PARAMS.append("l3")
    ecl = eph.get("eclipses", {})
    phase2 = ecl.get("phase2", 0.5)
    eccentric = args.eccentric == "yes" or (args.eccentric == "auto" and ecl.get("secondary_detected")
                                             and abs(phase2 - 0.5) > 0.01)
    ecosw0 = math.pi / 2 * (phase2 - 0.5) if eccentric else 0.0
    w1, w2 = ecl.get("halfwidth1", 1), ecl.get("halfwidth2", 1)
    esinw0 = float(np.clip((w2 - w1) / (w2 + w1), -0.5, 0.5)) if eccentric else 0.0
    if eccentric:
        PARAMS.extend(["ecosw", "esinw"])
    tag = args.tag
    distortion = args.distortion or ("roche" if period < 3 else "sphere")
    fit_q = args.fit_q == "yes" or (args.fit_q == "auto" and distortion == "roche")
    if fit_q:
        PARAMS.append("q")  # ellipsoidal amplitude depends on q for close (Roche) systems
    log(f"=== {name} {datetime.now(timezone.utc).isoformat()}  P={period:.6f} t0={t0:.5f} Teff1={teff1:.0f} "
        f"M1~{m1:.2f} q={q} distortion={distortion} nbins={len(times)}")

    rsum0, tr0, width = initial_guesses(eph, flux, times)
    B = bounds_for(rsum0, distortion)
    m = Model(times, flux, sig, period, t0, teff1, m1, q, distortion, args.ntri, l3=args.l3)
    m.l3_start, m.ecosw_start, m.esinw_start, m.q_start = args.l3, ecosw0, esinw0, q
    m.dt0_start = float(np.clip(ecl.get("phase1", 0.0) * period, -0.009, 0.009))
    if eccentric:
        log(f"  eccentric: secondary at phase {phase2:.4f} -> ecosw0={ecosw0:.4f}, esinw0={esinw0:.3f}")
    log(f"  sma={m.sma:.2f} Rsun, eclipse width {width:.3f} -> rsum0={rsum0:.3f}, teffratio0={tr0:.2f}")

    t_start = time.time()
    grid_cache = cdir / f"{name}_grid{tag}.json"
    if grid_cache.exists() and json.loads(grid_cache.read_text()).get("params") == PARAMS:
        g = json.loads(grid_cache.read_text())
        chi_grid, v_grid = g["chi2"], np.array(g["v"])
        top = [(c, np.array(v)) for c, v in g.get("top", [[chi_grid, g["v"]]])]
        log("  using cached grid result")
    else:
        (chi_grid, v_grid), top = grid_search(m, rsum0, tr0, B, log)
        if v_grid is not None:
            grid_cache.write_text(json.dumps({"params": PARAMS, "chi2": chi_grid, "v": list(map(float, v_grid)),
                                              "top": [[float(c), list(map(float, v))] for c, v in top]}))
    if v_grid is None:
        raise SystemExit("grid search produced no valid model")
    log(f"grid best chi2={chi_grid:.2f} ({time.time() - t_start:.0f}s)")
    opt_cache = cdir / f"{name}_localopt{tag}.json"
    if opt_cache.exists() and json.loads(opt_cache.read_text()).get("params") == PARAMS:
        o = json.loads(opt_cache.read_text())
        chi_pow, v_pow, nev = o["chi2"], np.array(o["v"]), o["evals"]
        log("  using cached local-optimisation result")
    else:
        # multi-start: local optimisation from each of the top distinct grid points
        chi_pow, v_pow, nev = np.inf, v_grid, 0
        for i, (c0, v0) in enumerate(top):
            c, v, n_ = powell(m, v0, B, args.powell_maxiter, log)
            nev += n_
            log(f"  start {i + 1}/{len(top)} (grid chi2 {c0:.1f}) -> chi2 {c:.2f}")
            if c < chi_pow:
                chi_pow, v_pow = c, v
        opt_cache.write_text(json.dumps({"params": PARAMS, "chi2": float(chi_pow), "v": list(map(float, v_pow)),
                                         "evals": nev}))
    log(f"local-opt (Nelder-Mead) chi2={chi_pow:.2f} after {nev} evals v={np.round(v_pow, 5).tolist()}")
    dof = max(len(times) - len(PARAMS) - 1, 1)
    scale = max(chi_pow / dof, 1.0)  # inflate errors so best fit has chi2_red ~ 1 for the posterior widths
    # emcee's red-blue move needs at least 2 x ndim walkers
    nwalkers = max(args.walkers, 2 * len(PARAMS) + 2)
    v_best, flat, chain, lp = run_emcee(m, v_pow, B, nwalkers, args.iters, scale, log)
    chi_mc, _ = m.evaluate(v_best)
    if chi_pow < chi_mc:
        v_best = v_pow
    np.save(cdir / f"{name}_emcee_chain{tag}.npy", chain)

    hi = Model(times, flux, sig, period, t0, teff1, m1, q, distortion, args.ntri_final, l3=args.l3)
    chi_hi, pred = hi.evaluate(v_best, label="final", keep=True)
    hi.b.save(str(cdir / f"{name}_phoebe_best{tag}.phoebe"))
    vb = dict(zip(PARAMS, v_best))
    incl, dt0, rsum, k, tr = (vb[n] for n in PARAMS[:5])
    l3_best = float(vb.get("l3", args.l3))
    ecc_best = math.hypot(vb.get("ecosw", 0.0), vb.get("esinw", 0.0))
    per0_best = math.degrees(math.atan2(vb.get("esinw", 0.0), vb.get("ecosw", 0.0))) % 360
    R1 = rsum / (1 + k) * hi.sma
    R2 = rsum * k / (1 + k) * hi.sma
    T2 = teff1 * tr
    phot = photometry_check(times, flux, pred, sig, t0 + dt0, period, len(PARAMS), phase2)
    bb = blackbody_check(name, R1, R2, teff1, T2, gaia_dist)
    post = {p: [float(np.percentile(flat[:, i], 16)), float(np.median(flat[:, i])), float(np.percentile(flat[:, i], 84))]
            for i, p in enumerate(PARAMS)}
    result = {
        "designation": name, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "period": period, "t0": t0 + dt0, "teff1": teff1, "teff2": T2, "q_assumed": q, "mass1_assumed": m1,
        "sma_rsun": hi.sma, "incl": incl, "rsum": rsum, "k": k, "teffratio": tr,
        "r1_frac": rsum / (1 + k), "r2_frac": rsum * k / (1 + k), "R1_rsun": R1, "R2_rsun": R2,
        "distortion": distortion, "l3_frac": l3_best, "l3_fitted": bool(args.fit_l3),
        "eccentric": bool(eccentric), "ecc": ecc_best, "per0_deg": per0_best,
        "ecosw": vb.get("ecosw", 0.0), "esinw": vb.get("esinw", 0.0), "params": list(PARAMS), "chi2_grid": chi_grid, "chi2_powell": chi_pow, "chi2_final_hires": chi_hi,
        "posterior_16_50_84": post, "photometry_check": phot,
        "blackbody_check": {k_: v_ for k_, v_ in bb.items() if k_ not in ("lam", "flux", "sigma", "model")},
        "runtime_s": time.time() - t_start,
    }
    result["q"] = float(vb.get("q", q))
    result["q_fitted"] = bool(fit_q)
    ratio = bb.get("distance_ratio_phot_over_gaia")
    if ratio and np.isfinite(ratio):
        # what the Gaia distance implies if the fitted geometry (R/a, T2/T1) is right and l3 = 0:
        # the system is 1/ratio^2 more luminous than assumed -> a scales by 1/ratio, total mass by 1/ratio^3
        mtot = m1 * (1 + result["q"])
        result["distance_scaling"] = {
            "luminosity_factor_vs_assumed": 1 / ratio**2,
            "implied_sma_rsun": hi.sma / ratio,
            "implied_total_mass_msun": mtot / ratio**3,
            "assumed_total_mass_msun": mtot,
            "note": "ratio << 1 means the system is brighter than an MS binary at the TIC Teff: Teff "
                    "underestimated (reddening), evolved components, or third light",
        }
    result["passes_all"] = bool(phot["pass_fit"] and phot["pass_no_eclipse_bias"] and bb.get("pass_photosphere")
                                and bb.get("pass_distance") is not False and bb.get("pass_excess_survives"))
    (cdir / f"{name}_phoebe_result{tag}.json").write_text(json.dumps(result, indent=2, default=float) + "\n")
    plot_fit(cdir / f"{name}_phoebe_fit{tag}.png", name, times, flux, sig, pred, t0 + dt0, period,
             f"chi2_red={phot['chi2_red']:.2f} rms={phot['rms']:.5f} i={incl:.2f} r1={R1:.2f} r2={R2:.2f} Rsun "
             f"T2/T1={tr:.3f}")
    if "lam" in bb:
        plot_sed(cdir / f"{name}_phoebe_sed_check{tag}.png", name, bb)
    log(json.dumps({k_: result[k_] for k_ in ("incl", "r1_frac", "r2_frac", "teffratio", "chi2_final_hires",
                                               "passes_all")}, default=float))
    log(f"photometry: {json.dumps(phot, default=float)}")
    log(f"blackbody: {json.dumps(result['blackbody_check'], default=float)}")


if __name__ == "__main__":
    main()
