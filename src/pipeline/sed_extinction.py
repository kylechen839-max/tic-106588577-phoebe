"""Extinction-aware binary-blackbody SED check for finished PHOEBE results.

The post-fit check in run_phoebe_candidate.py fits an unreddened binary blackbody to J, H, Ks, W1 only.
That cannot tell a hot reddened star from a cool unreddened one, which is why the TIC Teffs (reddening-
biased low) left every system "overluminous" at its Gaia distance. This script adds Gaia G/BP/RP and a
free extinction A_V (Wang & Chen 2019 law) and asks two questions per result:

1. Reddened photosphere at the PHOEBE Teff: fit scale + A_V to G, BP, RP, J, H, Ks, W1. Report chi2,
   A_V (compare with Gaia GSP-Phot A0 and TIC E(B-V)), d_phot/d_Gaia, and the W3/W4 excess against the
   reddened model.
2. Teff at the Gaia distance: scan T1 (T2/T1 and R/a held at the fitted values, sma rescaled with the
   main-sequence mass for each T1), with the flux scale fixed by the Gaia distance and A_V free.
   The best T1 is the Teff to pass to `run_phoebe_candidate.py --teff` for a self-consistent refit.

Gaia passbands are treated as monochromatic at their effective wavelengths and the stars as blackbodies,
so expect ~5% systematics; a 0.05 mag error floor is used for the optical bands, 0.03 mag for 2MASS/WISE.
Third light (l3) is ignored. Writes `sed_extinction` into each result JSON (it does not change
`passes_all`) and a markdown table to outputs/candidates/sed_extinction_summary.md.

Usage: ~/phoebe-env/bin/python src/pipeline/sed_extinction.py [result.json ...]
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_phoebe_candidate as r  # noqa: E402

GAIA = r.ROOT / "outputs" / "crossmatch" / "dd1_support_cache" / "gaia.csv"
TIC = r.ROOT / "outputs" / "crossmatch" / "dd1_support_cache" / "tic.csv"
OUT_MD = r.ROOT / "outputs" / "candidates" / "sed_extinction_summary.md"

# Vega zero points (Jy) and effective wavelengths (micron) for the Gaia DR3 bands (SVO filter service)
ZP = dict(r.ZP, G=3228.75, BP=3552.01, RP=2554.95)
LAM = dict(r.LAM, G=0.5822, BP=0.5036, RP=0.7620)
# A_lambda / A_V, Wang & Chen (2019)
EXT = {"BP": 1.002, "G": 0.789, "RP": 0.589, "J": 0.243, "H": 0.131, "K": 0.078,
       "W1": 0.039, "W2": 0.026, "W3": 0.040, "W4": 0.020}
PHOT = ["BP", "G", "RP", "J", "H", "K", "W1"]
ALL = PHOT + ["W2", "W3", "W4"]
AV_GRID = np.arange(0.0, 8.0001, 0.01)


def first_row(path, name):
    rows = sorted(r.read_rows(path, name), key=lambda x: float(x["angDist"]))
    return rows[0] if rows else None


def photometry(name):
    a = first_row(r.ALLWISE, name) or first_row(r.EXTRA_ALLWISE, name)
    g = first_row(GAIA, name)
    if a is None or g is None:
        return None, None
    mags = {}
    for b, c in zip(["J", "H", "K", "W1", "W2", "W3", "W4"], ["Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag"]):
        mags[b] = (r.fnum(a[c]), max(r.fnum(a["e_" + c]) or 0.1, 0.03))
    for b, c in (("G", "Gmag"), ("BP", "BPmag"), ("RP", "RPmag")):
        mags[b] = (r.fnum(g[c]), max(r.fnum(g["e_" + c]) or 0.05, 0.05))
    bands = [b for b in ALL if mags[b][0] is not None]
    f = np.array([ZP[b] * 10 ** (-0.4 * mags[b][0]) for b in bands])
    s = np.array([f[i] * 0.4 * math.log(10) * mags[b][1] for i, b in enumerate(bands)])
    t = first_row(TIC, name)
    priors = {"gaia_A0": r.fnum(g.get("A0")), "gaia_teff": r.fnum(g.get("Teff")),
              "tic_AV": 3.1 * r.fnum(t["E(B-V)"]) if t and r.fnum(t.get("E(B-V)")) is not None else None}
    return (bands, f, s), priors


def shape(bands, R1, R2, T1, T2):
    lam = np.array([LAM[b] for b in bands])
    return R1**2 * r.planck(T1, lam) + R2**2 * r.planck(T2, lam)


def fit_free_scale(bands, f, s, sh):
    """Best (chi2, A_V, scale) over the A_V grid, photosphere bands only."""
    ip = np.array([b in PHOT for b in bands])
    k = np.array([EXT[b] for b in bands])
    best = (np.inf, 0.0, 0.0)
    for av in AV_GRID:
        m = sh * 10 ** (-0.4 * av * k)
        w = s[ip] ** -2
        A = np.sum(w * m[ip] * f[ip]) / np.sum(w * m[ip] ** 2)
        c = float(np.sum(((f[ip] - A * m[ip]) / s[ip]) ** 2))
        if c < best[0]:
            best = (c, float(av), float(A))
    return best


def fit_fixed_scale(bands, f, s, sh, A):
    ip = np.array([b in PHOT for b in bands])
    k = np.array([EXT[b] for b in bands])
    chis = [float(np.sum(((f[ip] - A * sh[ip] * 10 ** (-0.4 * av * k[ip])) / s[ip]) ** 2)) for av in AV_GRID]
    i = int(np.argmin(chis))
    return chis[i], float(AV_GRID[i])


def analyse(res, phot, priors, gaia_dist):
    bands, f, s = phot
    R1, R2, T1, T2 = res["R1_rsun"], res["R2_rsun"], res["teff1"], res["teff2"]
    sh = shape(bands, R1, R2, T1, T2)
    chi2, av, A = fit_free_scale(bands, f, s, sh)
    k = np.array([EXT[b] for b in bands])
    model = A * sh * 10 ** (-0.4 * av * k)
    d_phot = r.RSUN_M * math.sqrt(math.pi * 1e26 / A) / r.PC_M if A > 0 else np.nan
    exc = {b: float((f[i] - model[i]) / s[i]) for i, b in enumerate(bands) if b in ("W3", "W4")}
    out = {"bands": bands, "chi2_phot": chi2, "n_phot": int(sum(b in PHOT for b in bands)), "A_V": av,
           "d_phot_pc": d_phot, "d_ratio": d_phot / gaia_dist if gaia_dist else np.nan,
           "excess_sigma": exc, **priors}
    # Teff at the Gaia distance: scale fixed, A_V free, geometry rescaled with MS mass
    if gaia_dist and np.isfinite(gaia_dist):
        q = res.get("q", res["q_assumed"])
        scan = []
        for T in np.arange(4000, 40001, 100.0):
            m1, _ = r.ms_mass_radius(T)
            a = r.kepler_sma(res["period"], m1 * (1 + q))
            sh_t = shape(bands, res["r1_frac"] * a, res["r2_frac"] * a, T, T * res["teffratio"])
            A_d = math.pi * 1e26 * (r.RSUN_M / (gaia_dist * r.PC_M)) ** 2
            c, av_t = fit_fixed_scale(bands, f, s, sh_t, A_d)
            scan.append((c, T, av_t))
        c, T, av_t = min(scan)
        out.update({"teff_at_gaia_dist": T, "A_V_at_gaia_dist": av_t, "chi2_at_gaia_dist": c})
    return out


def main():
    paths = [Path(p) for p in sys.argv[1:]] or sorted((r.ROOT / "outputs" / "candidates").glob("*/*_phoebe_result*.json"))
    lines = ["| result | T1 used | χ²_red LC | A_V fit | Gaia A0 | TIC A_V | χ² BP..W1 (7 bands) | d_phot/d_Gaia (reddened) | W3σ / W4σ | Teff at Gaia d | A_V there | χ² there |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|"]
    for p in paths:
        res = json.loads(p.read_text())
        name = res["designation"]
        phot, priors = photometry(name)
        if phot is None:
            print(f"{p.name}: no AllWISE/Gaia row")
            continue
        rows = r.read_rows(r.SCORED, name) or (r.read_rows(r.EXTRA_SCORED, name) if r.EXTRA_SCORED.exists() else [])
        row = rows[0] if rows else {}
        plx = r.fnum(row.get("plx"))
        gd = r.fnum(row.get("gaia_dist")) or (1000.0 / plx if plx and plx > 0 else None)
        sx = analyse(res, phot, priors, gd)
        res["sed_extinction"] = sx
        p.write_text(json.dumps(res, indent=2, default=float) + "\n")
        tag = p.stem.split("_phoebe_result")[1] or "(nominal)"
        fmt = lambda v, f="{:.2f}": f.format(v) if v is not None and np.isfinite(v) else "–"  # noqa: E731
        ex = sx["excess_sigma"]
        lines.append(f"| {name} {tag} | {res['teff1']:.0f} | {res['photometry_check']['chi2_red']:.1f} | {sx['A_V']:.2f} | "
                     f"{fmt(priors['gaia_A0'])} | {fmt(priors['tic_AV'])} | {sx['chi2_phot']:.1f} | {fmt(sx['d_ratio'])} | "
                     f"{fmt(ex.get('W3'), '{:.1f}')} / {fmt(ex.get('W4'), '{:.1f}')} | {fmt(sx.get('teff_at_gaia_dist'), '{:.0f}')} | "
                     f"{fmt(sx.get('A_V_at_gaia_dist'))} | {fmt(sx.get('chi2_at_gaia_dist'), '{:.1f}')} |")
        print(lines[-1])
    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
