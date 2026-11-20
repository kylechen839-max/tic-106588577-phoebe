"""Reddened-blackbody SED fit (Teff and A_V free) from Gaia G/BP/RP + 2MASS + WISE W1/W2.

  python src/crossmatch/sed_extinction_fit.py DESIGNATION OUTDIR

For stars with no catalogue Teff, where our J/H/Ks/W1 blackbody (no extinction) may mistake reddening
for a cool photosphere. Model: F_nu = s * B_nu(lambda_eff, T) * 10^(-0.4 A_V k_lambda), with
k = A_lambda/A_V from Wang & Chen (2019) for the Gaia, 2MASS and WISE bands. W3/W4 are left out (excess).
A blackbody at each band's effective wavelength is crude for the broad Gaia bands, so a 0.05 mag error
floor is used and the result is a temperature range, not a spectral type.
"""
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.ipac.irsa import Irsa

ROOT = Path(__file__).resolve().parents[2]
sys.stdout.reconfigure(encoding="utf-8")
NAME, OUT = sys.argv[1], Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
row = s[s.designation == NAME].iloc[0]
c = SkyCoord(row.ra, row.dec, unit="deg")

#        lambda_eff (um), Vega zero point (Jy), A_lambda/A_V (Wang & Chen 2019)
BANDS = {"BP": (0.511, 3552.0, 1.002), "G": (0.622, 3229.0, 0.789), "RP": (0.777, 2555.0, 0.589),
         "J": (1.235, 1594.0, 0.243), "H": (1.662, 1024.0, 0.131), "K": (2.159, 666.7, 0.078),
         "W1": (3.35, 309.54, 0.039), "W2": (4.60, 171.79, 0.026)}
FLOOR = 0.05

g = Gaia.launch_job(
    "SELECT source_id, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, phot_g_mean_flux_over_error, "
    "phot_bp_mean_flux_over_error, phot_rp_mean_flux_over_error, phot_bp_rp_excess_factor, parallax, "
    "parallax_error, ruwe, DISTANCE(POINT(ra, dec), POINT(%f, %f))*3600 AS sep FROM gaiadr3.gaia_source "
    "WHERE 1=CONTAINS(POINT(ra, dec), CIRCLE(%f, %f, %f)) ORDER BY sep" %
    (c.ra.deg, c.dec.deg, c.ra.deg, c.dec.deg, 3 / 3600)).get_results().to_pandas().iloc[0]
a = Irsa.query_region(c, catalog="allwise_p3as_psd", spatial="Cone", radius=3 * u.arcsec).to_pandas().iloc[0]
mag = {"BP": (g.phot_bp_mean_mag, 1.0857 / g.phot_bp_mean_flux_over_error),
       "G": (g.phot_g_mean_mag, 1.0857 / g.phot_g_mean_flux_over_error),
       "RP": (g.phot_rp_mean_mag, 1.0857 / g.phot_rp_mean_flux_over_error),
       "J": (a.j_m_2mass, a.j_msig_2mass), "H": (a.h_m_2mass, a.h_msig_2mass), "K": (a.k_m_2mass, a.k_msig_2mass),
       "W1": (a.w1mpro, a.w1sigmpro), "W2": (a.w2mpro, a.w2sigmpro)}
bands = [b for b in BANDS if np.isfinite(mag[b][0])]
lam = np.array([BANDS[b][0] for b in bands])
k = np.array([BANDS[b][2] for b in bands])
F = np.array([BANDS[b][1] * 10 ** (-0.4 * mag[b][0]) for b in bands])
em = np.array([max(mag[b][1] if np.isfinite(mag[b][1]) else 0.1, FLOOR) for b in bands])
sF = F * 0.4 * np.log(10) * em
print("photometry:", {b: round(float(mag[b][0]), 3) for b in bands})
print(f"Gaia: parallax {g.parallax:.3f} ± {g.parallax_error:.3f} mas, RUWE {g.ruwe:.1f}, BP/RP excess factor {g.phot_bp_rp_excess_factor:.2f}")


def bnu(l, T):
    return 1 / (l ** 3 * np.expm1(14387.77 / (l * T)))


Ts = np.exp(np.linspace(np.log(2500), np.log(45000), 260))
Avs = np.linspace(0, 15, 301)
chi = np.full((len(Avs), len(Ts)), np.nan)
for i, Av in enumerate(Avs):
    ext = 10 ** (-0.4 * Av * k)
    for j, T in enumerate(Ts):
        m = bnu(lam, T) * ext
        w = sF ** -2
        sc = np.sum(w * m * F) / np.sum(w * m * m)
        chi[i, j] = np.sum(((F - sc * m) / sF) ** 2)
i0, j0 = np.unravel_index(np.nanargmin(chi), chi.shape)
cmin = chi[i0, j0]
dof = len(bands) - 3
ok = chi <= cmin + 2.30  # 68% for 2 parameters
ok95 = chi <= cmin + 6.18
Tbest, Avbest = Ts[j0], Avs[i0]
Tr = (Ts[ok.any(0)].min(), Ts[ok.any(0)].max())
Tr95 = (Ts[ok95.any(0)].min(), Ts[ok95.any(0)].max())
Avr = (Avs[ok.any(1)].min(), Avs[ok.any(1)].max())
noext = chi[0].min(), Ts[np.nanargmin(chi[0])]
print(f"best: T = {Tbest:.0f} K, A_V = {Avbest:.2f}, chi2 = {cmin:.2f} (dof {dof})")
print(f"68% region: T {Tr[0]:.0f}-{Tr[1]:.0f} K, A_V {Avr[0]:.2f}-{Avr[1]:.2f};  95%: T {Tr95[0]:.0f}-{Tr95[1]:.0f} K")
print(f"no extinction (A_V = 0): best T = {noext[1]:.0f} K, chi2 = {noext[0]:.1f}")
for Tq in (3350, 5000, 8000, 12000, 20000, 30000):
    jj = np.argmin(np.abs(Ts - Tq))
    ii = np.nanargmin(chi[:, jj])
    print(f"  at T = {Tq:>5} K: best A_V = {Avs[ii]:.2f}, chi2 = {chi[ii, jj]:.1f}")

fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]
im = ax.pcolormesh(Ts, Avs, np.clip(chi - cmin, 0, 60), shading="auto", cmap="viridis_r")
ax.contour(Ts, Avs, chi - cmin, levels=[2.30, 6.18], colors=["w", "w"], linestyles=["-", "--"], linewidths=1)
ax.plot(Tbest, Avbest, "r*", ms=12)
ax.set_xscale("log"); ax.set_xlabel("blackbody temperature (K)"); ax.set_ylabel("A_V (mag)")
ax.set_title(f"{NAME}: Δχ² (solid 68%, dashed 95%)", fontsize=10)
fig.colorbar(im, ax=ax, label="Δχ² (clipped at 60)")
ax = axs[1]
gl = np.logspace(np.log10(0.35), np.log10(25), 300)
kk = np.interp(np.log(gl), np.log(lam), np.log(k))
ax.errorbar(lam, F * 1e3, yerr=sF * 1e3, fmt="o", color="k", label="Gaia + 2MASS + WISE W1/W2")
for b3, l3 in (("W3", 11.56), ("W4", 22.09)):
    m3 = a[f"{b3.lower()}mpro"]
    ax.plot(l3, {"W3": 31.674, "W4": 8.363}[b3] * 10 ** (-0.4 * m3) * 1e3, "s", mfc="none", mec="k", label=None)
for T, Av, col, lab in ((Tbest, Avbest, "C3", f"best: {Tbest:.0f} K, A_V {Avbest:.1f}"), (noext[1], 0, "C0", f"no extinction: {noext[1]:.0f} K")):
    m = bnu(lam, T) * 10 ** (-0.4 * Av * k); w = sF ** -2; sc = np.sum(w * m * F) / np.sum(w * m * m)
    ext = 10 ** (-0.4 * Av * np.exp(kk)) if Av else 1
    ax.plot(gl, sc * bnu(gl, T) * ext * 1e3, color=col, label=lab)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("wavelength (µm)"); ax.set_ylabel("F_ν (mJy)")
ax.set_title("SED (open squares: W3/W4, not fitted)", fontsize=10); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_sed_extinction_fit.png", dpi=110)
json.dump(dict(designation=NAME, T_best=float(Tbest), Av_best=float(Avbest), chi2=float(cmin), dof=dof,
               T68=[float(x) for x in Tr], T95=[float(x) for x in Tr95], Av68=[float(x) for x in Avr],
               T_noext=float(noext[1]), chi2_noext=float(noext[0]), parallax=float(g.parallax),
               parallax_error=float(g.parallax_error), ruwe=float(g.ruwe)),
          open(OUT / f"{NAME}_sed_extinction_fit.json", "w"), indent=1)
