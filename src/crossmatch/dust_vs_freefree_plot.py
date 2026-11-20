"""J0719: WISE photometry with the dust and free-free models of src/wise_excess_diagnostic.py.

Each model = Rayleigh-Jeans photosphere (A / lambda^2) + one excess component, both amplitudes fitted
by weighted least squares to W1-W4. Reproduces the chi2 in results/test/wise_excess/*_summary.md.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "excess_equations"
OUT.mkdir(parents=True, exist_ok=True)
sys.stdout.reconfigure(encoding="utf-8")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
plt.rcParams.update({"font.family": "Arial", "font.size": 11, "mathtext.fontset": "cm", "axes.edgecolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})

d = pd.read_csv(ROOT / "results" / "test" / "wise_excess" / "J071951.40-240400.6_wise_excess.csv")
lam, F, s = d.lambda_um.values, 1e3 * d.flux_jy.values, 1e3 * d.flux_sigma_jy.values  # mJy


def bnu(lam_um, T):
    x = 14387.77 / (lam_um * T)
    return 1 / (lam_um ** 3 * np.expm1(x))


def fit(basis2):
    X = np.vstack([lam ** -2.0, basis2(lam)]).T
    w = 1 / s
    coef, *_ = np.linalg.lstsq(X * w[:, None], F * w, rcond=None)
    chi2 = float(np.sum(((F - X @ coef) / s) ** 2))
    return coef, chi2


models = {
    "dust": (lambda l: bnu(l, 158.2), AQUA, "dust blackbody, 158 K"),
    "ff": (lambda l: l ** 0.1, ORANGE, r"optically thin free-free, $F_\nu \propto \nu^{-0.1}$"),
    "wind": (lambda l: l ** -0.6, YELLOW, r"stellar-wind free-free, $F_\nu \propto \nu^{+0.6}$"),
}
g = np.logspace(np.log10(2.5), np.log10(30), 300)
fig, ax = plt.subplots(figsize=(8.5, 5.6))
phot_only = np.sum(F[:2] * lam[:2] ** -2 / s[:2] ** 2) / np.sum(lam[:2] ** -4 / s[:2] ** 2)
ax.plot(g, phot_only / g ** 2, color=BLUE, lw=2, ls=(0, (4, 3)), label=r"stars only: $A/\lambda^2$ fitted to W1, W2")
res = {}
for k, (b, c, lab) in models.items():
    coef, chi2 = fit(b)
    res[k] = chi2
    ax.plot(g, coef[0] / g ** 2 + coef[1] * b(g), color=c, lw=2, label=f"stars + {lab}:  $\\chi^2$ = {chi2:.2f}" if chi2 < 1
            else f"stars + {lab}:  $\\chi^2$ = {chi2:.0f}")
ax.errorbar(lam, F, yerr=s, fmt="o", ms=8, color=INK, mfc=INK, mec=SURF, mew=1.5, zorder=5, capsize=0, elinewidth=1.5)
for l, f, n in zip(lam, F, ["W1", "W2", "W3", "W4"]):
    ax.annotate(n, (l, f), xytext=(7, 6), textcoords="offset points", fontsize=10, color=INK)
ax.annotate("W4 = 10.7× the stars\n(13σ)", (22.09, F[3]), xytext=(11.5, 28), fontsize=10, color=INK2,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([3, 5, 10, 20, 30]); ax.set_xticklabels(["3", "5", "10", "20", "30"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlim(2.5, 30); ax.set_ylim(0.8, 120)
ax.set_xlabel("Wavelength (µm)"); ax.set_ylabel(r"Flux density $F_\nu$ (mJy)")
ax.set_title("J0719: cool dust fits the WISE excess; ionised gas cannot", loc="left", fontsize=13, fontweight="bold",
             color=INK, pad=20)
ax.text(0, 1.015, "Each model = stellar photosphere + one excess component, fitted to all four WISE bands",
        transform=ax.transAxes, fontsize=9.5, color=INK2)
ax.legend(fontsize=9.5, loc="lower left")
fig.tight_layout(); fig.savefig(OUT / "J0719_dust_vs_freefree.png", dpi=200)
print({k: round(v, 2) for k, v in res.items()})
