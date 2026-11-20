"""Free-free slope range vs a cool dust blackbody, normalised at W3, with J0719's W3/W4 excess.

F_nu ∝ nu^alpha: alpha = -0.1 (optically thin) to +0.6 (ionised stellar wind) is the range the excess
diagnostic and sanity_checks.py test; +2 (fully optically thick, Rayleigh-Jeans) is the physical limit.
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
sys.stdout.reconfigure(encoding="utf-8")
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
ORANGE, AQUA, YELLOW, GREY = "#eb6834", "#1baf7a", "#eda100", "#a9a8a2"
plt.rcParams.update({"font.family": "Arial", "font.size": 11, "mathtext.fontset": "cm", "axes.edgecolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})

d = pd.read_csv(ROOT / "results" / "test" / "wise_excess" / "J071951.40-240400.6_wise_excess.csv").set_index("band")
L3, L4 = d.loc["W3", "lambda_um"], d.loc["W4", "lambda_um"]
ex3, ex4 = d.loc["W3", "excess_jy"], d.loc["W4", "excess_jy"]
e3, e4 = d.loc["W3", "flux_sigma_jy"], d.loc["W4", "flux_sigma_jy"]


def pl(lam, alpha):  # F_nu ∝ nu^alpha = lambda^-alpha, normalised at W3
    return (lam / L3) ** (-alpha)


def bb(lam, T):
    x = 14387.77 / (lam * T)
    f = 1 / (lam ** 3 * np.expm1(x))
    return f / (1 / (L3 ** 3 * np.expm1(14387.77 / (L3 * T))))


g = np.logspace(np.log10(3), np.log10(60), 400)
fig, ax = plt.subplots(figsize=(8.5, 5.6))
ax.fill_between(g, pl(g, 0.6), pl(g, -0.1), color=ORANGE, alpha=0.18, lw=0, zorder=1)
ax.plot(g, pl(g, -0.1), color=ORANGE, lw=2, label=r"free-free, optically thin: $\alpha = -0.1$")
ax.plot(g, pl(g, 0.6), color=YELLOW, lw=2, label=r"free-free, ionised stellar wind: $\alpha = +0.6$")
ax.plot(g, pl(g, 2.0), color=GREY, lw=1.5, ls=(0, (4, 3)), label=r"free-free limit, fully opaque: $\alpha = +2$")
ax.plot(g, bb(g, 158.2), color=AQUA, lw=2.5, label="dust blackbody, 158 K")
ax.plot(g, pl(g, -2.8), color=INK2, lw=1, ls=":", label=r"power law J0719 would need: $\alpha = -2.8$")
r = ex4 / ex3
rerr = r * np.hypot(e3 / ex3, e4 / ex4)
ax.errorbar([L3, L4], [1, r], yerr=[e3 / ex3, rerr], fmt="o", ms=9, color=INK, mec=SURF, mew=1.5, zorder=6,
            elinewidth=1.5, label="J0719 excess (W3, W4)")
ax.text(4.6, 0.5, "shaded: slopes tested\nfor ionised gas", fontsize=9, color=INK2, ha="center", va="top")
ax.annotate(f"W4/W3 excess = {r:.1f}", (L4, r), xytext=(26, 2.2), fontsize=10, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.7))
for lam, n in ((3.35, "W1"), (4.6, "W2"), (L3, "W3"), (L4, "W4")):
    ax.axvline(lam, color=GRID, lw=1.2, zorder=0)
    ax.text(lam, 30, n, ha="center", fontsize=9, color=INK2)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(3, 60); ax.set_ylim(0.01, 40)
ax.set_xticks([3, 5, 10, 20, 50]); ax.set_xticklabels(["3", "5", "10", "20", "50"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Wavelength (µm)"); ax.set_ylabel(r"Excess flux $F_\nu$, relative to W3")
ax.set_title("Ionised gas is too flat to make J0719's 12 → 22 µm jump", loc="left", fontsize=13, fontweight="bold",
             color=INK, pad=20)
ax.text(0, 1.015, r"Spectral shape of the excess only (stars subtracted); every curve set to 1 at 11.6 µm; $F_\nu \propto \nu^{\alpha}$",
        transform=ax.transAxes, fontsize=9.5, color=INK2)
ax.legend(fontsize=9, loc="lower right")
fig.tight_layout(); fig.savefig(OUT / "freefree_slope_range.png", dpi=200)
print(f"W4/W3 excess {r:.2f} ± {rerr:.2f}; free-free predicts {pl(L4, -0.1):.2f} (α=-0.1), {pl(L4, 0.6):.2f} (α=0.6), "
      f"{pl(L4, 2):.2f} (α=2); dust 158 K {bb(L4, 158.2):.2f}")
