"""Population plots for all 800 DD1.0 x eclipsing-binary cross-matches.

Reads outputs/crossmatch/dd1_candidates_scored.csv and writes PNGs to
outputs/crossmatch/population_plots/. Colours are the dataviz skill's validated default
categorical slots (fixed order), with tier C / "other" in neutral grey.

Usage: python src/crossmatch/population_plots.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "population_plots"
OUT.mkdir(parents=True, exist_ok=True)

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
GREY = "#a9a8a2"
TIER = {"A": (SLOT[0], 70, 3), "B": (SLOT[1], 22, 2), "C": (GREY, 9, 1)}  # colour, marker size, zorder
CLASSES = [  # (label, SIMBAD otypes)
    ("Eclipsing / spectroscopic binary", {"EclBin", "SB*", "**", "EllipVar", "RSCVnV*"}),
    ("Evolved / long-period variable", {"LongPeriodV*", "LongPeriodV*_Candidate", "Mira", "RGB*", "AGB*", "C*",
                                        "S*", "OH/IR*", "post-AGB*", "Supergiant", "RedSG*", "semi-regV*"}),
    ("Young stellar object", {"TTauri*", "OrionV*", "YSO", "YSO_Candidate", "Ae*", "HerbigHaro", "TTauri*_Candidate"}),
    ("Hot emission-line star", {"Be*", "EmLine*", "WolfRayet*", "Symbiotic*", "BlueSG*", "Em*", "PN"}),
]
CLASS_COL = SLOT[:4]

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.axisbelow": True, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False})

s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
s["cls"] = "Other / unclassified"
for lab, types in CLASSES:
    s.loc[s.simbad_otype.isin(types), "cls"] = lab
s["MG"] = np.where(s.plx_snr > 5, s.Gmag + 5 * np.log10(s.plx.clip(lower=1e-3) / 100), np.nan)
g = SkyCoord(s.ra.values * u.deg, s.dec.values * u.deg).galactic
s["l"] = g.l.wrap_at(180 * u.deg).deg
s["b"] = g.b.deg
N = len(s)


def by_tier(ax, x, y, **kw):
    for t in ("C", "B", "A"):
        m = (s.tier == t) & s[x].notna() & s[y].notna()
        c, ms, z = TIER[t]
        ax.scatter(s.loc[m, x], s.loc[m, y], s=ms, c=c, edgecolors=SURF, linewidths=0.6 if t == "A" else 0.3,
                   zorder=z, label=f"Tier {t} ({(s.tier == t).sum()})", **kw)


def by_class(ax, x, y):
    order = ["Other / unclassified"] + [c for c, _ in CLASSES]
    cols = dict(zip(order, [GREY] + CLASS_COL))
    for c in order:
        m = (s.cls == c) & s[x].notna() & s[y].notna()
        ax.scatter(s.loc[m, x], s.loc[m, y], s=14, c=cols[c], edgecolors=SURF, linewidths=0.3,
                   zorder=1 if c.startswith("Other") else 2, label=f"{c} ({(s.cls == c).sum()})")
    a = (s.tier == "A") & s[x].notna() & s[y].notna()
    ax.scatter(s.loc[a, x], s.loc[a, y], s=80, facecolors="none", edgecolors=INK, linewidths=1.2, zorder=4,
               label="Tier A (12)")


def title(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12, color=INK, fontweight="bold", pad=18)
    ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=9, color=INK2, va="bottom")


def save(fig, name):
    fig.tight_layout(); fig.savefig(OUT / name, dpi=200); plt.close(fig)


# 1. Gaia colour-magnitude diagram, by SIMBAD class
fig, ax = plt.subplots(figsize=(7.5, 6))
by_class(ax, "bp_rp", "MG")
ax.invert_yaxis(); ax.set_xlabel("Gaia BP − RP (mag)"); ax.set_ylabel("Absolute G (mag)")
title(ax, "What kinds of stars are the 800 matches?",
      f"Gaia CMD, parallax S/N > 5 ({s.MG.notna().sum()} of {N}); colour = SIMBAD type; circles = tier A")
ax.legend(fontsize=8, loc="lower left", markerscale=1.3)
save(fig, "01_cmd_by_simbad_class.png")

# 2. Same CMD, by tier
fig, ax = plt.subplots(figsize=(7.5, 6))
by_tier(ax, "bp_rp", "MG")
ax.invert_yaxis(); ax.set_xlabel("Gaia BP − RP (mag)"); ax.set_ylabel("Absolute G (mag)")
title(ax, "The checks keep hot, blue main-sequence stars",
      "Gaia CMD of all matches, coloured by candidate tier")
ax.legend(fontsize=9, loc="lower left", markerscale=1.2)
save(fig, "02_cmd_by_tier.png")

# 3. Galactic sky map
fig, ax = plt.subplots(figsize=(10, 4.6))
by_tier(ax, "l", "b")
ax.set_xlim(180, -180); ax.set_ylim(-90, 90)
ax.set_xlabel("Galactic longitude (°)"); ax.set_ylabel("Galactic latitude (°)")
title(ax, "Candidates crowd the Galactic plane",
      f"{(np.abs(s.b) < 10).mean():.0%} of matches lie within 10° of the plane (tier A: {(np.abs(s.b[s.tier == 'A']) < 10).mean():.0%})")
for lab, (l0, b0) in {"LMC": (-79.5, -32.9), "SMC": (-57.2, -44.3)}.items():
    ax.annotate(lab, (l0, b0), xytext=(l0 - 12, b0 - 10), fontsize=8, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.legend(fontsize=9, loc="lower left", ncol=3, markerscale=1.2)
save(fig, "03_galactic_map.png")

# 4. Dust temperature vs fractional IR luminosity
fig, ax = plt.subplots(figsize=(7.5, 6))
m = s.f_ir > 0
s["log_fir"] = np.where(m, np.log10(s.f_ir.where(m)), np.nan)
by_class(ax, "dust_temp", "log_fir")
ax.axhspan(np.log10(1e-6), np.log10(1e-2), color=SLOT[2], alpha=0.07, zorder=0)
ax.text(1450, np.log10(3e-3), "debris-disk range\n($L_{IR}/L_* < 10^{-2}$)", ha="right", va="top", fontsize=8, color=INK2)
deb = (s.f_ir > 0) & (s.f_ir < 1e-2)
ax.set_xlabel("Dust blackbody temperature (K)"); ax.set_ylabel("log fractional IR luminosity  $L_{IR}/L_*$")
title(ax, f"Tier A sits in the debris-disk regime ({(deb & (s.tier == 'A')).sum()} of 12); most matches do not",
      f"{deb.mean():.0%} of all 800 have $L_{{IR}}/L_* < 10^{{-2}}$; one-temperature dust fit to W2–W4 (40 K and 1500 K are grid limits)")
ax.legend(fontsize=8, loc="upper right", markerscale=1.3)
save(fig, "04_dust_temp_vs_fir.png")

# 5. Period distribution by tier
fig, ax = plt.subplots(figsize=(7.5, 4.6))
bins = np.logspace(-1.3, 3.3, 37)
p = s.period_d.where(s.period_d > 0)
ax.hist([p[s.tier == "C"].dropna(), p[s.tier == "B"].dropna(), p[s.tier == "A"].dropna()], bins=bins,
        stacked=True, color=[GREY, SLOT[1], SLOT[0]], edgecolor=SURF, linewidth=1,
        label=[f"Tier C ({p[s.tier == 'C'].notna().sum()})", f"Tier B ({p[s.tier == 'B'].notna().sum()})",
               f"Tier A ({p[s.tier == 'A'].notna().sum()})"])
ax.set_xscale("log"); ax.set_xlabel("Catalogue period (days)"); ax.set_ylabel("Number of matches")
lpv = (s.cls == "Evolved / long-period variable") & (p > 0)
title(ax, "Tier A periods span 1–17 days; tier C carries the long-period tail",
      f"{p.notna().sum()} matches with a catalogue period; tier A median {p[s.tier == 'A'].median():.1f} d, all matches {p.median():.1f} d")
ax.legend(fontsize=9)
save(fig, "05_period_histogram.png")

# 6. W4 excess vs photosphere temperature
fig, ax = plt.subplots(figsize=(7.5, 6))
s["log_w4"] = np.log10(s.W4_ratio.where(s.W4_ratio > 0))
by_tier(ax, "bb_teff", "log_w4")
ax.axhline(0, color=INK2, lw=0.8, zorder=0)
ax.set_xscale("log"); ax.set_xlabel("Photosphere blackbody temperature (K)")
ax.set_ylabel("log W4 flux / photosphere")
title(ax, "Tier A: hot stars with modest 22 µm excess",
      "W4 excess against the J/H/Ks/W1 photosphere fit; 0 = no excess; the column at 2500 K is the fit's lower limit")
ax.legend(fontsize=9, loc="upper right", markerscale=1.2)
save(fig, "06_w4_excess_vs_teff.png")

# 7. Funnel of checks (horizontal bars)
c = {k: s[k].astype(bool) if k in s else None for k in ("w3w4_quality_ok", "c_ir_excess", "c_photosphere",
                                                        "c_unconfused", "c_nearby", "c_not_evolved_or_wind")}
steps = [("All DD1.0 × EB matches", N)]
mask = np.ones(N, bool)
for k, lab in (("w3w4_quality_ok", "W3/W4 photometry quality"), ("c_ir_excess", "Significant IR excess"),
               ("c_unconfused", "No source confusion"), ("c_not_evolved_or_wind", "Not an evolved / wind star"),
               ("c_nearby", "Distance < 3 kpc")):
    if c[k] is not None:
        mask &= c[k].values
        steps.append((lab, int(mask.sum())))
steps += [("Tier A + B candidates", int(s.tier.isin(["A", "B"]).sum())), ("Tier A (high priority)", int((s.tier == "A").sum()))]
fig, ax = plt.subplots(figsize=(8, 4.2))
labs, vals = zip(*steps)
y = np.arange(len(vals))[::-1]
ax.barh(y, vals, height=0.62, color=[SLOT[0] if i == len(vals) - 1 else (SLOT[1] if i == len(vals) - 2 else GREY) for i in range(len(vals))],
        edgecolor=SURF, linewidth=2)
for yi, v in zip(y, vals):
    ax.text(v + N * 0.01, yi, f"{v}", va="center", fontsize=10, color=INK)
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=9); ax.set_xlim(0, N * 1.1); ax.grid(axis="y", visible=False)
ax.set_xlabel("Number of objects")
title(ax, "From 800 matches to 12 high-priority candidates", "Cumulative checks (in this order); tiers also weigh further checks")
save(fig, "07_selection_funnel.png")

# 8. Catalogue overlap
fig, ax = plt.subplots(figsize=(7, 4))
vc = s.n_surveys.value_counts().sort_index()
ab = s[s.tier.isin(["A", "B"])].n_surveys.value_counts().reindex(vc.index, fill_value=0)
ax.bar(vc.index - 0.18, vc.values, width=0.36, color=GREY, edgecolor=SURF, linewidth=2, label="All matches")
ax.bar(ab.index + 0.18, ab.values, width=0.36, color=SLOT[1], edgecolor=SURF, linewidth=2, label="Tier A + B")
ax.set_xlabel("Number of EB catalogues listing the object"); ax.set_ylabel("Number of objects"); ax.grid(axis="x", visible=False)
title(ax, "Most matches are listed by two or more EB catalogues",
      f"{(s.n_surveys >= 2).mean():.0%} of all matches, {(s[s.tier.isin(['A', 'B'])].n_surveys >= 2).mean():.0%} of tier A + B; VSX ingests other surveys, so not all are independent")
ax.legend(fontsize=9)
save(fig, "08_catalogue_overlap.png")

# summary numbers for the README
summ = {
    "n": N, "frac_plane": float((np.abs(s.b) < 10).mean()), "class_counts": s.cls.value_counts().to_dict(),
    "tierA_classes": s[s.tier == "A"].cls.value_counts().to_dict(), "funnel": steps,
    "lpv_median_P": float(p[lpv].median()), "eb_median_P": float(p[s.cls.str.startswith("Eclipsing")].median()),
    "frac_fir_debris": float(((s.f_ir > 0) & (s.f_ir < 1e-2)).mean()),
    "median_bbteff_A": float(s[s.tier == "A"].bb_teff.median()), "median_bbteff_all": float(s.bb_teff.median()),
    "yso": int(s.is_yso.sum()),
    "tierA_debris": int(((s.f_ir > 0) & (s.f_ir < 1e-2) & (s.tier == "A")).sum()),
    "tierA_P_median": float(p[s.tier == "A"].median()), "all_P_median": float(p.median()),
    "AB_multi": float((s[s.tier.isin(["A", "B"])].n_surveys >= 2).mean()),
}
for k, v in summ.items():
    print(k, v)
