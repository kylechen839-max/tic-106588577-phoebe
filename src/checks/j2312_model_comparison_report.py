"""J2312 at 22 um: figure report comparing the getsf-style (W4-only) and 12 um-map (method B) nebula models.

Reuses j2312_model_figure.py (method B) and j2312_getsf_figure.py (getsf-style, W = 120") and writes four figures
to outputs/j2312_contamination/report/:
  fig1_images.png  data, nebula models, data - nebula, residuals, on common color scales
  fig2_profiles.png  east-west cuts across the ridge and radial profiles after nebula removal
  fig3_flux.png  W4 at the star for every background model, with sigma
  fig4_split.png  point + extended split of (data - method B nebula): fit quality vs extended size
Usage: .venv-phoebe/bin/python src/checks/j2312_model_comparison_report.py
"""
import runpy
from pathlib import Path
import numpy as np
from scipy import ndimage
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
b = runpy.run_path(str(here / "j2312_model_figure.py"))
g = runpy.run_path(str(here / "j2312_getsf_figure.py"))
m = b["m"]; PIX = m.PIX
R, X, Y, k = b["R"], b["X"], b["Y"], b["k"]
OUT = m.D / "report"; OUT.mkdir(exist_ok=True)
F_PHOT, E_CAT = m.F_PHOT, m.E_CAT
ext = [-160, 160, -160, 160]
plt.rcParams.update({"font.size": 9})

z0 = np.nanmedian(b["s4"][(R > 140) & (R < 160)])
data = (b["s4"] - z0) * k
w3 = b["s3"] - b["coef3"][0] * b["P3"]
w3 = (w3 - np.nanmedian(w3[(R > 140) & (R < 160)])) * 1.998e-6 * 1e3
nebB = (b["neb"] - z0) * k
nebG = (g["neb"] - np.nanmedian(g["neb"][(g["R"] > 140) & (g["R"] < 160)])) * g["k"]
subB = (b["s4"] - b["neb"]) * k                      # data - nebula = star + whatever is left
subG = (g["s4"] - g["neb"]) * g["k"]
resB = np.where(R < 150, b["res"] * k, np.nan)
resG = np.where(g["R"] < 150, g["res"] * g["k"], np.nan)
F_B, F_G = b["F_B"], g["F"]

# ---------- fig 1: images ----------
v = np.nanpercentile(data[R < 150], [2, 99.7])
vz = np.nanpercentile(subB[R < 60], [1, 99.8])
lr = np.nanpercentile(np.abs(resG[np.isfinite(resG)]), 99)
fig, axs = plt.subplots(2, 4, figsize=(17, 8.6))
rows = [("getsf-style (22 µm only)", nebG, subG, resG, F_G),
        ("12 µm map (method B)", nebB, subB, resB, F_B)]
left = [(data, "22 µm data (W4)", "magma", v), (w3, "12 µm data, star removed (W3)\n= the template source", "magma",
                                                     np.nanpercentile(w3[R < 150], [2, 99.5]))]
for i, (name, neb, sub, res, F) in enumerate(rows):
    a, t, cm, lim = left[i]
    panels = [(a, t, cm, lim, 120),
              (neb, f"{name}\nnebula model", "magma", v, 120),
              (sub, f"data − nebula model (zoom)\nstar fit: {F:.1f} mJy = {F / F_PHOT:.2f}× stars alone", "magma", vz, 60),
              (res, "residual: data − nebula − star\n(red = light left unmodeled)", "coolwarm", (-lr, lr), 120)]
    for j, (img, title, cm, lim, zz) in enumerate(panels):
        ax = axs[i, j]
        im = ax.imshow(img, origin="lower", cmap=cm, vmin=lim[0], vmax=lim[1], extent=ext)
        ax.plot(0, 0, "+", color="c", ms=9, mew=1.4)
        ax.set_xlim(-zz, zz); ax.set_ylim(-zz, zz); ax.set_title(title, fontsize=9.5)
        if j == 2:
            ax.add_patch(plt.Circle((0, 0), 8.5, fill=False, ls="--", color="w", lw=0.8))
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=7)
    axs[i, 1].set_ylabel(name, fontsize=11, weight="bold")
for ax in axs.flat:
    ax.tick_params(labelsize=7)
fig.suptitle("J2312 at 22 µm: the two nebula models side by side (arcsec; color = mJy/pixel; dashed circle = point-source FWHM)",
             fontsize=12)
fig.tight_layout(); fig.savefig(OUT / "fig1_images.png", dpi=130); plt.close(fig)

# ---------- fig 2: cuts and radial profiles ----------
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
xs = np.arange(-90, 91, 4.0)
w4_nostar = (b["s4"] - b["star"] - z0) * k
def cut(a, y0, y1):
    return np.array([np.nanmean(a[(Y >= y0) & (Y < y1) & (X >= x - 2) & (X < x + 2)]) for x in xs])
for ax, (y0, y1), lab in [(axs[0], (15, 35), "15–35″ north"), (axs[1], (-35, -15), "15–35″ south")]:
    ax.plot(xs, cut(w4_nostar, y0, y1), "-", color="k", lw=2, label="22 µm data (star removed)")
    ax.plot(xs, cut(nebB, y0, y1), "-", color="C0", lw=1.6, label="12 µm-map model")
    ax.plot(xs, cut(nebG, y0, y1), "--", color="C1", lw=1.6, label="getsf-style model")
    ax.axvline(0, color="c", lw=0.8); ax.axhline(0, color="0.7", lw=0.6)
    ax.set_title(f"East–west cut across the ridge, {lab} of the star"); ax.set_xlabel("east–west offset (arcsec)")
    ax.set_ylabel("mJy / pixel"); ax.legend(fontsize=8)
ax = axs[2]
edges = np.arange(0, 47, 2.75); rc = 0.5 * (edges[1:] + edges[:-1])
prof = lambda a, RR: np.array([np.nanmean(a[(RR >= e0) & (RR < e1)]) for e0, e1 in zip(edges[:-1], edges[1:])])
pB, pG = prof(subB, R), prof(subG, g["R"])
psf = prof(b["P4"] * k, R)
ax.plot(rc, pB / pB[0], "o-", color="C0", label="data − 12 µm-map nebula")
ax.plot(rc, pG / pG[0], "s--", color="C1", label="data − getsf-style nebula")
ax.plot(rc, psf / psf[0], "-", color="k", lw=2, label="point source (field-star PSF)")
ax.axhline(0, color="0.7", lw=0.6)
fwB = g["t"].gfit(b["starplus"], b["hp"] + b["dx"], b["hp"] + b["dy"], 20)[3]
fwG = g["fw_sub"]; fwP = g["fw_psf"]
ax.set_title(f"Radial profile, peak-normalized. FWHM within 20″:\n12 µm map {fwB:.1f}″, getsf-style {fwG:.1f}″, point source {fwP:.1f}″",
             fontsize=9)
ax.set_xlabel("distance from the star (arcsec)"); ax.set_ylabel("relative brightness"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "fig2_profiles.png", dpi=130); plt.close(fig)

# ---------- fig 3: flux by model ----------
models = [("catalogue (flat sky)", 15.36, "0.6"), ("tilted plane", 14.55, "0.6"), ("quadratic", 12.30, "C2"),
          ("size separation", 12.40, "C2"), ("getsf-style", F_G, "C1"), ("12 µm map (method B)", F_B, "C0")]
fig, ax = plt.subplots(figsize=(8.5, 4.2))
y = np.arange(len(models))
ax.barh(y, [v for _, v, _ in models], xerr=E_CAT, color=[c for *_, c in models], capsize=3, alpha=0.9)
ax.axvline(F_PHOT, color="k", ls="--", lw=1.5)
ax.text(F_PHOT + 0.15, -0.6, f"stars alone {F_PHOT:.2f} mJy", fontsize=8)
for yi, (_, vv, _) in zip(y, models):
    ax.text(vv + 1.0, yi, f"{vv:.2f} mJy · {vv / F_PHOT:.2f}× · {(vv - F_PHOT) / E_CAT:.1f}σ", va="center", fontsize=8.5)
ax.set_yticks(y); ax.set_yticklabels([n for n, *_ in models]); ax.invert_yaxis()
ax.set_xlim(0, 24); ax.set_xlabel("W4 flux at the star (mJy); error bar = catalogue error 0.88 mJy")
ax.set_title("W4 at J2312 by background model (gray: nebula ignored; color: nebula modeled)")
fig.tight_layout(); fig.savefig(OUT / "fig3_flux.png", dpi=130); plt.close(fig)

# ---------- fig 4: point + extended split of (data - method B nebula) ----------
s4, neb, P4, mjy = b["s4"], b["neb"], b["P4"], b["k"]
yv = s4 - neb
M = R < 40
one = np.ones_like(s4)
noise = 1.4826 * np.nanmedian(np.abs(b["res"][(R > 30) & (R < 150)] - np.nanmedian(b["res"][(R > 30) & (R < 150)])))
def fit(cols, yy):
    A = np.c_[[c_[M] for c_ in cols]].T
    sol, *_ = np.linalg.lstsq(A, yy[M], rcond=None)
    return sol, float(np.sum((yy[M] - A @ sol) ** 2) / noise ** 2)
_, chi_pt = fit([P4, one], yv)
sizes = np.array([4, 6, 8, 10, 13, 16, 20, 25, 30])
fixed, free = [], []
for fe in sizes:
    E = ndimage.gaussian_filter(P4, fe / 2.3548 / PIX)
    s2, c2 = fit([E, one], yv - F_PHOT / mjy * P4); fixed.append((s2[0] * mjy, c2))
    s1, c1 = fit([P4, E, one], yv); free.append((s1[0] * mjy, s1[1] * mjy, c1))
fixed, free = np.array(fixed), np.array(free)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axs[0]
ok = free[:, 0] >= 0
ax.plot(sizes, fixed[:, 1] - chi_pt, "o-", color="C3", label="point fixed at stars alone + extended")
ax.plot(sizes[ok], free[ok, 2] - chi_pt, "s--", color="C0", label="point free (≥ 0) + extended")
ax.axhline(0, color="k", lw=1, label="point source only")
ax.set_xlabel("extended component size (FWHM, arcsec)"); ax.set_ylabel("Δχ² vs point source only (lower = better)")
ax.set_title("Is an extended component needed?\n(pixels treated as independent: Δχ² overstated)")
ax.legend(fontsize=8)
ax = axs[1]
labels = ["point only", "point = stars alone\n+ 8″ extended", "free point\n+ 13″ extended", "free point\n+ 20″ extended"]
i8, i13, i20 = (list(sizes).index(s) for s in (8, 13, 20))
pt = [fit([P4, one], yv)[0][0] * mjy, F_PHOT, free[i13, 0], free[i20, 0]]
ex = [0, fixed[i8, 0], free[i13, 1], free[i20, 1]]
xb = np.arange(4)
ax.bar(xb, [min(p, F_PHOT) for p in pt], color="0.55", label="stars alone (photosphere)")
ax.bar(xb, [max(p - F_PHOT, 0) for p in pt], bottom=[min(p, F_PHOT) for p in pt], color="C1", label="compact excess (disk-like)")
ax.bar(xb, ex, bottom=pt, color="C0", label="extended dust")
for xi, p, e in zip(xb, pt, ex):
    ax.text(xi, p + e + 0.25, f"{p + e:.1f}", ha="center", fontsize=8)
ax.set_xticks(xb); ax.set_xticklabels(labels, fontsize=8); ax.set_ylabel("W4 flux within 40″ (mJy)")
ax.set_title("Where the 22 µm flux goes (after the 12 µm-map nebula is removed)")
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(OUT / "fig4_split.png", dpi=130); plt.close(fig)
print("chi2 point only", round(chi_pt), "fixed", np.round(fixed, 2).tolist())
print(OUT)
