"""J2312: is the compact W3/W4 excess the same material as the surrounding 22 um halo?

Compares the 12/22 um flux ratio (F_nu) of
  (a) the diffuse halo around the star (annuli 25-150", above the 150-200" sky, and a pixel-level
      W3-vs-W4 regression), with
  (b) the compact excess at the star (16.5" aperture minus the binary photosphere), for the
      catalogue 50-70" sky and for the halo extrapolated under the star (outputs/j2312_contamination/w4_bgsub.py).
The images (SkyView AllWISE atlas, 10', 2"/px, DN) are calibrated with the star's own AllWISE
magnitude in the same aperture. Also fits a Gaussian width to the W3 and W4 source.
"""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "outputs" / "j2312_contamination"
N = "J231201.40+532028.6"
C = SkyCoord(348.0058494, 53.3413022, unit="deg")
PIX, RAP = 2.0, 16.5
ZP = {"W3": 31.674, "W4": 8.363}
MAG = {"W3": 8.062, "W4": 6.840}
# (label, W3 mag, W3 F/F*, W4 mag, W4 F/F*) from w4_bgsub_log.txt (free-Teff photosphere)
EXCESS_CASES = [("catalogue sky (50-70\")", 8.062, 1.257643, 6.840, 3.629307),
                ("halo extrapolated under star", 8.245, 1.06, 7.142, 2.75)]


def bb_ratio(T, l3=11.56, l4=22.09):
    f = lambda lam: (1 / lam ** 3) / np.expm1(14387.77 / (lam * T))
    return f(l3) / f(l4)


def load(band, fname):
    h = fits.open(D / f"{N}_{fname}_10arcmin.fits")[0]
    img = h.data.astype(float)
    x0, y0 = WCS(h.header).world_to_pixel(C)
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - x0, yy - y0) * PIX
    far = np.nanmedian(img[(r >= 150) & (r < 200)])
    near = np.nanmedian(img[(r >= 50) & (r < 70)])
    k = ZP[band] * 10 ** (-0.4 * MAG[band]) / np.nansum(img[r < RAP] - near)  # Jy per summed DN
    return dict(img=img, r=r, xx=xx, yy=yy, x0=x0, y0=y0, far=far, k=k)


def gauss_fwhm(d):
    m = d["r"] < 40
    xx, yy, img = d["xx"][m], d["yy"][m], d["img"][m]

    def resid(p):
        A, xc, yc, s, c0, cx, cy = p
        return A * np.exp(-((xx - xc) ** 2 + (yy - yc) ** 2) / (2 * s ** 2)) + c0 + cx * (xx - d["x0"]) + cy * (yy - d["y0"]) - img
    p0 = [img.max() - np.median(img), d["x0"], d["y0"], 4, np.median(img), 0, 0]
    p = least_squares(resid, p0).x
    return 2.3548 * p[3] * PIX, np.hypot(p[1] - d["x0"], p[2] - d["y0"]) * PIX


def main():
    w3, w4 = load("W3", "WISE12"), load("W4", "WISE22")
    lines = []
    rings = [(25, 40), (40, 55), (55, 70), (70, 100), (100, 150)]
    halo = []
    for a, b in rings:
        s3 = (np.nanmedian(w3["img"][(w3["r"] >= a) & (w3["r"] < b)]) - w3["far"]) * w3["k"]
        s4 = (np.nanmedian(w4["img"][(w4["r"] >= a) & (w4["r"] < b)]) - w4["far"]) * w4["k"]
        halo.append(s3 / s4)
        lines.append(f"halo {a}-{b}\": F12/F22 = {s3 / s4:.2f}")
    m = (w3["r"] >= 25) & (w3["r"] < 150)
    a3 = (w3["img"][m] - w3["far"]) * w3["k"]
    a4 = (w4["img"][m] - w4["far"]) * w4["k"]
    q = np.quantile(a4, np.linspace(0, 1, 11))
    xb = [np.median(a4[(a4 >= q[i]) & (a4 < q[i + 1])]) for i in range(10)]
    yb = [np.median(a3[(a4 >= q[i]) & (a4 < q[i + 1])]) for i in range(10)]
    slope = np.polyfit(xb, yb, 1)[0]
    lines.append(f"halo pixel regression (25-150\", 10 W4 bins): dF12/dF22 = {slope:.2f}")
    exc = []
    for lab, m3, R3, m4, R4 in EXCESS_CASES:
        F3, F4 = ZP["W3"] * 1e3 * 10 ** (-0.4 * m3), ZP["W4"] * 1e3 * 10 ** (-0.4 * m4)
        e3, e4 = F3 * (1 - 1 / R3), F4 * (1 - 1 / R4)
        s3 = F3 * 0.4 * np.log(10) * 0.03  # pipeline 0.03 mag error floor (W3 SNR 51)
        s4 = F4 * 0.4 * np.log(10) * 0.062
        rat = e3 / e4
        srat = rat * np.hypot(s3 / e3, s4 / e4)
        Ts = np.linspace(80, 400, 3201)
        T = Ts[np.argmin(np.abs(bb_ratio(Ts) - rat))]
        exc.append((lab, rat, srat, T))
        lines.append(f"compact excess, {lab}: W3 {e3:.2f}±{s3:.2f} mJy, W4 {e4:.2f}±{s4:.2f} mJy, F12/F22 = {rat:.2f}±{srat:.2f} (blackbody {T:.0f} K)")
    for band, d in (("W3", w3), ("W4", w4)):
        fw, off = gauss_fwhm(d)
        lines.append(f"{band} source: Gaussian FWHM {fw:.1f}\", centroid offset {off:.1f}\"")
    (D / "j2312_excess_colour_log.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    fig, ax = plt.subplots(figsize=(6.4, 4))
    Ts = np.linspace(60, 600, 400)
    ax.plot(Ts, bb_ratio(Ts), "k-", lw=1, label="single-temperature blackbody")
    ax.axhspan(min(halo + [slope]), max(halo + [slope]), color="C0", alpha=0.25, label="surrounding 22 µm halo (25–150″)")
    for (lab, rat, srat, T), col in zip(exc, ("C3", "C1")):
        ax.errorbar(T, rat, yerr=srat, fmt="o", color=col, capsize=3, label=f"compact excess, {lab}")
    ax.set_yscale("log"); ax.set_xlim(60, 600); ax.set_ylim(0.01, 3)
    ax.set_xlabel("dust temperature (K)"); ax.set_ylabel("F$_ν$(12 µm) / F$_ν$(22 µm)")
    ax.set_title("J2312: colour of the excess vs. the nebula around it", fontsize=10)
    ax.legend(fontsize=7, loc="lower right"); ax.grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(D / "j2312_excess_colour.png", dpi=150)


if __name__ == "__main__":
    main()
