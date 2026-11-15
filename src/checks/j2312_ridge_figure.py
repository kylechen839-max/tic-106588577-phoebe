"""J2312: where the nebular ridge is. Annotated W3/W4 images with the star's point source removed, the residual
the getsf-style model leaves, and east-west brightness cuts across the ridge above and below the star.
Output: outputs/j2312_contamination/j2312_ridge_figure.png
"""
import runpy
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
b = runpy.run_path(str(here / "j2312_model_figure.py"))     # method B: 12 um template
g = runpy.run_path(str(here / "j2312_getsf_figure.py"))     # getsf-style
m = b["m"]; PIX = m.PIX
k = b["k"]; R = b["R"]; X = b["X"]; Y = b["Y"]
z0 = np.nanmedian(b["s4"][(R > 140) & (R < 160)])
w4_nostar = (b["s4"] - b["star"] - z0) * k                    # W4 with the fitted point source removed
s3 = b["s3"]; P3 = b["P3"]; c3 = b["coef3"]
w3_nostar = s3 - c3[0] * P3
w3_nostar = (w3_nostar - np.nanmedian(w3_nostar[(R > 140) & (R < 160)])) * 1.998e-6 * 1e3   # mJy/pixel
gres = np.where(g["R"] < 150, g["res"] * g["k"], np.nan)
neb_b = (b["neb"] - z0) * k
neb_g = (g["neb"] - np.nanmedian(g["neb"][(g["R"] > 140) & (g["R"] < 160)])) * g["k"]

ext = [-160, 160, -160, 160]
fig, axs = plt.subplots(2, 3, figsize=(15.5, 10))


def show(ax, a, title, cmap="magma", lim=None):
    if lim is None:
        lim = np.nanpercentile(a[R < 150], [2, 99.5])
    im = ax.imshow(a, origin="lower", cmap=cmap, vmin=lim[0], vmax=lim[1], extent=ext)
    ax.plot(0, 0, "+", color="c", ms=10, mew=1.5)
    ax.set_title(title, fontsize=10); ax.set_xlim(-120, 120); ax.set_ylim(-150, 120)
    ax.set_xticks([-100, 0, 100]); ax.set_yticks([-100, 0, 100]); ax.tick_params(labelsize=8)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=7)
    return lim


def ridge_marks(ax, col="w"):
    ax.annotate("", xy=(-4, 45), xytext=(-55, 85), arrowprops=dict(arrowstyle="->", color=col, lw=1.5))
    ax.annotate("", xy=(-8, -60), xytext=(-75, -40), arrowprops=dict(arrowstyle="->", color=col, lw=1.5))
    ax.text(-115, 92, "ridge", color=col, fontsize=11, weight="bold")
    ax.text(-118, -36, "ridge", color=col, fontsize=11, weight="bold")
    ax.plot([-30, 30, 30, -30, -30], [-35, -35, -15, -15, -35], color="c", lw=0.8, ls="--")
    ax.plot([-30, 30, 30, -30, -30], [15, 15, 35, 35, 15], color="c", lw=0.8, ls="--")


lim4 = show(axs[0, 0], w4_nostar, "22 µm (W4) with the star's point source removed\n(this is the nebula around J2312)")
ridge_marks(axs[0, 0])
show(axs[0, 1], w3_nostar, "12 µm (W3) with the star removed\nsame ridge, sharper")
ridge_marks(axs[0, 1])
show(axs[0, 2], neb_b, "Nebula model from the 12 µm map\n(captures the ridge)", lim=lim4)
show(axs[1, 0], neb_g, "getsf-style nebula model\n(too smooth: misses the ridge)", lim=lim4)
lr = np.nanpercentile(np.abs(gres), 99)
show(axs[1, 1], gres, "getsf-style residual\n(red = light the model left in: the ridge)", cmap="coolwarm", lim=(-lr, lr))
ridge_marks(axs[1, 1], col="k")
for ax in axs.flat[:5]:
    ax.set_xlabel("arcsec (east ← → west, as on the sky grid)", fontsize=7)

# east-west cuts in the two dashed boxes (avoiding the star itself)
ax = axs[1, 2]
xs = np.arange(-90, 91, 4.0)
def cut(a, y0, y1):
    return np.array([np.nanmean(a[(Y >= y0) & (Y < y1) & (X >= x - 2) & (X < x + 2)]) for x in xs])
for (y0, y1), ls in [((15, 35), "-"), ((-35, -15), "--")]:
    ax.plot(xs, cut(w4_nostar, y0, y1), ls, color="C3", lw=2, label=f"22 µm data, {y0}″ to {y1}″")
    ax.plot(xs, cut(neb_b, y0, y1), ls, color="C0", lw=1.3, label=f"12 µm-map model, {y0}″ to {y1}″")
    ax.plot(xs, cut(neb_g, y0, y1), ls, color="0.5", lw=1.3, label=f"getsf-style model, {y0}″ to {y1}″")
ax.axvline(0, color="c", lw=0.8)
ax.set_xlabel("east–west offset from the star (arcsec)"); ax.set_ylabel("brightness (mJy / pixel)")
ax.set_title("Cuts across the ridge, 15–35″ north and south of the star\n(star excluded)", fontsize=10)
ax.legend(fontsize=7)
fig.suptitle("J2312: the north–south nebular ridge running through the star", fontsize=13)
fig.tight_layout()
out = m.D / "j2312_ridge_figure.png"
fig.savefig(out, dpi=140)
print(out)
