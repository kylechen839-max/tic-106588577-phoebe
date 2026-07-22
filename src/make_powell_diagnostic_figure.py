import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import phoebe


bundle_path = "outputs/tic_106588577_final.phoebe"
out_path = "outputs/powell_diagnostic_lightcurve.png"

b = phoebe.load(bundle_path)

period = float(b.get_value("period@binary@orbit@component"))
t0 = float(b.get_value("t0_supconj@binary@orbit@component"))

obs_t = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
obs_f = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
model_f = np.asarray(
    b.get_value("fluxes@lc01@phoebe01@after_powell@lc@model"),
    dtype=float,
)

phase = ((obs_t - t0) / period) % 1.0
phase = np.where(phase > 0.8, phase - 1.0, phase)
order = np.argsort(phase)
resid = obs_f - model_f

fig, (ax, rx) = plt.subplots(
    2,
    1,
    figsize=(9, 7),
    sharex=True,
    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
)

ax.scatter(
    phase,
    obs_f,
    s=9,
    color="0.25",
    alpha=0.35,
    label="Observed lc01",
    linewidths=0,
)
ax.plot(phase[order], model_f[order], color="#2563eb", lw=2.4, label="After Powell")

rx.axhline(0, color="0.35", lw=1)
rx.scatter(phase, resid, s=8, color="0.25", alpha=0.35, linewidths=0)
rx.plot(phase[order], resid[order], color="#2563eb", lw=1.0, alpha=0.85)

ax.set_ylabel("Flux")
rx.set_ylabel("Obs - model")
rx.set_xlabel("Phase relative to primary eclipse")
ax.set_xlim(-0.30, 0.80)
ax.legend(loc="best", frameon=False)
ax.set_title("TIC 106588577: Powell model and residuals")
fig.tight_layout()
fig.savefig(out_path, dpi=220, bbox_inches="tight")

print(out_path)
print("period", period)
print("t0", t0)
print("chi2", b.calculate_chi2(model="after_powell", dataset="lc01"))
print("resid rms", float(np.sqrt(np.mean(resid**2))))
