import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import phoebe


bundle_path = "outputs/tic_106588577_final.phoebe"
out_path = "outputs/emcee_diagnostic_lightcurve.png"

b = phoebe.load(bundle_path)

period = float(b.get_value("period@binary@orbit@component"))
t0 = float(b.get_value("t0_supconj@binary@orbit@component"))

obs_t = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
obs_f = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)

phase = ((obs_t - t0) / period) % 1.0
phase = np.where(phase > 0.7, phase - 1.0, phase)
order = np.argsort(phase)

powell_f = np.asarray(
    b.get_value("fluxes@lc01@phoebe01@after_powell@lc@model"),
    dtype=float,
)
emcee_f = np.asarray(
    b.get_value("fluxes@lc01@ellcbackend@emcee_posterior_models@lc@model"),
    dtype=float,
)
if emcee_f.ndim == 1:
    emcee_f = emcee_f[None, :]

median_f = np.nanmedian(emcee_f, axis=0)
resid = obs_f - median_f

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
    color="0.2",
    alpha=0.35,
    label="Observed lc01",
    linewidths=0,
)
ax.plot(phase[order], powell_f[order], color="#2563eb", lw=2.0, label="After Powell")

for i, f in enumerate(emcee_f):
    ax.plot(
        phase[order],
        f[order],
        color="#f97316",
        lw=1.2,
        alpha=0.35,
        label="emcee samples" if i == 0 else None,
    )
ax.plot(phase[order], median_f[order], color="#dc2626", lw=2.2, label="emcee median")

rx.axhline(0, color="0.35", lw=1)
rx.scatter(phase, resid, s=8, color="0.2", alpha=0.35, linewidths=0)
rx.plot(phase[order], resid[order], color="#dc2626", lw=1.0, alpha=0.85)

ax.set_ylabel("Flux")
rx.set_ylabel("Obs - model")
rx.set_xlabel("Phase relative to primary eclipse")
ax.set_xlim(-0.30, 0.80)
ax.legend(loc="best", frameon=False)
ax.set_title("TIC 106588577: Powell and emcee posterior model comparison")
fig.tight_layout()
fig.savefig(out_path, dpi=220, bbox_inches="tight")

print(out_path)
print("posterior curves", emcee_f.shape[0])
print("period", period)
print("t0", t0)
print("resid rms", float(np.sqrt(np.mean(resid**2))))
