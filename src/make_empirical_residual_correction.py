import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import phoebe
from scipy.interpolate import UnivariateSpline


REPO_DIR = Path(__file__).resolve().parents[1]
BUNDLE_PATH = REPO_DIR / "results" / "final" / "tic_106588577_multisector_controlled_best.phoebe"
OUTPUT_DIR = REPO_DIR / "outputs" / "empirical_residual_correction"
MODEL = "controlled_search_best"
DATASET = "lc01"
COMPUTE = "phoebe01"


def phase_values(times, period, t0):
    phase = ((times - t0) / period) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def binned_medians(phase, residuals, n_bins=160):
    bins = np.linspace(-0.2, 0.8, n_bins + 1)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (phase >= lo) & (phase < hi)
        if not np.any(mask):
            continue
        center = 0.5 * (lo + hi)
        median = np.median(residuals[mask])
        scatter = 1.4826 * np.median(np.abs(residuals[mask] - median))
        rows.append((center, median, max(scatter, 1e-5), int(mask.sum())))
    return np.asarray(rows, dtype=float)


def trend_metrics(phase, residuals):
    primary = np.abs(phase) < 0.045
    secondary = np.abs(phase - 0.5) < 0.035
    shoulders = ((np.abs(phase) >= 0.045) & (np.abs(phase) < 0.11)) | (
        (np.abs(phase - 0.5) >= 0.035) & (np.abs(phase - 0.5) < 0.09)
    )

    def mean(mask):
        return float(np.mean(residuals[mask]))

    def rms(mask):
        return float(np.sqrt(np.mean(residuals[mask] ** 2)))

    def mean_abs(mask):
        return float(np.mean(np.abs(residuals[mask])))

    return {
        "rms": float(np.sqrt(np.mean(residuals**2))),
        "mean_abs": float(np.mean(np.abs(residuals))),
        "max_abs": float(np.max(np.abs(residuals))),
        "primary_core_mean": mean(primary),
        "secondary_core_mean": mean(secondary),
        "primary_core_rms": rms(primary),
        "secondary_core_rms": rms(secondary),
        "shoulder_mean_abs": mean_abs(shoulders),
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    b = phoebe.load(str(BUNDLE_PATH))
    period = float(b.get_value("period@binary@orbit@component"))
    t0 = float(b.get_value("t0_supconj@binary@orbit@component"))
    times = np.asarray(b.get_value(f"times@{DATASET}@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value(f"fluxes@{DATASET}@lc@dataset"), dtype=float)

    try:
        model_flux = np.asarray(
            b.get_value(f"fluxes@{DATASET}@{COMPUTE}@{MODEL}@lc@model"),
            dtype=float,
        )
    except Exception:
        b.run_compute(compute=COMPUTE, model=MODEL, overwrite=True)
        model_flux = np.asarray(
            b.get_value(f"fluxes@{DATASET}@{COMPUTE}@{MODEL}@lc@model"),
            dtype=float,
        )

    phase = phase_values(times, period, t0)
    raw_residual = observed - model_flux
    bins = binned_medians(phase, raw_residual)

    # The smoothing value is deliberately nonzero: this is meant to remove the
    # coherent phase-scale trend, not interpolate every noisy data point.
    weights = 1.0 / bins[:, 2]
    spline = UnivariateSpline(bins[:, 0], bins[:, 1], w=weights, s=len(bins) * 0.55)
    correction = spline(phase)
    corrected_model = model_flux + correction
    corrected_residual = observed - corrected_model

    order = np.argsort(phase)
    correction_grid = np.linspace(-0.2, 0.8, 600)
    correction_curve = spline(correction_grid)
    corrected_bins = binned_medians(phase, corrected_residual)

    fig, (ax, rx) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    ax.scatter(phase, observed, s=8, color="0.25", alpha=0.45, linewidths=0)
    ax.plot(phase[order], model_flux[order], color="#dc2626", lw=1.8, alpha=0.55, label="PHOEBE")
    ax.plot(
        phase[order],
        corrected_model[order],
        color="#2563eb",
        lw=2.0,
        label="PHOEBE + smooth residual correction",
    )
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylabel("Flux")
    ax.set_title("Empirical residual-corrected diagnostic on multi-sector binned curve")

    rx.axhline(0, color="0.45", lw=1)
    rx.scatter(phase, corrected_residual, s=8, color="0.25", alpha=0.42, linewidths=0)
    rx.plot(
        corrected_bins[:, 0],
        corrected_bins[:, 1],
        color="#2563eb",
        lw=1.8,
        label="binned corrected residual median",
    )
    rx.legend(loc="upper right", fontsize=8)
    rx.set_xlabel("Phase")
    rx.set_ylabel("Obs-model")
    ax.set_xlim(-0.2, 0.8)
    fig.savefig(OUTPUT_DIR / "empirical_residual_corrected_diagnostic.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    np.savetxt(
        OUTPUT_DIR / "empirical_residual_corrected_curve.txt",
        np.column_stack([times, phase, observed, model_flux, correction, corrected_model, corrected_residual]),
        header="time phase observed_flux phoebe_flux empirical_correction corrected_model_flux corrected_residual",
    )
    np.savetxt(
        OUTPUT_DIR / "empirical_residual_correction_spline.txt",
        np.column_stack([correction_grid, correction_curve]),
        header="phase empirical_correction",
    )

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundle": str(BUNDLE_PATH.relative_to(REPO_DIR)),
        "source_model": MODEL,
        "method": "smooth spline fit to binned PHOEBE residual medians",
        "note": "Empirical diagnostic correction, not a physical PHOEBE binary solution.",
        "raw_phoebe_residuals": trend_metrics(phase, raw_residual),
        "corrected_residuals": trend_metrics(phase, corrected_residual),
        "n_phase_bins": int(len(bins)),
        "spline_smoothing": float(len(bins) * 0.55),
    }
    with open(OUTPUT_DIR / "empirical_residual_correction_result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
