#!/usr/bin/env python3
"""Evaluate J1122 PHOEBE model with separate per-sector light-curve datasets."""

from __future__ import annotations

import json
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import phoebe

import run_j1122_detached_search as base
import run_j1122_morphology_continuation as morph


OUTPUT_DIR = base.REPO_DIR / "outputs" / "j1122_sector_dataset_model"
RAW_LC = (
    base.REPO_DIR
    / "outputs"
    / "j1122_tess_multisector_s10_s11_s37"
    / f"{base.TARGET}_LightCurve_multisector.txt"
)


def load_raw():
    data = np.genfromtxt(
        RAW_LC,
        comments="#",
        names=["time", "flux", "err", "sector"],
    )
    return data


def binned_sector_rows(data, sector, bins=120, min_points=3):
    mask = data["sector"].astype(int) == sector
    time = data["time"][mask]
    flux = data["flux"][mask]
    err = data["err"][mask]
    phase = base.phase_values(time)
    edges = np.linspace(-0.2, 0.8, bins + 1)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = (phase >= lo) & (phase < hi)
        if np.sum(in_bin) < min_points:
            continue
        med = float(np.nanmedian(flux[in_bin]))
        scatter = 1.4826 * float(np.nanmedian(np.abs(flux[in_bin] - med)))
        if not np.isfinite(scatter) or scatter <= 0:
            scatter = float(np.nanmedian(err[in_bin]))
        sigma = max(scatter / np.sqrt(np.sum(in_bin)), 0.0004)
        bin_phase = 0.5 * (lo + hi)
        rows.append(
            (
                base.T0 + bin_phase * base.PERIOD,
                med,
                sigma,
                bin_phase,
                np.sum(in_bin),
            )
        )
    return np.asarray(rows, dtype=float)


def build_sector_bundle(sector_rows):
    b = phoebe.default_binary()
    for sector, rows in sector_rows.items():
        b.add_dataset(
            "lc",
            times=rows[:, 0],
            fluxes=rows[:, 1],
            sigmas=rows[:, 2],
            dataset=f"lc_s{sector}",
            passband="TESS:T",
        )
    b.set_value_all("atm", "blackbody")
    b.set_value_all("ld_mode_bol", "manual")
    b.set_value_all("ld_func_bol", "linear")
    b.set_value_all("ld_coeffs_bol", [0.5])
    b.set_value_all("ld_mode", "manual")
    b.set_value_all("ld_func", "linear")
    b.set_value_all("ld_coeffs", [0.5])
    b.set_value("period@binary", base.PERIOD)
    b.set_value("t0_supconj@binary", base.T0)
    b.set_value("ecc@binary", 0.0)
    b.set_value("q@binary", 1.0)
    b.set_value("sma@binary", 20.0)
    b.set_value("teff@primary", 6000.0)
    b.set_value("irrad_frac_refl_bol@primary", 0.5)
    b.set_value("irrad_frac_refl_bol@secondary", 0.5)
    for dataset in sector_rows:
        b.set_value(f"pblum_mode@lc_s{dataset}", "dataset-scaled")
    b.set_value("distortion_method@primary@phoebe01", "sphere")
    b.set_value("distortion_method@secondary@phoebe01", "sphere")
    b.set_value_all("ntriangles", 1800)
    return b


def set_candidate(b, candidate):
    base.set_candidate(b, candidate)
    for dataset in b.datasets:
        if str(dataset).startswith("lc_s"):
            b.set_value(f"ld_coeffs@primary@{dataset}", [candidate["ld_primary"]])
            b.set_value(f"ld_coeffs@secondary@{dataset}", [candidate["ld_secondary"]])
    b.set_value("ecc@binary", candidate["ecc"])
    b.set_value("per0@binary", candidate["per0"])


def collect_residuals(b, model):
    rows = []
    for dataset in b.datasets:
        if not str(dataset).startswith("lc_s"):
            continue
        times = np.asarray(b.get_value(f"times@{dataset}@lc@dataset"), dtype=float)
        observed = np.asarray(b.get_value(f"fluxes@{dataset}@lc@dataset"), dtype=float)
        predicted = np.asarray(b.get_value(f"fluxes@{dataset}@phoebe01@{model}@lc@model"), dtype=float)
        phase = base.phase_values(times)
        residual = observed - predicted
        sector = int(str(dataset).replace("lc_s", ""))
        rows.append(np.column_stack((phase, observed, predicted, residual, np.full_like(phase, sector))))
    return np.vstack(rows)


def metrics_from_rows(rows, chi2):
    phase, observed, predicted, residual, _ = rows.T
    metrics = base.residual_metrics(base.T0 + phase * base.PERIOD, observed, predicted, base.T0)
    return {
        "chi2": float(chi2),
        **metrics,
        "objective": float(metrics["rms"] + 0.7 * metrics["trend_score"]),
    }


def plot_rows(rows, result, output_path):
    phase, observed, predicted, residual, sector = rows.T
    order = np.argsort(phase)
    fig, (ax, rx) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    for s in sorted(set(sector.astype(int))):
        mask = sector.astype(int) == s
        ax.scatter(phase[mask], observed[mask], s=14, alpha=0.55, linewidths=0, label=f"S{s}")
        rx.scatter(phase[mask], residual[mask], s=12, alpha=0.45, linewidths=0)
    ax.plot(phase[order], predicted[order], color="#dc2626", lw=1.4, label="PHOEBE")
    ax.legend(loc="best", fontsize=9)
    ax.set_ylabel("Flux")
    ax.set_title(
        f"{base.TARGET} sector datasets: chi2={result['chi2']:.1f}, "
        f"rms={result['rms']:.6f}, trend={result['trend_score']:.6f}"
    )
    rx.axhline(0, color="0.45", lw=1)
    rx.set_ylabel("Obs-model")
    rx.set_xlabel("Phase")
    rx.set_xlim(-0.2, 0.8)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def evaluate_candidate(b, candidate, label):
    set_candidate(b, candidate)
    b.run_compute(compute="phoebe01", model=label, overwrite=True, progressbar=False)
    rows = collect_residuals(b, label)
    chi2 = 0.0
    for dataset in b.datasets:
        if str(dataset).startswith("lc_s"):
            chi2 += float(b.calculate_chi2(model=label, dataset=dataset))
    result = {
        "label": label,
        "candidate": dict(candidate),
        **metrics_from_rows(rows, chi2),
    }
    plot_rows(rows, result, OUTPUT_DIR / f"{label}_diagnostic.png")
    b.save(str(OUTPUT_DIR / f"{label}.phoebe"))
    return result


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    raw = load_raw()
    sectors = sorted(set(raw["sector"].astype(int)))
    sector_rows = {sector: binned_sector_rows(raw, sector) for sector in sectors}
    for sector, rows in sector_rows.items():
        np.savetxt(
            OUTPUT_DIR / f"sector_{sector}_binned.txt",
            rows,
            header="Time Flux Flux_Err Phase N",
        )
        print(f"sector {sector}: {len(rows)} bins", flush=True)

    b = build_sector_bundle(sector_rows)
    targeted = json.loads(
        (base.REPO_DIR / "outputs" / "j1122_targeted_eclipse_grid" / "j1122_targeted_eclipse_grid_result.json").read_text()
    )
    candidates = [
        ("sector_start", morph.load_start()),
        ("sector_best_chi2", targeted["best_final_chi2"]["candidate"]),
        ("sector_best_objective", targeted["best_final_objective"]["candidate"]),
    ]
    results = []
    for label, cand in candidates:
        cand = dict(cand)
        cand.setdefault("ld_primary", 0.5)
        cand.setdefault("ld_secondary", 0.5)
        cand.setdefault("ecc", 0.0)
        cand.setdefault("per0", 0.0)
        result = evaluate_candidate(b, cand, label)
        results.append(result)
        print(
            f"{label}: chi2={result['chi2']:.2f} rms={result['rms']:.6f} "
            f"primary_abs={result['primary_mean_abs']:.6f}",
            flush=True,
        )

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": base.TARGET,
        "raw_lc": str(RAW_LC.relative_to(base.REPO_DIR)),
        "sector_bins": {str(k): int(len(v)) for k, v in sector_rows.items()},
        "results": results,
        "best_chi2": min(results, key=lambda r: r["chi2"]),
        "best_objective": min(results, key=lambda r: r["objective"]),
        "status": "success",
    }
    with (OUTPUT_DIR / "j1122_sector_dataset_model_result.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
