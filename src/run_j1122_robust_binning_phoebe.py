#!/usr/bin/env python3
"""Build robust J1122 phase-bin variants and evaluate PHOEBE eclipse bias."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_j1122_detached_search as base
import run_j1122_morphology_continuation as morph


OUTPUT_DIR = base.REPO_DIR / "outputs" / "j1122_robust_binning_phoebe"
RAW_LC = (
    base.REPO_DIR
    / "outputs"
    / "j1122_tess_multisector_s10_s11_s37"
    / f"{base.TARGET}_LightCurve_multisector.txt"
)


def load_raw():
    return np.genfromtxt(RAW_LC, comments="#", names=["time", "flux", "err", "sector"])


def phase(time):
    return base.phase_values(time)


def robust_scatter(values):
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    scatter = 1.4826 * np.nanmedian(np.abs(values - med))
    if not np.isfinite(scatter) or scatter <= 0:
        scatter = np.nanstd(values)
    return float(scatter) if np.isfinite(scatter) and scatter > 0 else 0.0004


def bin_raw(data, label, sector_filter=None, bins=220, min_points=3, require_sector_count=None):
    phases = phase(data["time"])
    sectors = data["sector"].astype(int)
    use = np.isfinite(data["flux"])
    if sector_filter is not None:
        use &= sector_filter(data, phases, sectors)
    edges = np.linspace(-0.2, 0.8, bins + 1)
    rows = []
    sector_detail = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = use & (phases >= lo) & (phases < hi)
        sector_meds = []
        sector_counts = []
        for sector in sorted(set(sectors[in_bin])):
            smask = in_bin & (sectors == sector)
            if np.sum(smask) >= min_points:
                sector_meds.append(float(np.nanmedian(data["flux"][smask])))
                sector_counts.append((int(sector), int(np.sum(smask))))
        if not sector_meds:
            continue
        if require_sector_count is not None and len(sector_meds) < require_sector_count:
            continue
        bin_phase = 0.5 * (lo + hi)
        flux = float(np.nanmedian(sector_meds))
        inter_sector = robust_scatter(sector_meds)
        sigma = max(inter_sector / max(np.sqrt(len(sector_meds)), 1.0), 0.0004)
        rows.append((base.T0 + bin_phase * base.PERIOD, flux, sigma, bin_phase, sum(n for _, n in sector_counts), len(sector_meds)))
        sector_detail.append({"phase": bin_phase, "sector_counts": sector_counts, "sector_medians": sector_meds})
    rows = np.asarray(rows, dtype=float)
    path = OUTPUT_DIR / f"{label}_lightcurve.txt"
    np.savetxt(path, rows, header="Time Flux Flux_Err Phase N NSector", fmt="%.8e")
    with (OUTPUT_DIR / f"{label}_bin_details.json").open("w", encoding="utf-8") as handle:
        json.dump(sector_detail, handle, indent=2)
        handle.write("\n")
    return path, rows


def sigma_clip_rows(label, rows, clip_sigma=4.0):
    flux = rows[:, 1]
    phase_values = rows[:, 3]
    continuum = (
        ((phase_values > 0.08) & (phase_values < 0.38))
        | ((phase_values > 0.58) & (phase_values < 0.78))
        | ((phase_values > -0.18) & (phase_values < -0.08))
    )
    center = np.nanmedian(flux[continuum])
    scatter = robust_scatter(flux[continuum])
    eclipse = (np.abs(phase_values) < 0.05) | (np.abs(phase_values - 0.5) < 0.05)
    keep = eclipse | (np.abs(flux - center) < clip_sigma * scatter)
    clipped = rows[keep]
    path = OUTPUT_DIR / f"{label}_lightcurve.txt"
    np.savetxt(path, clipped, header="Time Flux Flux_Err Phase N NSector", fmt="%.8e")
    with (OUTPUT_DIR / f"{label}_clip_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "source_bins": int(len(rows)),
                "kept_bins": int(len(clipped)),
                "removed_bins": int(np.sum(~keep)),
                "clip_sigma": clip_sigma,
                "continuum_median": float(center),
                "continuum_robust_scatter": float(scatter),
                "removed_phases": phase_values[~keep].tolist(),
            },
            handle,
            indent=2,
        )
        handle.write("\n")
    return path, clipped


def make_variants(data):
    variants = {
        "sector_consensus_min1": bin_raw(data, "sector_consensus_min1", require_sector_count=1),
        "sector_consensus_min2": bin_raw(data, "sector_consensus_min2", require_sector_count=2),
        "sector_consensus_min3": bin_raw(data, "sector_consensus_min3", require_sector_count=3),
        "no_sector37": bin_raw(data, "no_sector37", sector_filter=lambda d, p, s: s != 37, require_sector_count=1),
        "no_sector10": bin_raw(data, "no_sector10", sector_filter=lambda d, p, s: s != 10, require_sector_count=1),
        "no_sector11": bin_raw(data, "no_sector11", sector_filter=lambda d, p, s: s != 11, require_sector_count=1),
        "primary_no_sector37": bin_raw(
            data,
            "primary_no_sector37",
            sector_filter=lambda d, p, s: ~((np.abs(p) < 0.035) & (s == 37)),
            require_sector_count=1,
        ),
    }
    variants["sector_consensus_min3_clip4"] = sigma_clip_rows(
        "sector_consensus_min3_clip4",
        variants["sector_consensus_min3"][1],
        clip_sigma=4.0,
    )
    variants["sector_consensus_min3_clip3"] = sigma_clip_rows(
        "sector_consensus_min3_clip3",
        variants["sector_consensus_min3"][1],
        clip_sigma=3.0,
    )
    return variants


def build_bundle_from_rows(rows, ntriangles=1800):
    b = base.build_bundle(rows[:, 0], rows[:, 1], rows[:, 2])
    b.set_value_all("ntriangles", ntriangles)
    return b


def candidate_set():
    targeted = json.loads(
        (base.REPO_DIR / "outputs" / "j1122_targeted_eclipse_grid" / "j1122_targeted_eclipse_grid_result.json").read_text()
    )
    start = morph.load_start()
    for key in ["ld_primary", "ld_secondary", "ecc", "per0"]:
        start.setdefault(key, morph.START[key])
    return {
        "start": start,
        "targeted_best_chi2": targeted["best_final_chi2"]["candidate"],
        "targeted_best_objective": targeted["best_final_objective"]["candidate"],
    }


def evaluate_variant(label, rows, candidates):
    b = build_bundle_from_rows(rows)
    results = []
    for cname, candidate in candidates.items():
        model = f"{label}_{cname}"
        candidate = dict(candidate)
        candidate.setdefault("ld_primary", 0.5)
        candidate.setdefault("ld_secondary", 0.5)
        candidate.setdefault("ecc", 0.0)
        candidate.setdefault("per0", 0.0)
        result = morph.evaluate(b, candidate, model, keep_model=True)
        result["variant"] = label
        result["candidate_name"] = cname
        results.append(result)
        base.save_plot(b, result, OUTPUT_DIR, model)
        print(
            f"{label} {cname}: chi2={result['chi2']:.2f} rms={result['rms']:.6f} "
            f"primary_abs={result['primary_mean_abs']:.6f} trend={result['trend_score']:.6f}",
            flush=True,
        )
    best = min(results, key=lambda r: r["objective"])
    b.save(str(OUTPUT_DIR / f"{label}_best.phoebe"))
    return results, best


def plot_variant_overview(variant_rows):
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, (_, rows) in variant_rows.items():
        ax.plot(rows[:, 3], rows[:, 1], lw=1.2, label=label)
    ax.set_xlim(-0.08, 0.08)
    ax.set_xlabel("Phase")
    ax.set_ylabel("Flux")
    ax.set_title(f"{base.TARGET}: primary eclipse robust binning comparison")
    ax.legend(fontsize=7, ncols=2)
    fig.savefig(OUTPUT_DIR / "robust_binning_primary_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_raw()
    variants = make_variants(data)
    plot_variant_overview(variants)
    candidates = candidate_set()
    all_results = []
    best_by_variant = {}
    for label, (_, rows) in variants.items():
        results, best = evaluate_variant(label, rows, candidates)
        all_results.extend(results)
        best_by_variant[label] = best
    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": base.TARGET,
        "raw_lc": str(RAW_LC.relative_to(base.REPO_DIR)),
        "variants": {label: {"path": str(path.relative_to(base.REPO_DIR)), "n_bins": int(len(rows))} for label, (path, rows) in variants.items()},
        "best_by_variant": best_by_variant,
        "best_overall": min(all_results, key=lambda r: r["objective"]),
        "all_results": all_results,
        "status": "success",
    }
    with (OUTPUT_DIR / "j1122_robust_binning_phoebe_result.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
