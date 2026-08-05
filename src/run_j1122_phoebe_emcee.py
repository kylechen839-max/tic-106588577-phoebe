import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import phoebe


REPO_DIR = Path(__file__).resolve().parents[1]
TARGET = "J112238.89-592027.5"
INPUT_LC = (
    REPO_DIR
    / "outputs"
    / "j1122_tess_multisector_s10_s11_s37"
    / f"{TARGET}_LightCurve_multisector_binned.txt"
)
INPUT_LC = Path(os.environ.get("J1122_INPUT_LC", str(INPUT_LC)))
if not INPUT_LC.is_absolute():
    INPUT_LC = REPO_DIR / INPUT_LC
OUTPUT_DIR = REPO_DIR / "outputs" / "j1122_phoebe_emcee"
PERIOD = 4.790978516308605
T0 = 1573.0869


def save_autofig(plot_result, path):
    fig = plot_result[0] if isinstance(plot_result, tuple) else plot_result
    if hasattr(fig, "save"):
        fig.save(str(path))
    else:
        plt.savefig(path, dpi=220, bbox_inches="tight")
        plt.close("all")


def residual_metrics(b, model):
    residuals = b.calculate_residuals(model=model, dataset="lc01")
    values = np.asarray(residuals.value if hasattr(residuals, "value") else residuals)
    return {
        "sum_squares": float(np.sum(values**2)),
        "rms": float(np.sqrt(np.mean(values**2))),
        "mean_abs": float(np.mean(np.abs(values))),
        "max_abs": float(np.max(np.abs(values))),
        "n_points": int(values.size),
    }


def build_bundle():
    raw = np.genfromtxt(INPUT_LC).T
    times = raw[0]
    fluxes = raw[1] / np.nanmedian(raw[1])
    sigmas = raw[2] / np.nanmedian(raw[1])
    sigmas = np.where(np.isfinite(sigmas) & (sigmas > 0), sigmas, np.nanmedian(sigmas[sigmas > 0]))

    b = phoebe.default_binary()
    b.add_dataset(
        "lc",
        times=times,
        fluxes=fluxes,
        sigmas=sigmas,
        dataset="lc01",
        passband="TESS:T",
    )
    b.set_value_all("atm", "blackbody")
    b.set_value_all("ld_mode_bol", "manual")
    b.set_value_all("ld_func_bol", "linear")
    b.set_value_all("ld_coeffs_bol", [0.5])
    b.set_value_all("ld_mode", "manual")
    b.set_value_all("ld_func", "linear")
    b.set_value_all("ld_coeffs", [0.5])
    b.set_value("period@binary", PERIOD)
    b.set_value("t0_supconj@binary", T0)
    b.set_value("sma@binary", 20.0)
    b.set_value("incl@binary", 88.0)
    b.set_value("requiv@primary", 2.0)
    b.set_value("requiv@secondary", 1.5)
    b.set_value("teff@primary", 18000.0)
    b.set_value("teff@secondary", 12000.0)
    b.set_value("gravb_bol@primary", 0.9)
    b.set_value("gravb_bol@secondary", 0.9)
    b.set_value("irrad_frac_refl_bol@primary", 1.0)
    b.set_value("irrad_frac_refl_bol@secondary", 1.0)
    b.set_value("pblum_mode@lc01", "dataset-scaled")
    return b


def run_lc_geometry(b):
    b.add_solver(
        "estimator.lc_geometry",
        solver="j1122_lcgeo_solver",
        analytical_model="polyfit",
        overwrite=True,
    )
    b.run_solver(
        solver="j1122_lcgeo_solver",
        solution="j1122_lcgeo_solution",
        overwrite=True,
    )
    adopt = [
        "t0_supconj@binary@orbit@component",
        "ecc@binary@orbit@component",
        "per0@binary@orbit@component",
        "requivsumfrac@binary@orbit@component",
    ]
    try:
        b.flip_constraint("requivsumfrac@binary", solve_for="requiv@primary")
    except Exception:
        pass
    b.adopt_solution(solution="j1122_lcgeo_solution", adopt_parameters=adopt)
    b.run_compute(compute="phoebe01", model="j1122_lcgeo_model", overwrite=True)


def run_powell(b):
    b.add_solver("optimizer.powell", solver="j1122_powell_solver", overwrite=True)
    b.set_value(
        "fit_parameters@j1122_powell_solver",
        [
            "incl@binary",
            "requivsumfrac@binary",
            "t0_supconj@binary",
        ],
    )
    b.set_value("maxiter@j1122_powell_solver", int(os.environ.get("J1122_POWELL_MAXITER", "80")))
    b.set_value("xtol@j1122_powell_solver", 1e-3)
    b.set_value("ftol@j1122_powell_solver", 1e-3)
    b.run_solver(
        solver="j1122_powell_solver",
        solution="j1122_powell_solution",
        overwrite=True,
    )
    b.adopt_solution(solution="j1122_powell_solution")
    b.run_compute(compute="phoebe01", model="j1122_after_powell", overwrite=True)


def create_baseline_without_powell(b):
    b.run_compute(compute="phoebe01", model="j1122_after_lcgeo", overwrite=True)


def run_emcee_probe(b):
    params = ["incl@binary", "t0_supconj@binary"]
    b.add_distribution(
        {
            "incl@binary": phoebe.gaussian_around(float(os.environ.get("J1122_EMCEE_INCL_SIGMA", "0.01"))),
            "t0_supconj@binary": phoebe.gaussian_around(float(os.environ.get("J1122_EMCEE_T0_SIGMA", "0.0005"))),
        },
        distribution="j1122_ball_around_powell",
        overwrite_all=True,
    )
    b.add_solver(
        "sampler.emcee",
        solver="j1122_emcee_solver",
        init_from="j1122_ball_around_powell",
        compute="phoebe01",
        overwrite=True,
    )
    nwalkers = int(os.environ.get("J1122_EMCEE_NWALKERS", "8"))
    niters = int(os.environ.get("J1122_EMCEE_NITERS", "20"))
    b.run_solver(
        solver="j1122_emcee_solver",
        solution="j1122_emcee_solution",
        nwalkers=nwalkers,
        niters=niters,
        progress_every_niters=5,
        overwrite=True,
    )
    samples = np.asarray(
        b.filter(context="solution", solution="j1122_emcee_solution", qualifier="samples").get_value(),
        dtype=float,
    )
    lnprobs = np.asarray(
        b.filter(context="solution", solution="j1122_emcee_solution", qualifier="lnprobabilities").get_value(),
        dtype=float,
    )
    fitted_twigs = list(
        b.filter(context="solution", solution="j1122_emcee_solution", qualifier="fitted_twigs").get_value()
    )
    finite = np.isfinite(lnprobs)
    if not np.any(finite):
        raise RuntimeError("No finite emcee samples found for J1122.")
    best_flat = int(np.nanargmax(np.where(finite, lnprobs, -np.inf)))
    best_iter, best_walker = np.unravel_index(best_flat, lnprobs.shape)
    for twig, value in zip(fitted_twigs, samples[best_iter, best_walker, :]):
        b.set_value(twig, float(value))
    b.run_compute(compute="phoebe01", model="j1122_emcee_best", overwrite=True)
    return {
        "params": params,
        "nwalkers": nwalkers,
        "niters": niters,
        "best_iter": int(best_iter),
        "best_walker": int(best_walker),
        "best_lnprob": float(lnprobs[best_iter, best_walker]),
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": TARGET,
        "input_lc": str(INPUT_LC.relative_to(REPO_DIR) if INPUT_LC.is_relative_to(REPO_DIR) else INPUT_LC),
        "period": PERIOD,
        "t0": T0,
        "status": "started",
    }
    b = build_bundle()
    try:
        b.run_compute(compute="phoebe01", model="j1122_initial", overwrite=True)
        run_lc_geometry(b)
        if os.environ.get("J1122_SKIP_POWELL", "false").lower() in ("1", "true", "yes"):
            create_baseline_without_powell(b)
            baseline_model = "j1122_after_lcgeo"
        else:
            run_powell(b)
            baseline_model = "j1122_after_powell"
        result["baseline_model"] = baseline_model
        result["baseline_chi2"] = float(b.calculate_chi2(model=baseline_model, dataset="lc01"))
        result["baseline_residuals"] = residual_metrics(b, baseline_model)
        result["emcee"] = run_emcee_probe(b)
        result["emcee_best_chi2"] = float(b.calculate_chi2(model="j1122_emcee_best", dataset="lc01"))
        result["emcee_best_residuals"] = residual_metrics(b, "j1122_emcee_best")
        result["status"] = "success"
        save_autofig(
            b.plot(model=[baseline_model, "j1122_emcee_best"], dataset="lc01", x="phases:primary", show=False),
            OUTPUT_DIR / "j1122_phoebe_emcee_comparison.png",
        )
        save_autofig(
            b.plot(model="j1122_emcee_best", dataset="lc01", x="phases:primary", y="residuals", show=False),
            OUTPUT_DIR / "j1122_emcee_residuals.png",
        )
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = repr(exc)
    finally:
        bundle_name = "j1122_phoebe_emcee.phoebe" if result["status"] == "success" else "j1122_phoebe_failed.phoebe"
        b.save(str(OUTPUT_DIR / bundle_name))
        with (OUTPUT_DIR / "j1122_phoebe_emcee_result.json").open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
