import json
import os
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import phoebe


BASELINE_CHI2 = 127974.04617631363


def solution_value(b, qualifier, solution="emcee_solution"):
    return b.filter(
        context="solution",
        solution=solution,
        qualifier=qualifier,
    ).get_value()


def safe_flip_constraint(b, constraint, solve_for):
    try:
        b.flip_constraint(constraint, solve_for=solve_for)
        print(f"Flipped {constraint} to solve for {solve_for}")
    except Exception as exc:
        print(f"Skipped flipping {constraint}: {exc}")


def parse_list(name, default):
    value = os.environ.get(name)
    if not value:
        return default
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_sigmas(params):
    sigmas = {
        "incl@binary": float(os.environ.get("EMCEE_INCL_SIGMA", "0.002")),
        "t0_supconj@binary": float(os.environ.get("EMCEE_T0_SIGMA", "0.00001")),
        "requiv@primary": float(os.environ.get("EMCEE_REQUIV_PRIMARY_SIGMA", "0.001")),
        "requivsumfrac@binary": float(
            os.environ.get("EMCEE_REQUIVSUMFRAC_SIGMA", "0.0002")
        ),
        "teffratio@binary": float(os.environ.get("EMCEE_TEFFRATIO_SIGMA", "0.0002")),
    }
    missing = [param for param in params if param not in sigmas]
    if missing:
        raise ValueError(f"No sigma configured for {missing}")
    return {param: sigmas[param] for param in params}


def apply_secondary_fill_factor(b):
    secondary_fill_factor = float(
        os.environ.get("EMCEE_SECONDARY_FILL_FACTOR", "1.0")
    )
    if secondary_fill_factor >= 1.0:
        return secondary_fill_factor

    semidetached_constraints = b.filter(
        context="constraint",
        constraint_func="semidetached",
    ).twigs
    for constraint_twig in semidetached_constraints:
        b.remove_constraint(constraint_twig)

    secondary_requiv_max = b.get_value("requiv_max@secondary@star@component")
    b.set_value(
        "requiv@secondary@star@component",
        secondary_fill_factor * secondary_requiv_max,
    )
    print(
        "Sampling with detached secondary fill factor "
        f"{secondary_fill_factor:.8f}"
    )
    return secondary_fill_factor


def phase_data(b):
    period = float(b.get_value("period@binary@orbit@component"))
    t0 = float(b.get_value("t0_supconj@binary@orbit@component"))
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    fluxes = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    phase = ((times - t0) / period) % 1.0
    phase = np.where(phase > 0.8, phase - 1.0, phase)
    return phase, fluxes


def save_best_plot(b, output_dir, model="emcee_best"):
    phase, obs_fluxes = phase_data(b)
    model_fluxes = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    residuals = obs_fluxes - model_fluxes
    order = np.argsort(phase)

    fig, (ax, rx) = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    ax.scatter(phase, obs_fluxes, s=9, color="0.25", alpha=0.35, linewidths=0)
    ax.plot(phase[order], model_fluxes[order], color="#dc2626", lw=2.4)
    rx.axhline(0, color="0.35", lw=1)
    rx.scatter(phase, residuals, s=8, color="0.25", alpha=0.35, linewidths=0)
    rx.plot(phase[order], residuals[order], color="#dc2626", lw=1.0)
    ax.set_ylabel("Flux")
    rx.set_ylabel("Obs - model")
    rx.set_xlabel("Phase relative to primary eclipse")
    ax.set_xlim(-0.30, 0.80)
    ax.set_title("TIC 106588577: best emcee sample model")
    fig.savefig(
        os.path.join(output_dir, "emcee_best_diagnostic_lightcurve.png"),
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(script_dir)
    output_dir = os.environ.get("OUTPUT_DIR", os.path.join(repo_dir, "outputs"))
    bundle_path = os.environ.get(
        "POWELL_BUNDLE_PATH",
        os.path.join(repo_dir, "bundles", "tic_106588577_powell_ready.phoebe"),
    )
    os.makedirs(output_dir, exist_ok=True)

    b = phoebe.load(bundle_path)
    b.run_compute(compute="phoebe01", model="after_powell_check", overwrite=True)
    after_powell_chi2 = b.calculate_chi2(model="after_powell_check", dataset="lc01")
    secondary_fill_factor = apply_secondary_fill_factor(b)

    emcee_params = parse_list(
        "EMCEE_PARAMS",
        ["incl@binary", "t0_supconj@binary"],
    )
    if "requivsumfrac@binary" in emcee_params:
        safe_flip_constraint(
            b=b,
            constraint="requivsumfrac@binary",
            solve_for="requiv@primary",
        )
    if "teffratio@binary" in emcee_params:
        safe_flip_constraint(
            b=b,
            constraint="teffratio@binary",
            solve_for="teff@secondary",
        )

    sigmas = parse_sigmas(emcee_params)
    print("Experiment params:", emcee_params)
    print("Experiment sigmas:", sigmas)

    b.add_distribution(
        {
            param: phoebe.gaussian_around(sigma)
            for param, sigma in sigmas.items()
        },
        distribution="ball_around_powell",
        overwrite_all=True,
    )
    b.add_solver(
        "sampler.emcee",
        solver="emcee_solver",
        init_from="ball_around_powell",
        compute="phoebe01",
        overwrite=True,
    )
    nwalkers = int(os.environ.get("EMCEE_NWALKERS", "16"))
    niters = int(os.environ.get("EMCEE_NITERS", "100"))
    burnin = int(os.environ.get("EMCEE_BURNIN", "0"))

    b.run_solver(
        solver="emcee_solver",
        solution="emcee_solution",
        nwalkers=nwalkers,
        niters=niters,
        progress_every_niters=int(os.environ.get("EMCEE_PROGRESS", "25")),
        overwrite=True,
    )

    samples = np.asarray(solution_value(b, "samples"), dtype=float)
    lnprobs = np.asarray(solution_value(b, "lnprobabilities"), dtype=float)
    fitted_twigs = list(solution_value(b, "fitted_twigs"))

    finite = np.isfinite(lnprobs)
    if burnin > 0:
        finite[:burnin, :] = False
    if not np.any(finite):
        try:
            failed_samples = solution_value(b, "failed_samples")
        except Exception as exc:
            failed_samples = {"diagnostic_error": str(exc)}
        failure_path = os.path.join(output_dir, "emcee_experiment_failure.json")
        with open(failure_path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                    "params": emcee_params,
                    "sigmas": sigmas,
        "nwalkers": nwalkers,
        "niters": niters,
        "secondary_fill_factor": secondary_fill_factor,
        "finite_lnprob_count": int(np.sum(finite)),
                    "lnprob_min": float(np.nanmin(lnprobs)),
                    "lnprob_max": float(np.nanmax(lnprobs)),
                    "failed_samples": failed_samples,
                },
                handle,
                indent=2,
                default=str,
            )
            handle.write("\n")
        b.save(os.path.join(output_dir, "tic_106588577_emcee_failed.phoebe"))
        print(f"Saved failure diagnostics: {failure_path}")
        raise RuntimeError("No finite emcee samples found.")

    best_flat_index = int(np.nanargmax(np.where(finite, lnprobs, -np.inf)))
    best_iter, best_walker = np.unravel_index(best_flat_index, lnprobs.shape)
    best_values = samples[best_iter, best_walker, :]
    best_lnprob = float(lnprobs[best_iter, best_walker])

    for twig, value in zip(fitted_twigs, best_values):
        b.set_value(twig, float(value))
        print(f"Best sample {twig} = {float(value)}")

    b.run_compute(compute="phoebe01", model="emcee_best", overwrite=True)
    emcee_best_chi2 = b.calculate_chi2(model="emcee_best", dataset="lc01")
    save_best_plot(b, output_dir, model="emcee_best")

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "bundle_path": bundle_path,
        "params": emcee_params,
        "sigmas": sigmas,
        "nwalkers": nwalkers,
        "niters": niters,
        "secondary_fill_factor": secondary_fill_factor,
        "burnin": burnin,
        "best_iter": int(best_iter),
        "best_walker": int(best_walker),
        "best_lnprob": best_lnprob,
        "after_powell_chi2": float(after_powell_chi2),
        "emcee_best_chi2": float(emcee_best_chi2),
        "beats_after_powell": bool(emcee_best_chi2 < after_powell_chi2),
        "target_after_powell_chi2": BASELINE_CHI2,
    }
    result_path = os.path.join(output_dir, "emcee_experiment_result.json")
    with open(result_path, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")

    b.save(os.path.join(output_dir, "tic_106588577_emcee_best.phoebe"))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
