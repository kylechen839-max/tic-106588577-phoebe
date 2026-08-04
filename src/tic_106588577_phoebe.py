import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import phoebe
from phoebe import u  # units


def safe_flip_constraint(b, constraint, solve_for):
    """
    Flip a PHOEBE constraint only if possible.
    If it was already flipped in a previous run, skip cleanly.
    """
    try:
        b.flip_constraint(constraint, solve_for=solve_for)
        print(f"Flipped {constraint} to solve for {solve_for}")
    except Exception as e:
        print(f"Skipped flipping {constraint}: {e}")


def save_plot(plot_result, filename):
    fig = plot_result[0] if isinstance(plot_result, tuple) else plot_result
    if hasattr(fig, "save"):
        fig.save(filename)
    else:
        plt.savefig(filename, dpi=200, bbox_inches="tight")
        plt.close("all")
    print(f"Saved plot: {filename}")


def build_initial_bundle(data_path):
    logger = phoebe.logger(clevel="WARNING")

    raw_data = np.genfromtxt(data_path).T
    times = raw_data[0]
    raw_fluxes = raw_data[1]
    raw_sigmas = raw_data[2]

    # Calculate relative flux and sigma to pass to dataset.
    median_flux = np.median(raw_fluxes)
    fluxes = raw_fluxes / median_flux
    sigmas = raw_sigmas / median_flux

    b = phoebe.default_binary(semidetached="secondary")
    b.add_dataset(
        "lc",
        times=times,
        fluxes=fluxes,
        sigmas=sigmas,
        dataset="lc01",
        passband="TESS:T",
    )

    # Approximations and guesses.
    b.set_value_all("atm", "blackbody")
    b.set_value_all("ld_mode_bol", "manual")
    b.set_value_all("ld_func_bol", "linear")
    b.set_value_all("ld_coeffs_bol", [0.5])
    b.set_value_all("ld_mode", "manual")
    b.set_value_all("ld_func", "linear")
    b.set_value_all("ld_coeffs", [0.5])
    b.set_value(qualifier="teff", component="primary", value=17000)
    b.set_value("gravb_bol@primary", 0.9)
    b.set_value("gravb_bol@secondary", 0.9)
    b.set_value("irrad_frac_refl_bol@primary", 1.0)
    b.set_value("irrad_frac_refl_bol@secondary", 1.0)

    b.set_value("pblum_mode", "dataset-scaled")

    # Known information.
    b["period@binary"] = float(os.environ.get("PERIOD_DAYS", 0.5055155277660734 * 2))
    b["teff@primary"] = 21000
    b["teff@secondary"] = 15138
    b["requiv@primary"] = 2

    # Calculating semi-major axis.
    sol = 1.98 * 10**30
    sma = (
        ((6.67 * 10**-11) * ((8 + 6) * sol) / (4 * (np.pi**2))) * (86400**2)
    ) ** (1 / 3)
    b["sma@binary"] = sma / (6.957 * 10**8)

    print(b["requiv@secondary@component"])
    print(b["requiv_max@secondary@component"])

    b.run_compute()
    return b


def run_lc_geometry(b, output_dir):
    print("Running light-curve geometry estimator...")

    b.add_solver(
        "estimator.lc_geometry",
        solver="lcgeo_solver",
        analytical_model="polyfit",
        overwrite=True,
    )
    b.run_solver(solver="lcgeo_solver", solution="lcgeo_solution", overwrite=True)

    safe_flip_constraint(
        b=b, constraint="requivsumfrac", solve_for="requiv@primary"
    )
    safe_flip_constraint(b=b, constraint="teffratio", solve_for="teff@secondary")

    print(b.filter(context="solution", solution="lcgeo_solution"))

    adopt_list = [
        "t0_supconj@binary@orbit@component",
        "ecc@binary@orbit@component",
        "per0@binary@orbit@component",
        "requivsumfrac@binary@orbit@component",
    ]

    b.adopt_solution(solution="lcgeo_solution", adopt_parameters=adopt_list)
    b.run_compute(model="lcgeo_model", overwrite=True)

    plot_result = b.plot(
        dataset="lc01",
        model="lcgeo_model",
        x="phases:primary",
        show=False,
    )
    save_plot(plot_result, os.path.join(output_dir, "lcgeo_model.png"))

    print("t0_supconj:", b.get_value("t0_supconj@binary@orbit"))
    print("ecc:", b.get_value("ecc@binary"))
    print("per0:", b.get_value("per0@binary"))
    print("requivsumfrac:", b.get_value("requivsumfrac@binary"))


def run_first_powell(b, output_dir):
    print("Running first Powell optimizer...")

    safe_flip_constraint(
        b=b, constraint="requiv@primary", solve_for="requivsumfrac@binary"
    )
    safe_flip_constraint(
        b=b, constraint="teff@secondary", solve_for="teffratio@binary"
    )

    b.set_value("pblum_mode@lc01", "dataset-scaled")

    b.add_solver("optimizer.powell", solver="powell_solver", overwrite=True)

    b.set_value(
        "fit_parameters@powell_solver",
        [
            "incl@binary",
            "requiv@primary",
            "t0_supconj@binary",
        ],
    )

    b.set_value("maxiter@powell_solver", 200)
    b.set_value("xtol@powell_solver", 1e-3)
    b.set_value("ftol@powell_solver", 1e-3)

    b.run_solver("powell_solver", solution="powell_solution", overwrite=True)

    b.adopt_solution(solution="powell_solution")
    b.run_compute(model="after_powell", overwrite=True)

    plot_result = b.plot(
        model="after_powell",
        x="phases:primary",
        linestyle={"model": "solid"},
        legend=True,
        show=False,
    )
    save_plot(plot_result, os.path.join(output_dir, "after_powell.png"))

    print("Inclination:", b.get_value("incl@binary@orbit@component"))
    print("Temperature of Primary Star (in K):", b.get_value("teff@primary"))
    print(
        "Temperature of Secondary Star (in K):",
        b.get_value("teff@secondary@star"),
    )
    print(
        "Equivalent Radius of Secondary (Sol):",
        b.get_value("requiv@secondary@star@component"),
    )
    print(
        "Equivalent Radius of Primary (Sol):",
        b.get_value("requiv@primary@star@component"),
    )
    print("Semi Major Axis (Sol):", b.get_value("sma@binary"))


def run_second_powell(b):
    print("Running second Powell optimizer...")

    safe_flip_constraint(
        b=b, constraint="requivsumfrac@binary", solve_for="requiv@primary"
    )
    safe_flip_constraint(
        b=b, constraint="teffratio@binary", solve_for="teff@secondary"
    )

    b.set_value(
        "fit_parameters@powell_solver",
        [
            "requivsumfrac@binary",
            "incl@binary",
            "teffratio@binary",
            "t0_supconj@binary",
        ],
    )

    b.run_solver("powell_solver", solution="powell_solution_2", overwrite=True)
    b.adopt_solution(solution="powell_solution_2")
    b.run_compute(model="after_powell_2", overwrite=True)


def compare_powell_models(b, output_dir):
    print("Comparing Powell models...")

    plot_result = b.plot(
        model=["after_powell_2", "after_powell"],
        x="phases:primary",
        linestyle={"model": "solid"},
        legend=True,
        show=False,
    )
    save_plot(plot_result, os.path.join(output_dir, "powell_model_comparison.png"))

    plot_result = b.plot(
        model=["after_powell_2", "after_powell"],
        x="phases:primary",
        linestyle={"model": "solid"},
        y="residuals",
        legend=True,
        show=False,
    )
    save_plot(plot_result, os.path.join(output_dir, "powell_residual_comparison.png"))

    print(
        "Powell 1 Chi-Squared:",
        b.calculate_chi2(model="after_powell", dataset="lc01"),
    )
    print(
        "Powell 2 Chi-Squared:",
        b.calculate_chi2(model="after_powell_2", dataset="lc01"),
    )

    r1 = b.calculate_residuals(model="after_powell")
    r2 = b.calculate_residuals(model="after_powell_2")

    r1_vals = np.asarray(r1.value if hasattr(r1, "value") else r1)
    r2_vals = np.asarray(r2.value if hasattr(r2, "value") else r2)

    print("Powell 1 RMS:", np.sqrt(np.mean(r1_vals**2)))
    print("Powell 2 RMS:", np.sqrt(np.mean(r2_vals**2)))

    print("Powell 1 mean abs:", np.mean(np.abs(r1_vals)))
    print("Powell 2 mean abs:", np.mean(np.abs(r2_vals)))


def prepare_powell_bundle(data_path, output_dir, powell_bundle_path):
    b = build_initial_bundle(data_path)

    plot_result = b.plot(dataset="lc01", x="phases:primary", show=False)
    save_plot(plot_result, os.path.join(output_dir, "initial_lightcurve.png"))

    run_lc_geometry(b, output_dir)
    run_first_powell(b, output_dir)
    run_second_powell(b)
    compare_powell_models(b, output_dir)

    b.save(powell_bundle_path)
    print(f"Saved Powell-ready bundle: {powell_bundle_path}")
    return b


def run_emcee(b, output_dir):
    print("Running emcee solver...")

    secondary_fill_factor = float(
        os.environ.get("EMCEE_SECONDARY_FILL_FACTOR", "1.0")
    )
    if secondary_fill_factor < 1.0:
        semidetached_constraints = b.filter(
            context="constraint", constraint_func="semidetached"
        ).twigs
        for constraint_twig in semidetached_constraints:
            b.remove_constraint(constraint_twig)

        secondary_requiv_max = b.get_value(
            "requiv_max@secondary@star@component"
        )
        b.set_value(
            "requiv@secondary@star@component",
            secondary_fill_factor * secondary_requiv_max,
        )
        print(
            "Sampling with a detached secondary at "
            f"{secondary_fill_factor:.6f} of its critical Roche radius."
        )

    default_emcee_params = [
        "incl@binary",
        "requivsumfrac@binary",
        "teffratio@binary",
        "t0_supconj@binary",
    ]
    emcee_params_env = os.environ.get("EMCEE_PARAMS")
    emcee_params = (
        [param.strip() for param in emcee_params_env.split(",") if param.strip()]
        if emcee_params_env
        else default_emcee_params
    )
    print("emcee parameters:", emcee_params)

    if "requivsumfrac@binary" in emcee_params:
        safe_flip_constraint(
            b=b, constraint="requivsumfrac@binary", solve_for="requiv@primary"
        )
    if "teffratio@binary" in emcee_params:
        safe_flip_constraint(
            b=b, constraint="teffratio@binary", solve_for="teff@secondary"
        )

    incl_sigma = float(os.environ.get("EMCEE_INCL_SIGMA", "0.2"))
    requivsumfrac_sigma = float(os.environ.get("EMCEE_REQUIVSUMFRAC_SIGMA", "0.01"))
    teffratio_sigma = float(os.environ.get("EMCEE_TEFFRATIO_SIGMA", "0.01"))
    t0_sigma = float(os.environ.get("EMCEE_T0_SIGMA", "0.0005"))

    emcee_sigmas = {
        "incl@binary": incl_sigma,
        "requivsumfrac@binary": requivsumfrac_sigma,
        "teffratio@binary": teffratio_sigma,
        "t0_supconj@binary": t0_sigma,
    }
    print("emcee starting widths:", emcee_sigmas)

    missing_sigmas = [param for param in emcee_params if param not in emcee_sigmas]
    if missing_sigmas:
        raise ValueError(f"No emcee starting width configured for: {missing_sigmas}")

    b.add_distribution(
        {
            param: phoebe.gaussian_around(emcee_sigmas[param])
            for param in emcee_params
        },
        distribution="ball_around_powell",
    )

    emcee_compute = os.environ.get("EMCEE_COMPUTE", "phoebe01")

    if emcee_compute == "ellcbackend":
        import ellc

        b.add_compute("ellc", compute="ellcbackend", overwrite=True)
    else:
        print(f"Using existing compute backend for emcee: {emcee_compute}")

    b.add_solver(
        "sampler.emcee",
        solver="emcee_solver",
        init_from="ball_around_powell",
        compute=emcee_compute,
        overwrite=True,
    )

    nwalkers = int(os.environ.get("EMCEE_NWALKERS", "32"))
    niters = int(os.environ.get("EMCEE_NITERS", "2000"))
    progress_every_niters = int(os.environ.get("EMCEE_PROGRESS", "100"))

    b.run_solver(
        solver="emcee_solver",
        solution="emcee_solution",
        nwalkers=nwalkers,
        niters=niters,
        progress_every_niters=progress_every_niters,
        overwrite=True,
    )

    sample_num = int(os.environ.get("EMCEE_SAMPLE_NUM", "20"))

    if sample_num > 0:
        b.run_compute(
            compute=emcee_compute,
            sample_from="emcee_solution",
            sample_num=sample_num,
            sample_mode="3-sigma",
            model="emcee_posterior_models",
            overwrite=True,
        )

        plot_result = b.plot(
            model="emcee_posterior_models",
            dataset="lc01",
            x="phases:primary",
            xlim=(-0.30, 0.80),
            show=False,
        )
        save_plot(plot_result, os.path.join(output_dir, "emcee_posterior_models.png"))
    else:
        print("Skipping posterior model compute because EMCEE_SAMPLE_NUM=0.")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.environ.get(
        "LIGHTCURVE_PATH",
        os.path.join(script_dir, "J071951.40-240400.6_LightCurve.txt"),
    )
    output_dir = os.environ.get("OUTPUT_DIR", os.path.join(script_dir, "outputs"))
    powell_bundle_path = os.environ.get(
        "POWELL_BUNDLE_PATH",
        os.path.join(script_dir, "tic_106588577_powell_ready.phoebe"),
    )
    use_powell_bundle = os.environ.get("USE_POWELL_BUNDLE", "auto").lower()

    os.makedirs(output_dir, exist_ok=True)

    print(f"Using light curve: {data_path}")
    print(f"Writing outputs to: {output_dir}")
    print(f"Powell bundle path: {powell_bundle_path}")

    if use_powell_bundle in ("1", "true", "yes") and not os.path.exists(
        powell_bundle_path
    ):
        raise FileNotFoundError(
            f"USE_POWELL_BUNDLE is true, but no bundle exists at {powell_bundle_path}"
        )

    if use_powell_bundle not in ("0", "false", "no") and os.path.exists(
        powell_bundle_path
    ):
        print("Loading existing Powell-ready bundle; skipping lc_geometry and Powell.")
        b = phoebe.load(powell_bundle_path)
    else:
        print("No Powell-ready bundle found; running lc_geometry and Powell first.")
        b = prepare_powell_bundle(data_path, output_dir, powell_bundle_path)

    if os.environ.get("SKIP_EMCEE", "false").lower() in ("1", "true", "yes"):
        print("Skipping emcee because SKIP_EMCEE=true.")
    else:
        run_emcee(b, output_dir)

    bundle_path = os.path.join(output_dir, "tic_106588577_final.phoebe")
    b.save(bundle_path)
    print(f"Saved final bundle: {bundle_path}")


if __name__ == "__main__":
    main()
