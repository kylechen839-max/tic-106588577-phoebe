import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


TARGET_NAME = "J071951.40-240400.6"
SEARCH_STRING = "07 19 51.40 -24 04 00.6"
PERIOD_DAYS = 0.5055155277660734 * 2
T0_SUPCONJ = 1492.5578767608158


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract and combine TESSCut light curves for TIC 106588577's "
            "source position using the Lightkurve workflow from the notebook."
        )
    )
    parser.add_argument("--target-name", default=TARGET_NAME)
    parser.add_argument("--search-string", default=SEARCH_STRING)
    parser.add_argument("--output-dir", default="outputs/tess_multisector")
    parser.add_argument("--cutout-size", type=int, default=20)
    parser.add_argument("--target-threshold", type=float, default=10.0)
    parser.add_argument("--reference-pixel", default="12,12")
    parser.add_argument("--background-y", default="18:24")
    parser.add_argument("--background-x", default="1:12")
    parser.add_argument("--pca-components", type=int, default=5)
    parser.add_argument(
        "--sector",
        action="append",
        type=int,
        help=(
            "TESS sector number to include. Repeat for multiple sectors. "
            "If omitted, all search results are processed."
        ),
    )
    parser.add_argument(
        "--known-sectors",
        default="7,33,34,61,87,88",
        help=(
            "Comma-separated sector list from the saved notebook output. "
            "Used when --sector is omitted and --skip-broad-search is set."
        ),
    )
    parser.add_argument(
        "--skip-broad-search",
        action="store_true",
        help="Search each known sector directly instead of first running a broad TESScut search.",
    )
    parser.add_argument(
        "--sigma-clip",
        type=float,
        default=6.0,
        help="Reject points more than this many MAD-scaled sigma from each sector median.",
    )
    parser.add_argument("--period", type=float, default=PERIOD_DAYS)
    parser.add_argument("--t0", type=float, default=T0_SUPCONJ)
    parser.add_argument(
        "--phase-bin-count",
        type=int,
        default=0,
        help="If positive, also save a phase-binned PHOEBE input with this many bins.",
    )
    return parser.parse_args()


def parse_pair(value):
    left, right = value.split(",", 1)
    return int(left), int(right)


def parse_slice(value):
    start, stop = value.split(":", 1)
    return slice(int(start), int(stop))


def parse_int_list(value):
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def sector_from_mission(mission):
    match = re.search(r"Sector\s+(\d+)", str(mission))
    return int(match.group(1)) if match else None


def finite_lightcurve(lc):
    good = np.isfinite(lc.time.value) & np.isfinite(lc.flux.value)
    if getattr(lc, "flux_err", None) is not None:
        good &= np.isfinite(lc.flux_err.value)
    return lc[good]


def mad_clip_mask(flux, sigma_clip):
    median = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - median))
    if not np.isfinite(mad) or mad == 0:
        return np.isfinite(flux)
    robust_sigma = 1.4826 * mad
    return np.abs(flux - median) < sigma_clip * robust_sigma


def build_background_mask(shape, y_slice, x_slice):
    mask = np.zeros(shape, dtype=bool)
    y0 = max(y_slice.start or 0, 0)
    y1 = min(y_slice.stop or shape[0], shape[0])
    x0 = max(x_slice.start or 0, 0)
    x1 = min(x_slice.stop or shape[1], shape[1])
    mask[y0:y1, x0:x1] = True
    return mask


def process_pixel_file(pf, args, sector):
    import lightkurve as lk

    reference_pixel = parse_pair(args.reference_pixel)
    background_y = parse_slice(args.background_y)
    background_x = parse_slice(args.background_x)

    target_mask = pf.create_threshold_mask(
        threshold=args.target_threshold,
        reference_pixel=reference_pixel,
    )
    background_mask = build_background_mask(
        pf.shape[1:],
        background_y,
        background_x,
    )

    lc_raw = finite_lightcurve(pf.to_lightcurve(aperture_mask=target_mask))
    design_matrix = (
        lk.DesignMatrix(pf.flux[:, background_mask], name="background_regressors")
        .pca(args.pca_components)
        .append_constant()
    )
    corrected = lk.RegressionCorrector(lc_raw).correct(design_matrix)
    corrected = finite_lightcurve(corrected)

    time = np.asarray(corrected.time.value, dtype=float)
    flux = np.asarray(corrected.flux.value, dtype=float)
    if getattr(corrected, "flux_err", None) is not None:
        flux_err = np.asarray(corrected.flux_err.value, dtype=float)
    else:
        flux_err = np.full_like(flux, np.nanstd(flux))

    keep = mad_clip_mask(flux, args.sigma_clip)
    time = time[keep]
    flux = flux[keep]
    flux_err = flux_err[keep]

    sector_median = np.nanmedian(flux)
    norm_flux = flux / sector_median
    norm_flux_err = flux_err / sector_median
    sectors = np.full_like(time, sector, dtype=float)
    return np.column_stack((time, norm_flux, norm_flux_err, sectors))


def phase_values(time, period, t0):
    phase = ((time - t0) / period) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def save_binned_lightcurve(data, args):
    time, flux, flux_err, _ = data.T
    phase = phase_values(time, args.period, args.t0)
    bins = np.linspace(-0.2, 0.8, args.phase_bin_count + 1)
    rows = []
    for lower, upper in zip(bins[:-1], bins[1:]):
        mask = (phase >= lower) & (phase < upper)
        if not np.any(mask):
            continue
        bin_phase = 0.5 * (lower + upper)
        bin_flux = np.nanmedian(flux[mask])
        scatter = 1.4826 * np.nanmedian(np.abs(flux[mask] - bin_flux))
        if not np.isfinite(scatter) or scatter == 0:
            scatter = np.nanmedian(flux_err[mask])
        bin_err = scatter / np.sqrt(np.sum(mask))
        bin_time = args.t0 + bin_phase * args.period
        rows.append((bin_time, bin_flux, bin_err, bin_phase, np.sum(mask)))

    binned = np.asarray(rows, dtype=float)
    output_path = os.path.join(
        args.output_dir,
        f"{args.target_name}_LightCurve_multisector_binned.txt",
    )
    np.savetxt(
        output_path,
        binned,
        header="Time\tFlux\tFlux_Err\tPhase\tN",
        fmt=["%.8e", "%.8e", "%.8e", "%.8e", "%.0f"],
        delimiter="\t",
    )
    return output_path, binned


def save_diagnostic(data, output_dir, target_name, period, t0, binned=None):
    time, flux, _, sector = data.T
    phase = phase_values(time, period, t0)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    for sector_number in sorted(set(sector.astype(int))):
        mask = sector.astype(int) == sector_number
        axes[0].scatter(
            time[mask],
            flux[mask],
            s=8,
            alpha=0.55,
            linewidths=0,
            label=f"S{sector_number}",
        )
        axes[1].scatter(
            phase[mask],
            flux[mask],
            s=8,
            alpha=0.55,
            linewidths=0,
            label=f"S{sector_number}",
        )
    axes[0].set_xlabel("BTJD")
    axes[0].set_ylabel("Sector-normalized flux")
    axes[1].set_xlabel("Phase relative to primary eclipse")
    axes[1].set_ylabel("Sector-normalized flux")
    axes[1].set_xlim(-0.30, 0.80)
    if binned is not None:
        axes[1].plot(
            binned[:, 3],
            binned[:, 1],
            color="black",
            lw=1.8,
            label="binned median",
        )
    axes[0].legend(ncols=3, fontsize=8)
    fig.suptitle(f"{target_name}: combined TESS sectors")
    fig.savefig(
        os.path.join(output_dir, f"{target_name}_multisector_diagnostic.png"),
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    import lightkurve as lk

    selected = []
    if args.skip_broad_search:
        sectors = args.sector or parse_int_list(args.known_sectors)
        for sector in sectors:
            selected.append((None, sector))
    else:
        search_results = lk.search_tesscut(args.search_string)
        for index, row in enumerate(search_results.table):
            sector = sector_from_mission(row["mission"])
            if args.sector and sector not in args.sector:
                continue
            selected.append((index, sector))

    if not selected:
        raise RuntimeError("No matching TESSCut sectors found.")

    all_data = []
    for index, sector in selected:
        if index is None:
            print(f"Searching, downloading, and processing sector {sector}")
            sector_results = lk.search_tesscut(args.search_string, sector=sector)
            if len(sector_results) == 0:
                print(f"  no TESScut result for sector {sector}; skipping")
                continue
            pf = sector_results[0].download(cutout_size=args.cutout_size)
        else:
            print(f"Downloading and processing sector {sector} from search index {index}")
            pf = search_results[index].download(cutout_size=args.cutout_size)
        sector_data = process_pixel_file(pf, args, sector)
        print(f"  kept {len(sector_data)} points")
        all_data.append(sector_data)

    if not all_data:
        raise RuntimeError("No sectors were successfully extracted.")

    combined = np.vstack(all_data)
    combined = combined[np.argsort(combined[:, 0])]

    output_path = os.path.join(args.output_dir, f"{args.target_name}_LightCurve_multisector.txt")
    np.savetxt(
        output_path,
        combined,
        header="Time\tFlux\tFlux_Err\tSector",
        fmt=["%.8e", "%.8e", "%.8e", "%.0f"],
        delimiter="\t",
    )
    binned = None
    if args.phase_bin_count > 0:
        binned_path, binned = save_binned_lightcurve(combined, args)
        print(f"Saved binned light curve: {binned_path}")
    save_diagnostic(
        combined,
        args.output_dir,
        args.target_name,
        args.period,
        args.t0,
        binned=binned,
    )
    print(f"Saved combined light curve: {output_path}")


if __name__ == "__main__":
    main()
