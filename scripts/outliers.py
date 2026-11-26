"""Select candidate dates and create cumulative outlier-removal inputs."""

import os
import json
import glob
import argparse
import logging
import re
import shutil
from pathlib import Path
import itertools
import logging
import geopandas as gpd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timezone

import phenology as ph
from functions import extract_phenology_metrics, set_logging, verify_arg_file, parse_args, write_provenance, smooth_cubic_spline


def get_candidates(
    lake_id: str,
    data_dir: str | Path,
    drop_max: int = 5,
    z_thresh: float = 3,
    variable: str = "chla",
    qa_variable: str = "lwlr_quality_flags",
    qa_filter: bool = True,
) -> list[pd.Timestamp]:
    """Select candidate observation dates for a lake.

    A date is eligible if at least one pixel exceeds z_thresh standard
    deviations above its own temporal mean. Eligible dates are ranked
    by their spatial observation range.

    Parameters
    ----------
    lake_id : str
        Lake identifier.
    data_dir : str or Path
        Folder containing the original extracted lake NetCDF file.
    drop_max : int
        Maximum number of candidate dates to return.
    z_thresh : float
        Strict lower threshold for unusually high z-scores.
    variable : str
        Observation variable in the NetCDF file.
    qa_variable : str
        Quality-flag variable.
    qa_filter : bool
        If True, keep only observations with quality flag 0.

    Returns
    -------
    list[pd.Timestamp]
        Candidate dates in cumulative removal order.
    """
    if not isinstance(drop_max, int) or isinstance(drop_max, bool):
        raise TypeError("drop_max must be an integer.")

    if drop_max < 1:
        raise ValueError("drop_max must be at least 1.")

    if not np.isfinite(z_thresh) or z_thresh < 0:
        raise ValueError("z_thresh must be finite and non-negative.")

    input_file = Path(data_dir) / f"{lake_id}.nc"

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    logging.info("Selecting candidate dates for lake %s", lake_id)

    with xr.open_dataset(
        input_file,
        mask_and_scale=False,
        decode_times=False,
    ) as dataset:
        # The extracted file stores time as Unix seconds.
        time = pd.DatetimeIndex(
            pd.to_datetime(dataset["time"].values, unit="s")
        )

        observations = dataset[variable].transpose("time", ...)
        data = observations.values.astype(float)

        data[(data == -9999) | ~np.isfinite(data)] = np.nan

        if qa_filter:
            flags = dataset[qa_variable].transpose(
                *observations.dims
            ).values

            if flags.shape != data.shape:
                raise ValueError(
                    "Observation and quality-flag shapes differ."
                )

            data[flags != 0] = np.nan

    if len(time) == 0:
        logging.info(
            "No observation dates available for lake %s",
            lake_id,
        )
        return []

    if time.hasnans or time.has_duplicates:
        raise ValueError(
            "Observation timestamps must be valid and unique."
        )

    # One row per observation date and one column per pixel.
    flat = data.reshape(len(time), -1)

    # Exclude pixels without valid observations.
    usable_pixels = np.any(np.isfinite(flat), axis=0)
    flat = flat[:, usable_pixels]

    if flat.shape[1] == 0:
        logging.info(
            "No valid observations available for lake %s",
            lake_id,
        )
        return []

    # Calculate temporal statistics separately for each pixel.
    mean = np.nanmean(flat, axis=0)
    standard_deviation = np.nanstd(flat, axis=0)

    sd_safe = np.where(
        standard_deviation > 0,
        standard_deviation,
        np.nan,
    )

    z = (flat - mean) / sd_safe

    # A date qualifies if at least one pixel exceeds the threshold.
    eligible = np.any(z > z_thresh, axis=1)

    if not eligible.any():
        logging.info(
            "No candidate dates exceed z_thresh=%s for lake %s",
            z_thresh,
            lake_id,
        )
        return []

    # Eligible dates contain at least one valid observation.
    eligible_values = flat[eligible]

    daily_range = (
        np.nanmax(eligible_values, axis=1)
        - np.nanmin(eligible_values, axis=1)
    )

    candidates = pd.DataFrame({
        "date": time[eligible],
        "range": daily_range,
    })

    ranked = candidates.sort_values(
        "range",
        ascending=False,
        kind="stable",
    )

    candidate_dates = list(ranked["date"].head(drop_max))

    logging.info(
        "Selected %d of at most %d candidate dates for lake %s: %s",
        len(candidate_dates),
        drop_max,
        lake_id,
        ", ".join(
            date.strftime("%Y-%m-%d")
            for date in candidate_dates
        ),
    )

    return candidate_dates


def create_updated_input_files(
    input_file: str | Path,
    candidate_dates: list[pd.Timestamp],
    output_dir: str | Path,
    variable: str = "chla",
) -> list[Path]:
    """Create cumulative input variants k1 to kN.

    The original file remains unchanged and serves as the baseline.
    Variant kN masks the first N candidate dates by setting the
    observation variable to -9999 across all pixels.

    Existing variant input files are overwritten. Existing phenology
    fits are not updated or removed.

    Parameters
    ----------
    input_file : str or Path
        Original extracted lake NetCDF file.
    candidate_dates : list[pd.Timestamp]
        Candidate dates in removal order.
    output_dir : str or Path
        Folder receiving k1_<filename>.nc to kN_<filename>.nc.
    variable : str
        Observation variable to mask.

    Returns
    -------
    list[Path]
        Created variant paths, excluding the original baseline.
    """
    input_file = Path(input_file)
    output_dir = Path(output_dir)

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    candidate_timestamps = pd.DatetimeIndex(
        [pd.Timestamp(date) for date in candidate_dates]
    )

    if (
        candidate_timestamps.hasnans
        or candidate_timestamps.has_duplicates
    ):
        raise ValueError(
            "Candidate dates must be valid and unique."
        )

    if len(candidate_timestamps) == 0:
        logging.info(
            "No candidate dates supplied; no variants created."
        )
        return []

    # Validate every candidate before creating or overwriting files.
    with netCDF4.Dataset(str(input_file), "r") as dataset:
        dataset.set_auto_mask(False)

        time = pd.DatetimeIndex(
            pd.to_datetime(
                dataset.variables["time"][:],
                unit="s",
            )
        )

        if variable not in dataset.variables:
            raise KeyError(
                f"Variable not found in input file: {variable}"
            )

        dimensions = dataset.variables[variable].dimensions

        if not dimensions or dimensions[0] != "time":
            raise ValueError(
                f"The first dimension of {variable} must be time."
            )

    candidate_indices = []

    for candidate_timestamp in candidate_timestamps:
        matches = np.flatnonzero(time == candidate_timestamp)

        if len(matches) != 1:
            raise ValueError(
                "Candidate date does not match exactly one time step: "
                f"{candidate_timestamp}"
            )

        candidate_indices.append(int(matches[0]))

    output_files = [
        output_dir / f"k{k}_{input_file.name}"
        for k in range(1, len(candidate_indices) + 1)
    ]

    for output_file in output_files:
        if output_file.resolve() == input_file.resolve():
            raise ValueError(
                "A variant must not overwrite the original file."
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    created_files = []

    for number_removed, output_file in enumerate(
        output_files,
        start=1,
    ):
        if output_file.exists():
            logging.info(
                "Overwriting variant input: %s",
                output_file,
            )

        # Each variant starts from the unchanged original file.
        shutil.copy2(input_file, output_file)

        with netCDF4.Dataset(str(output_file), "r+") as dataset:
            observations = dataset.variables[variable]

            for index in candidate_indices[:number_removed]:
                observations[index, ...] = -9999

            # Store the removal settings with the variant.
            dataset.setncattr("outlier_source_file", input_file.name)
            dataset.setncattr("outlier_removed_count", number_removed)
            dataset.setncattr("outlier_removed_dates", ", ".join(date.isoformat() for date in candidate_timestamps[:number_removed]))

        created_files.append(output_file)

        logging.info(
            "Created %s: %d candidate date(s) masked [%s]",
            output_file,
            number_removed,
            ", ".join(
                date.strftime("%Y-%m-%d")
                for date in candidate_timestamps[:number_removed]
            ),
        )

    return created_files


def extract_outliers(lake, p,  threads=1) -> list[Path]:
    """Prepare cumulative outlier-removal inputs for one lake.

    Parameters
    ----------
    lake : dict or GeoDataFrame row
        Lake metadata containing 'id'.
    p : dict
        Pipeline parameters:
        - out_folder: base output directory.
        - variable: observation variable and variable-folder name.
        - drop_max: maximum candidate count, default 5.
        - z_thresh: z-score threshold, default 3.
        - qa: quality-flag variable.
        - qa_filter: whether to accept only quality flag 0.
    threads : int
        Retained for compatibility with pipeline calls; not used here.
    batch_size : int
        Retained for compatibility with pipeline calls; not used here.

    Returns
    -------
    list[Path]
        Created variant input files. The baseline is not copied.
    """
    lake_id = str(lake["id"])
    variable = p["variable"]
    base_dir = Path(p["out_folder"])

    extract_dir = base_dir / "extract" / variable
    input_file = extract_dir / f"{lake_id}.nc"
    outliers_dir = extract_dir / "outliers"
    phenology_outliers_dir = (
        base_dir / "phenology" / variable / "outliers"
    )

    if not input_file.is_file():
        logging.info(
            "Input file not available for lake %s: %s",
            lake_id,
            input_file,
        )
        return []

    logging.info(
        "Starting outlier input preparation for lake %s",
        lake_id,
    )

    candidate_dates = get_candidates(
        lake_id=lake_id,
        data_dir=extract_dir,
        drop_max=p.get("drop_max", 5),
        z_thresh=p.get("z_thresh", 3),
        variable=variable,
        qa_variable=p.get("qa", "lwlr_quality_flags"),
        qa_filter=p.get("qa_filter", True),
    )

    if not candidate_dates:
        logging.info(
            "No candidate dates selected for lake %s; "
            "no variants created. Existing variants and fits "
            "are left unchanged.",
            lake_id,
        )
        return []

    # Manual override
    candidate_dates[0] = pd.Timestamp("2002-11-06")
    candidate_dates[1] = pd.Timestamp("2008-12-13")

    # Prepare the input and future fit-output directories.
    outliers_dir.mkdir(parents=True, exist_ok=True)
    phenology_outliers_dir.mkdir(parents=True, exist_ok=True)

    created_files = create_updated_input_files(
        input_file=input_file,
        candidate_dates=candidate_dates,
        output_dir=outliers_dir,
        variable=variable,
    )

    logging.info(
        "Completed outlier input preparation for lake %s: "
        "%d variants created",
        lake_id,
        len(created_files),
    )

    logging.info("Baseline remains unchanged: %s",input_file)

    logging.info(
        "Calculate matching phenology fits separately in %s. "
        "Existing fits are not updated by input preparation.",
        phenology_outliers_dir,
    )

    return created_files

def phenology_outliers(lake, p, threads=1, batch_size=100):
    """Extract phenology metrics for a single lake from time series data.

    Reads extracted time series, applies a smoothing cubic spline, and extracts
    phenology metrics (green-up/green-down onset, mid, advanced) for each pixel.
    Batch results are saved as .npy checkpoints so the run can be resumed if
    interrupted. The final netCDF4 is written in one sequential pass after all
    batches complete.

    Parameters
    ----------
    lake : dict or GeoDataFrame row
        Lake metadata containing at minimum 'id'.
    p : dict
        Parameters dict with keys:
        - out_folder : str - base output directory
        - variable : str - variable name (e.g. 'chla_mean')
        - qa : str - QA variable name
        - qa_filter : bool - whether to filter on QA flag
        - spline_min_phase_length : int - minimum phase length in days
        - spline_min_relative_amplitude : float - minimum relative amplitude (0-1)
        - spline_min_phase_data : int - minimum data points per phase
        - spline_subs_peak_win_size : int - window size for substantial peak check (days)
        - spline_subs_peak_ampl_frac : float - amplitude fraction for substantial peak check
        - spline_data_gap_size : int - minimum gap size to flag (days)
        - spline_data_gap_size_buffer : int - buffer around data gaps (days)
    batch_size : int
        Number of pixels per batch dispatched to each worker.
    """
    global _VAR_DATA, _QA_DATA
    
    e_outliers_dir = Path(p["out_folder"]) /  "extract"  / p["variable"] / "outliers"
    p_outliers_dir = Path(p["out_folder"]) / "phenology" / p["variable"] / "outliers"
    
    p_outliers_dir.mkdir(parents=True, exist_ok=True)

    input_files = list(Path(e_outliers_dir).glob(f"k*_{lake['id']}.nc"))
    kf = [(int(re.search(r"k(\d+)_", fp.name).group(1)), fp) for fp in input_files]

    for idx, input_file in kf:
        checkpoint_dir = Path(p_outliers_dir) / "checkpoints" / "_".join([str(idx), str(lake["id"])]) / f"bs{batch_size}"
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
                
        logging.info(f"Processing input file {input_file.name}")
        logging.info(f"Input file: {input_file}")
        output_file = p_outliers_dir / input_file.name

        if output_file.is_file():
            logging.info(f"File {output_file.name} already exists")
            continue

        output_file_temp = output_file.with_name(f"{output_file.stem}_tmp{output_file.suffix}")

        if os.path.isfile(input_file):
            logging.info(f"Starting outlier phenology computation for {lake['id']}")
        else:
            logging.info(f"Input outlier file not available for {lake['id']}, k = {idx}")
            return

        with netCDF4.Dataset(input_file) as nc:
            valid_coords = np.argwhere(nc.variables["summary"][:, :] > 2)
            lat = np.array(nc.variables["lat"])
            lon = np.array(nc.variables["lon"])
            t = np.array([datetime.utcfromtimestamp(int(ts)).toordinal() + 366 for ts in np.array(nc.variables["time"])])
            _VAR_DATA = np.array(nc.variables[p["variable"]])
            _QA_DATA = np.array(nc.variables[p["qa"]]) if p["qa_filter"] else None

        smooth_x_axis = np.arange(t.min(), t.max() + 1, 1)

        coord_list = [(int(x), int(y)) for x, y in valid_coords]
        batches = [coord_list[i:i + batch_size] for i in range(0, len(coord_list), batch_size)]

        # Find already completed batches for restart
        completed = {int(os.path.splitext(os.path.basename(f))[0].split('_')[1])
                    for f in glob.glob(os.path.join(checkpoint_dir, 'batch_*.npy'))}
        logging.info(f"Lake {lake['id']}: {len(completed)}/{len(batches)} batches already complete")

        if threads > 1:
            with ProcessPoolExecutor(max_workers=threads, initializer=ph._worker_init) as executor:
                futures = {executor.submit(compute_pixel_batch, batch, smooth_x_axis, p, t): i
                        for i, batch in enumerate(batches) if i not in completed}
                for future in tqdm(as_completed(futures), total=len(batches),
                                initial=len(completed), desc=str(lake['id']), unit='batch'):
                    batch_idx = futures[future]
                    try:
                        batch_result = future.result()
                    except Exception as e:
                        print(f"Worker failed on batch {batch_idx}: {e}", flush=True)
                        continue
                    checkpoint = [(x, y, result['smoothing'], pheno)
                                for x, y, result, pheno in batch_result if result is not None]
                    np.save(os.path.join(checkpoint_dir, f'batch_{batch_idx:04d}.npy'),
                            checkpoint, allow_pickle=True)
        else:
            for i, batch in enumerate(tqdm(batches, desc=str(lake['id']), unit='batch')):
                if i in completed:
                    continue
                batch_result = compute_pixel_batch(batch, smooth_x_axis, p, t)
                checkpoint = [(x, y, result['smoothing'], pheno)
                            for x, y, result, pheno in batch_result if result is not None]
                np.save(os.path.join(checkpoint_dir, f'batch_{i:04d}.npy'),
                        checkpoint, allow_pickle=True)

        _VAR_DATA = None
        _QA_DATA = None

        ph._write_from_checkpoints(checkpoint_dir, output_file_temp, output_file, lat, lon, p, lake["id"])
        shutil.rmtree(checkpoint_dir)


def compute_pixel_batch(coord_batch, smooth_x_axis, p, t):
    """Compute spline and phenology metrics for a batch of pixels.

    Reads pixel data from the module-level _VAR_DATA and _QA_DATA arrays
    inherited from the parent process via fork copy-on-write, then performs
    spline fitting and phenology extraction. Defined at module level so it
    can be pickled by ProcessPoolExecutor.

    Parameters
    ----------
    coord_batch : list of (x, y)
        Pixel coordinates to process.
    smooth_x_axis : ndarray
        Regular daily grid used for spline evaluation.
    p : dict
        Parameters dict (see phenology docstring).
    t : ndarray
        Pre-computed time ordinal array for all time steps.

    Returns
    -------
    list of tuple
        Each element is (x, y, result, pheno) where result and pheno are None
        if the spline conditions are not satisfied for that pixel.
    """
    batch_results = []
    for x, y in coord_batch:
        values = _VAR_DATA[:, x, y]
        mask = values != -9999
        if p["qa_filter"]:
            mask = (_QA_DATA[:, x, y] == 0) & mask
        values = values[mask]
        time = t[mask]

        if len(time) < 2:
            batch_results.append((x, y, None, None))
            continue
        try:
            result = smooth_cubic_spline(
                time, values, p["spline_min_phase_length"], p["spline_min_relative_amplitude"],
                p["spline_min_phase_data"], smoothing_change=1e-12, max_iterations=1e5,
                smooth_x_axis=smooth_x_axis
            )
            if not result['conditions_satisfied']:
                batch_results.append((x, y, None, None))
                continue
            pheno = extract_phenology_metrics(
                result['x_pks'], result['y_pks'], result['x_trgs'], result['y_trgs'],
                time, result['smooth_x_axis'], result['smooth_y_data'],
                p["spline_subs_peak_win_size"], p["spline_subs_peak_ampl_frac"],
                p["spline_data_gap_size"], p["spline_data_gap_size_buffer"]
            )
            batch_results.append((x, y, result, pheno))
        except Exception as e:
            logging.warning("Pixel (%d, %d) failed: %s", x, y, e)
            batch_results.append((x, y, None, None))

    return batch_results


def main(args, log=False, parallel="lake", threads=1, batch_size=100, args_file=None) -> list[Path]:
    set_logging(log)
    args = parse_args(args)
    run_id = datetime.now(timezone.utc).isoformat().replace(":", "-")

    logging.info("Reading lake shapefile from: {}".format(args["shapefile"]))
    gdf = gpd.read_file(args["shapefile"])
    if args["start_index"] and args["end_index"]:
        gdf = gdf.iloc[args["start_index"]:args["end_index"]]
    elif args["lakes"]:
        gdf = gdf[gdf['id'].isin(args["lakes"])]

    files = glob.glob(f"{args['images']}/**/*.nc", recursive=True)
    files.sort()

    lakes = [row for i, (_, row) in enumerate(gdf.iterrows())]

    if args["outliers"]:
        logging.info("Extracting outliers {} from {}".format(args["variable"], args["out_folder"]))
        if parallel == "pixels" or threads == 1:
            for lake in lakes:
                e_files = extract_outliers(lake, args, threads=threads)
        else:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                e_files = executor.map(extract_outliers, lakes, itertools.repeat(args), itertools.repeat(files))
        if args["provenance"]:
            write_provenance(args["out_folder"], "outliers", args, args_file=args_file,
                                extra={"lakes": [int(lake["id"]) for lake in lakes], "threads": threads, "parallel": parallel}, run_id=run_id)

        logging.info("Fitting phenology to outliers {} from {}".format(args["variable"], args["out_folder"]))
        if parallel == "pixels" or threads == 1:
            for lake in lakes:
                phenology_outliers(lake, args, threads=threads)
        else:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                executor.map(phenology_outliers, lakes, itertools.repeat(args), itertools.repeat(files))
        if args["provenance"]:
            write_provenance(args["out_folder"], "outliers", args, args_file=args_file,
                                extra={"lakes": [int(lake["id"]) for lake in lakes], "threads": threads, "parallel": parallel}, run_id=run_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run outlier computation')
    parser.add_argument('--file', '-f', type=verify_arg_file, help='Name of the argument file in /args')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    parser.add_argument('--threads', '-t', help="Number of threads", default=1)
    parser.add_argument('--parallel', '-p', help="Run parallelisation on lakes or pixels", choices=['lakes', 'pixels'], default='lake')
    parser.add_argument('--batch-size', '-b', help="Number of pixels per batch for phenology I/O", type=int, default=100)
    args = parser.parse_args()
    with open(args.file) as f:
        file_args = json.load(f)
    main(file_args, log=args.logs, threads=int(args.threads), parallel=args.parallel,
            batch_size=args.batch_size, args_file=os.path.splitext(os.path.basename(args.file))[0])