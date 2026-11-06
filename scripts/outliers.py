"""Select candidate dates and create cumulative outlier-removal inputs."""

import logging
import shutil
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr


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
            dataset.setncattr(
                "outlier_source_file",
                input_file.name,
            )
            dataset.setncattr(
                "outlier_removed_count",
                number_removed,
            )
            dataset.setncattr(
                "outlier_removed_dates",
                ", ".join(
                    date.isoformat()
                    for date in candidate_timestamps[:number_removed]
                ),
            )

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


def outliers(lake, p, threads=1, batch_size=100) -> list[Path]:
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

    logging.info(
        "Baseline remains unchanged: %s",
        input_file,
    )

    logging.info(
        "Calculate matching phenology fits separately in %s. "
        "Existing fits are not updated by input preparation.",
        phenology_outliers_dir,
    )

    return created_files