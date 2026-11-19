"""Shared utilities for generating synthetic lake inputs."""

import json
import os
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
from shapely.geometry import box

# Resolve the repository root from the location of synthetic/common.py.
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def project_path(value):
    """Resolve relative paths against the repository root."""
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_defaults():
    """Load the JSON stored next to this module."""
    path = Path(__file__).with_name("default.json")
    with path.open(encoding="utf-8") as file:
        return json.load(file)

def load_scenario_metadata(lake_id, p):
    """Load one synthetic scenario's metadata from synth_lakes.xlsx.

    The sheet `synth_descriptions` is the authoritative source for the
    scenario ID, short name, name, tag and description.
    """
    metadata = p["metadata"]
    excel_path = project_path(metadata["excel_file"])
    sheet_name = metadata["sheet"]

    if not excel_path.exists():
        raise FileNotFoundError(
            f"Scenario metadata file not found: {excel_path}"
        )

    metadata_df = pd.read_excel(
        excel_path,
        sheet_name=sheet_name,
        engine="openpyxl",
    )

    required_columns = {
        "id",
        "short_name",
        "name",
        "tag",
        "Description",
    }
    missing_columns = required_columns.difference(metadata_df.columns)
    if missing_columns:
        raise ValueError(
            f"Missing metadata columns in '{sheet_name}': "
            f"{sorted(missing_columns)}"
        )

    metadata_df["id"] = pd.to_numeric(
        metadata_df["id"], errors="raise"
    ).astype("int64")

    duplicates = metadata_df["id"].duplicated()
    if duplicates.any():
        duplicate_ids = metadata_df.loc[duplicates, "id"].unique().tolist()
        raise ValueError(
            f"Duplicate scenario IDs in '{sheet_name}': {duplicate_ids}"
        )

    matches = metadata_df[metadata_df["id"] == int(lake_id)]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one metadata row for lake ID {lake_id}, "
            f"but found {len(matches)}."
        )

    row = matches.iloc[0]

    metadata_fields = {
        "id": int(row["id"]),
        "short_name": str(row["short_name"]).strip(),
        "name": str(row["name"]).strip(),
        "tag": str(row["tag"]).strip(),
        "description": str(row["Description"]).strip(),
    }

    missing_values = [
        field
        for field, value in metadata_fields.items()
        if value in ("", "nan", "None")
    ]
    if missing_values:
        raise ValueError(
            f"Missing metadata values for lake ID {lake_id}: "
            f"{missing_values}"
        )

    return metadata_fields


def create_grid_and_time(p):
    """Create daily UTC dates and geographic pixel centres."""
    grid = p["grid"]

    n_lat = int(grid["n_lat"])
    n_lon = int(grid["n_lon"])
    step = float(grid["pixel_size_degrees"])

    if n_lat < 1 or n_lon < 1 or step <= 0:
        raise ValueError("Grid dimensions and pixel size must be positive.")

    dates = pd.date_range(
        p["time"]["start"],
        p["time"]["end"],
        freq="D",
        tz="UTC",
    )
    if len(dates) < 3:
        raise ValueError("At least three daily observations are required.")

    south = float(grid["lat_min"])
    west = float(grid["lon_min"])
    north = south + n_lat * step
    east = west + n_lon * step

    if not (-90 <= south < north <= 90):
        raise ValueError("Invalid latitude bounds.")
    if not (-180 <= west < east <= 180):
        raise ValueError("Invalid longitude bounds.")

    lat = south + (np.arange(n_lat) + 0.5) * step
    lon = west + (np.arange(n_lon) + 0.5) * step
    geometry = box(west, south, east, north)

    return dates, lat, lon, geometry


def create_sinusoid(dates, p):
    """Create one sinusoidal cycle per calendar year."""
    baseline = p["baseline"]

    if baseline["type"] != "annual_sinusoid":
        raise ValueError("Only annual_sinusoid is implemented here.")

    amplitude = float(baseline["amplitude_ug_l"])
    if amplitude < 0:
        raise ValueError("Amplitude must be non-negative.")

    year_days = np.where(dates.is_leap_year, 366.0, 365.0)
    year_fraction = (dates.dayofyear.to_numpy() - 1) / year_days
    phase_fraction = float(baseline["phase_months"]) / 12.0

    return (
        float(baseline["mean_ug_l"])
        + amplitude * np.sin(
            2 * np.pi * (year_fraction - phase_fraction)
        )
    )

def create_lon_phase_gradient(dates, lat, lon, p):
    """Create an annual sinusoid delayed from west to east.

    The gradient is expressed in months per longitude pixel.
    All latitude rows share the same noise-free signal.
    """
    baseline = p["baseline"]

    if baseline["type"] != "annual_sinusoid":
        raise ValueError(
            "Only annual_sinusoid is implemented for longitude "
            "phase gradients."
        )

    mean = float(baseline["mean_ug_l"])
    amplitude = float(baseline["amplitude_ug_l"])
    base_phase_months = float(baseline["phase_months"])
    phase_gradient = float(
        p["spatial"]["phase_shift_months_per_pixel_lon"]
    )

    parameters = [mean, amplitude, base_phase_months, phase_gradient]
    if not np.isfinite(parameters).all():
        raise ValueError("Sinusoid parameters must be finite.")
    if amplitude < 0:
        raise ValueError("Amplitude must be non-negative.")
    if phase_gradient < 0:
        raise ValueError(
            "Phase gradient must be non-negative for a west-to-east delay."
        )

    year_days = np.where(dates.is_leap_year, 366.0, 365.0)
    year_fraction = (dates.dayofyear.to_numpy() - 1) / year_days

    lon_pixel = np.arange(len(lon), dtype=np.float64)
    phase_months = base_phase_months + phase_gradient * lon_pixel

    reference = mean + amplitude * np.sin(
        2 * np.pi
        * (
            year_fraction[:, None, None]
            - phase_months[None, None, :] / 12.0
        )
    )

    return np.broadcast_to(
        reference,
        (len(dates), len(lat), len(lon)),
    ).copy()

def create_lat_phase_gradient(dates, lat, lon, p):
    """Create an annual sinusoid delayed from south to north.

    The gradient is expressed in months per latitude pixel.
    All longitude columns share the same noise-free signal.
    """
    baseline = p["baseline"]

    if baseline["type"] != "annual_sinusoid":
        raise ValueError(
            "Only annual_sinusoid is implemented for latitude "
            "phase gradients."
        )

    mean = float(baseline["mean_ug_l"])
    amplitude = float(baseline["amplitude_ug_l"])
    base_phase_months = float(baseline["phase_months"])
    phase_gradient = float(
        p["spatial"]["phase_shift_months_per_pixel_lat"]
    )

    parameters = [mean, amplitude, base_phase_months, phase_gradient]
    if not np.isfinite(parameters).all():
        raise ValueError("Sinusoid parameters must be finite.")
    if amplitude < 0:
        raise ValueError("Amplitude must be non-negative.")
    if phase_gradient < 0:
        raise ValueError(
            "Phase gradient must be non-negative for a south-to-north delay."
        )

    year_days = np.where(dates.is_leap_year, 366.0, 365.0)
    year_fraction = (dates.dayofyear.to_numpy() - 1) / year_days

    lat_pixel = np.arange(len(lat), dtype=np.float64)
    phase_months = base_phase_months + phase_gradient * lat_pixel

    reference = mean + amplitude * np.sin(
        2 * np.pi
        * (
            year_fraction[:, None, None]
            - phase_months[None, :, None] / 12.0
        )
    )

    return np.broadcast_to(
        reference,
        (len(dates), len(lat), len(lon)),
    ).copy()

def create_lon_amplitude_gradient(dates, lat, lon, p):
    """Create an annual sinusoid with amplitude increasing west to east.

    The increase is expressed in ug/L per longitude pixel.
    Mean concentration and phase remain spatially constant.
    """
    baseline = p["baseline"]

    if baseline["type"] != "annual_sinusoid":
        raise ValueError(
            "Only annual_sinusoid is implemented for longitude "
            "amplitude gradients."
        )

    mean = float(baseline["mean_ug_l"])
    base_amplitude = float(baseline["amplitude_ug_l"])
    phase_months = float(baseline["phase_months"])
    amplitude_gradient = float(
        p["spatial"]["amplitude_increase_per_pixel_lon"]
    )

    parameters = [
        mean,
        base_amplitude,
        phase_months,
        amplitude_gradient,
    ]
    if not np.isfinite(parameters).all():
        raise ValueError("Sinusoid parameters must be finite.")
    if base_amplitude < 0:
        raise ValueError("Amplitude must be non-negative.")
    if amplitude_gradient < 0:
        raise ValueError(
            "Amplitude gradient must be non-negative "
            "for a west-to-east increase."
        )

    lon_pixel = np.arange(len(lon), dtype=np.float64)
    amplitude = base_amplitude + amplitude_gradient * lon_pixel

    if np.any(amplitude > mean):
        raise ValueError(
            "Amplitude exceeds mean concentration and would produce "
            "negative reference concentrations. Reduce the amplitude "
            "gradient or increase the mean."
        )

    year_days = np.where(dates.is_leap_year, 366.0, 365.0)
    year_fraction = (dates.dayofyear.to_numpy() - 1) / year_days

    reference = mean + amplitude[None, None, :] * np.sin(
        2 * np.pi
        * (
            year_fraction[:, None, None]
            - phase_months / 12.0
        )
    )

    return np.broadcast_to(
        reference,
        (len(dates), len(lat), len(lon)),
    ).copy()

def create_lat_amplitude_gradient(dates, lat, lon, p):
    """Create an annual sinusoid with amplitude increasing south to north.

    The increase is expressed in ug/L per latitude pixel.
    Mean concentration and phase remain spatially constant.
    """
    baseline = p["baseline"]

    if baseline["type"] != "annual_sinusoid":
        raise ValueError(
            "Only annual_sinusoid is implemented for latitude "
            "amplitude gradients."
        )

    mean = float(baseline["mean_ug_l"])
    base_amplitude = float(baseline["amplitude_ug_l"])
    phase_months = float(baseline["phase_months"])
    amplitude_gradient = float(
        p["spatial"]["amplitude_increase_per_pixel_lat"]
    )

    parameters = [
        mean,
        base_amplitude,
        phase_months,
        amplitude_gradient,
    ]
    if not np.isfinite(parameters).all():
        raise ValueError("Sinusoid parameters must be finite.")
    if base_amplitude < 0:
        raise ValueError("Amplitude must be non-negative.")
    if amplitude_gradient < 0:
        raise ValueError(
            "Amplitude gradient must be non-negative "
            "for a south-to-north increase."
        )

    lat_pixel = np.arange(len(lat), dtype=np.float64)
    amplitude = base_amplitude + amplitude_gradient * lat_pixel

    if np.any(amplitude > mean):
        raise ValueError(
            "Amplitude exceeds mean concentration and would produce "
            "negative reference concentrations. Reduce the amplitude "
            "gradient or increase the mean."
        )

    year_days = np.where(dates.is_leap_year, 366.0, 365.0)
    year_fraction = (dates.dayofyear.to_numpy() - 1) / year_days

    reference = mean + amplitude[None, :, None] * np.sin(
        2 * np.pi
        * (
            year_fraction[:, None, None]
            - phase_months / 12.0
        )
    )

    return np.broadcast_to(
        reference,
        (len(dates), len(lat), len(lon)),
    ).copy()

def create_daily_sections(n_dates):
    """Split a daily time series into three nearly equal sections."""
    if n_dates < 3:
        raise ValueError("At least three daily observations are required.")

    first_change = n_dates // 3
    second_change = (2 * n_dates) // 3

    return first_change, second_change


def broadcast_reference_to_grid(reference, lat, lon):
    """Broadcast a daily reference series across the full pixel grid."""
    reference = np.asarray(reference, dtype=np.float64)

    if reference.ndim != 1:
        raise ValueError("Reference must be a one-dimensional daily series.")

    return np.broadcast_to(
        reference[:, None, None],
        (len(reference), len(lat), len(lon)),
    ).copy()

def add_gaussian_noise(
    values,
    dates,
    p,
    std_ug_l=None,
    std_start_ug_l=None,
    std_end_ug_l=None,
):
    """Add independent Gaussian noise to a generated chl-a scenario.

    Supports either a constant standard deviation or a standard deviation
    that increases linearly over the daily time axis.
    """
    values = np.asarray(values, dtype=np.float64)
    n_dates = len(dates)

    if std_ug_l is not None:
        if std_start_ug_l is not None or std_end_ug_l is not None:
            raise ValueError(
                "Use either std_ug_l or std_start_ug_l and "
                "std_end_ug_l, not both."
            )

        std_start = float(std_ug_l)
        std_end = float(std_ug_l)
        noise_type = "gaussian_constant"

    elif std_start_ug_l is not None and std_end_ug_l is not None:
        if std_ug_l is not None:
            raise ValueError(
                "Use either std_ug_l or std_start_ug_l and "
                "std_end_ug_l, not both."
            )

        std_start = float(std_start_ug_l)
        std_end = float(std_end_ug_l)
        noise_type = "gaussian_linear_variance"

    else:
        raise ValueError(
            "Noise requires either std_ug_l or both "
            "std_start_ug_l and std_end_ug_l."
        )

    if not 0 <= std_start <= std_end:
        raise ValueError(
            "Noise standard deviations must satisfy "
            "0 <= std_start <= std_end."
        )

    if noise_type == "gaussian_constant":
        std = np.full(n_dates, std_start, dtype=np.float64)
    else:
        variance = np.linspace(
            std_start**2,
            std_end**2,
            n_dates,
        )
        std = np.sqrt(variance)

    rng = np.random.default_rng(p["noise"]["random_seed"])
    noise = rng.normal(size=values.shape) * std[:, None, None]

    return values + noise, noise_type

def update_shapefile(path, lake, geometry):
    """Add or replace one lake while keeping other lake entries."""
    row = gpd.GeoDataFrame(
        {
            "id": [int(lake["id"])],
            "short_name": [lake["short_name"]],
            "name": [lake["name"]],
            "Type": ["Synthetic"],
        },
        geometry=[geometry],
        crs="EPSG:4326",
    )

    if path.exists():
        existing = gpd.read_file(path)

        if "id" not in existing.columns or existing.crs is None:
            raise ValueError("Existing shapefile needs an id field and CRS.")

        existing = existing.to_crs(row.crs)
        existing = existing[
            existing["id"].astype(int) != int(lake["id"])
        ]
        row = gpd.GeoDataFrame(
            pd.concat([existing, row], ignore_index=True),
            geometry="geometry",
            crs=row.crs,
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    row.to_file(
        path,
        driver="ESRI Shapefile",
        index=False,
        encoding="UTF-8",
    )


def write_netcdf(path, lake, dates, lat, lon, values, p):
    """Write the extract-like input consumed by phenology.py."""
    variables = p["variables"]

    with netCDF4.Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("time", len(dates))
        nc.createDimension("lat", len(lat))
        nc.createDimension("lon", len(lon))

        time = nc.createVariable("time", "i8", ("time",))
        time.units = "seconds since 1970-01-01 00:00:00 UTC"
        time.calendar = "standard"
        time[:] = np.array(
            [int(date.timestamp()) for date in dates],
            dtype=np.int64,
        )

        latitude = nc.createVariable("lat", "f8", ("lat",))
        latitude.units = "degrees_north"
        latitude[:] = lat

        longitude = nc.createVariable("lon", "f8", ("lon",))
        longitude.units = "degrees_east"
        longitude[:] = lon

        data = nc.createVariable(
            variables["value"],
            "f4",
            ("time", "lat", "lon"),
            fill_value=float(variables["value_fill"]),
            zlib=True,
        )
        data.units = "ug L-1"
        data.long_name = "Synthetic chlorophyll-a concentration"
        data[:] = values.astype(np.float32)

        qa = nc.createVariable(
            variables["qa"],
            "i2",
            ("time", "lat", "lon"),
            fill_value=int(variables["qa_fill"]),
            zlib=True,
        )
        qa.flag_values = np.array([0, 1, 2, 3], dtype=np.int16)
        qa.flag_meanings = "good fair poor no_data"
        qa[:] = 0

        summary = nc.createVariable(
            "summary", "i4", ("lat", "lon")
        )
        summary[:] = np.full(
            (len(lat), len(lon)), len(dates), dtype=np.int32
        )
        summary.description = "Number of valid observations per pixel"

        nc.lake_id = int(lake["id"])
        nc.title = lake["name"]
        nc.variable = variables["value"]
        nc.qa = variables["qa"]
        nc.synthetic = "true"
        nc.generator_parameters = json.dumps(p)
        nc.scenario_metadata = json.dumps(lake)
        nc.phase_convention = (
            "Calendar-year sine; positive phase_months delays the "
            "cycle by phase_months/12 of a year."
        )

def save_results(lake, dates, lat, lon, geometry, values, reference, p):
    """Save a complete-data scenario and its reference outputs."""
    expected_shape = (len(dates), len(lat), len(lon))
    values = np.asarray(values, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)

    if values.shape != expected_shape:
        raise ValueError(f"Expected data shape {expected_shape}.")

    if reference.ndim == 1:
        if reference.shape != (len(dates),):
            raise ValueError(
                "1D reference must contain one value per date."
            )
        reference_kind = "time_series"
    elif reference.ndim == 3:
        if reference.shape != expected_shape:
            raise ValueError(
                "3D reference must match the generated data shape."
            )
        reference_kind = "spatial_grid"
    else:
        raise ValueError(
            "Reference must have shape (time,) or (time, lat, lon)."
        )

    if not np.isfinite(values).all() or not np.isfinite(reference).all():
        raise ValueError("This version supports finite, complete data only.")

    if np.any(values == float(p["variables"]["value_fill"])):
        raise ValueError("Generated values collide with the fill value.")

    if int(p["variables"]["qa_fill"]) in (0, 1, 2, 3):
        raise ValueError(
            "QA fill value must differ from valid QA flags."
        )

    output = p["output"]
    lake_id = int(lake["id"])

    nc_path = (
        project_path(output["out_folder"])
        / "extract"
        / p["variables"]["value"]
        / f"{lake_id}.nc"
    )
    reference_dir = (
        project_path(output["reference_folder"]) / str(lake_id)
    )
    shapefile_path = project_path(output["shapefile"])

    targets = [
        nc_path,
        reference_dir / "reference.csv",
        reference_dir / "reference.nc",
        reference_dir / "resolved_params.json",
        reference_dir / "preview.png",
    ]

    if not output.get("overwrite", False):
        existing = [str(path) for path in targets if path.exists()]
        if existing:
            raise FileExistsError(
                "Outputs already exist; no files were written:\n"
                + "\n".join(existing)
            )

    nc_path.parent.mkdir(parents=True, exist_ok=True)
    reference_dir.mkdir(parents=True, exist_ok=True)

    update_shapefile(shapefile_path, lake, geometry)
    write_netcdf(nc_path, lake, dates, lat, lon, values, p)

    if reference_kind == "time_series":
        pd.DataFrame({
            "date": dates,
            "reference_chla_ug_l": reference,
        }).to_csv(reference_dir / "reference.csv", index=False)
    else:
        reference_path = reference_dir / "reference.nc"

        with netCDF4.Dataset(
            reference_path, "w", format="NETCDF4"
        ) as nc:
            nc.createDimension("time", len(dates))
            nc.createDimension("lat", len(lat))
            nc.createDimension("lon", len(lon))

            time = nc.createVariable("time", "i8", ("time",))
            time.units = "seconds since 1970-01-01 00:00:00 UTC"
            time.calendar = "standard"
            time[:] = np.array(
                [int(date.timestamp()) for date in dates],
                dtype=np.int64,
            )

            latitude = nc.createVariable("lat", "f8", ("lat",))
            latitude.units = "degrees_north"
            latitude[:] = lat

            longitude = nc.createVariable("lon", "f8", ("lon",))
            longitude.units = "degrees_east"
            longitude[:] = lon

            data = nc.createVariable(
                "reference_chla_ug_l",
                "f4",
                ("time", "lat", "lon"),
                zlib=True,
            )
            data.units = "ug L-1"
            data.long_name = (
                "Noise-free synthetic chlorophyll-a reference"
            )
            data[:] = reference.astype(np.float32)

            nc.lake_id = int(lake["id"])
            nc.title = lake["name"]
            nc.reference_kind = "spatial_grid"

    with (reference_dir / "resolved_params.json").open(
        "w", encoding="utf-8"
    ) as file:
        json.dump(
            {"lake": lake, "parameters": p},
            file,
            indent=2,
        )

    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(
        dates,
        values[:, 0, 0],
        label="Pixel (0, 0)",
        color="tab:blue",
        alpha=0.6,
    )

    if reference_kind == "time_series":
        ax.plot(
            dates,
            reference,
            label="Noise-free reference",
            color="black",
            linewidth=1,
        )
    else:
        ax.plot(
            dates,
            reference[:, 0, 0],
            label="Reference, southwest pixel (0, 0)",
            color="black",
            linewidth=1,
        )
        ax.plot(
            dates,
            reference[:, 0, -1],
            label=f"Reference, southeast pixel (0, {len(lon) - 1})",
            color="tab:red",
            linewidth=1,
        )
        ax.plot(
            dates,
            reference[:, -1, 0],
            label=f"Reference, northwest pixel ({len(lat) - 1}, 0)",
            color="tab:green",
            linestyle="--",
            linewidth=1,
        )

    ax.set_xlabel("Date")
    ax.set_ylabel("Chlorophyll-a (µg/L)")
    ax.set_title(lake["name"])
    ax.legend()
    fig.tight_layout()
    fig.savefig(reference_dir / "preview.png", dpi=150)
    plt.close(fig)

    print(f"NetCDF: {nc_path}")
    print(f"Shapefile: {shapefile_path}")
    print(f"Reference outputs: {reference_dir}")
    print(f"Reference kind: {reference_kind}")
    print(f"Shape (time, lat, lon): {values.shape}")
    print(f"Range: {values.min():.3f} to {values.max():.3f} µg/L")