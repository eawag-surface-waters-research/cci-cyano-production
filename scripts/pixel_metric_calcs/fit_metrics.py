import os
from pathlib import Path
import logging
import netCDF4
import numpy as np
import pandas as pd
import datetime
import multiprocessing
from functools import partial
from shapely.geometry import Point
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from csaps import csaps
import warnings
from contextlib import contextmanager
import functions as f
from pixel_metric_calcs.base import PixelCalcBase

_GLOBALS = {}

METRICS = [
    "values_per_pixel",
    "r2",
    "MAD",
    "RMSE",
    "correlation",
]

def _calculate_pixel_metrics(coord,start_year,end_year):
    i, j = coord

    smoothing_all = _GLOBALS["smoothing_all"]
    values_all = _GLOBALS["values_all"]
    qa_all = _GLOBALS["qa_all"]
    t_all = _GLOBALS["t_all"]
    years_all = _GLOBALS["years_all"]

    smoothing = float(smoothing_all[i, j])

    values = values_all[:, i, j]
    qa_values = qa_all[:, i, j]

    mask = (
        (values != -9999)
        & (qa_values == 0)
    )

    values_m = values[mask]
    time_m = t_all[mask]
    years_m = years_all[mask]

    function_start, function_end = f.define_year_range(start= start_year, end = end_year, years = years_m)

    mask_sub = (years_m>= function_start) & (years_m <=function_end)
    valid = np.isfinite(y_true) & np.isfinite(y_pred)

    if valid.sum()<1:
            warnings.warn(f"Check data for lat, lon indices:{(i, j)}, perhaps smoothing parameter is nan or duplicates in time axis.")

    combined_mask = valid & mask_sub

    results = {
        "values_per_pixel": np.nan,
        "r2": np.nan,
        "MAD": np.nan,
        "RMSE": np.nan,
        "correlation": np.nan,
    }

    results["values_per_pixel"] = len(values_m)

    if len(values_m) < 2:
        return (i, j), results

    try:
        y_pred = csaps(
            time_m,
            values_m,
            time_m,
            smooth=smoothing,
        )
    except Exception as e:
        warnings.warn(f"csaps failed at {(i,j)}: {e}")
        return (i, j), results

    valid = (
        np.isfinite(values_m)
        & np.isfinite(y_pred)
    )

    if valid.sum() < 2:
        return (i, j), results

    y_true = values_m[combined_mask]
    y_pred = y_pred[combined_mask]

    results["r2"] = r2_score(
        y_true,
        y_pred,
    )

    results["MAD"] = np.median(
        np.abs(y_true - y_pred)
    )

    results["RMSE"] = np.sqrt(
        mean_squared_error(
            y_true,
            y_pred,
        )
    )

    try:
        results["correlation"], _ = pearsonr(
            y_true,
            y_pred,
        )
    except Exception:
        pass

    return (i, j), results

def _init_worker(p_path, e_path):
    """Initialise per-process globals for multiprocessing metric computation.

    Called once per worker process by multiprocessing.Pool. Loads the parameter
    and extract NetCDF datasets into module-level _GLOBALS so they are reused
    across all pixel-level calls within the same worker.

    Parameters
    ----------
    p_path : str
        Path to the phenology parameter NetCDF file.
    e_path : str
        Path to the extract NetCDF file containing satellite observations.
    """
    nc_p = netCDF4.Dataset(p_path)
    nc_e = netCDF4.Dataset(e_path)

    variable = getattr(nc_e, "variable")
    qa_name = getattr(nc_e, "qa")

    time_raw = nc_e.variables["time"][:]
    t_all = f.unix_to_datenum(time_raw)

    # Convert once, not per pixel
    time_dt = np.array(f.datenum_to_datetime(t_all))
    years_all = np.array([d.year for d in time_dt])

    smoothing_all = np.asarray(nc_p.variables["smoothing_parameter"][:])
    values_all = np.asarray(nc_e.variables[variable][:])
    qa_all = np.asarray(nc_e.variables[qa_name][:])

    lats = nc_e.variables["lat"][:]
    lons = nc_e.variables["lon"][:]

    _GLOBALS["nc_p"] = nc_p
    _GLOBALS["nc_e"] = nc_e
    _GLOBALS["variable"] = getattr(nc_e, "variable")
    _GLOBALS["qa"] = getattr(nc_e, "qa")
    _GLOBALS["t_all"] = t_all
    _GLOBALS["years_all"] = years_all
    _GLOBALS["smoothing_all"] = smoothing_all
    _GLOBALS["values_all"] = values_all
    _GLOBALS["qa_all"] = qa_all
    _GLOBALS["lats"] = lats
    _GLOBALS["lons"] = lons

class FitMetric(PixelCalcBase):
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def build_path(self,):
        """Return the cache path for one lake and one time window."""

        out_dir = os.path.join(
            self.out_folder,
            "calculated_values",
            "metrics",
            self.variable,
        )
        os.makedirs(out_dir, exist_ok=True)
        filename = f"ID{self.lakeID}_{self.start_year}_{self.end_year}.nc"
        self.save_fp = os.path.join(out_dir, filename)

    def read_input(self):
        self.p_path = Path(self.p_path).resolve()
        self.e_path = Path(self.e_path).resolve()
        if not self.p_path.is_file():
            raise FileNotFoundError(f"Invalid extract data_path {self.p_path}")

        if not self.e_path.is_file():
            raise FileNotFoundError(f"Invalid extract data_path {self.e_path}")
        return self.p_path, self.e_path

    def calculate(self):
        """Compute one spline-fit metric or observation count for a single pixel.

        Designed as a multiprocessing worker; reads all data from module-level
        _GLOBALS populated by _init_worker. For 'values_per_pixel' only the valid
        observation count is returned; for all other metrics the stored smoothing
        parameter is used to refit the spline and the requested statistic is computed.

        Parameters
        ----------
        coord : tuple of int
            (i, j) grid index pair identifying the pixel.
        start : int
            First year of the evaluation window (0 = earliest available year).
        end : int
            Last year of the evaluation window (9999 = latest available year).
        metrics_to_compute : list of str
            Single-element list naming the metric. One of:
            ['values_per_pixel'], ['r2'], ['MAD'], ['RMSE'], ['correlation'].

        Returns
        -------
        tuple
            ((i, j), metric_value) where metric_value is an int for
            'values_per_pixel' and a float (or np.nan) for all others.
        """
        with multiprocessing.Pool(
            initializer=_init_worker,
            initargs=(self.p_path, self.e_path,self.start_year,self.end_year),
            processes=min(
                10,
                os.cpu_count() or 4,
            ),
        ) as pool:

            self.metric_data = {
                metric: {}
                for metric in METRICS
            }

            for coord, metrics in pool.imap_unordered(
                _calculate_pixel_metrics,
                self.valid_coords,
            ):
                for metric_name, value in metrics.items():
                    self.metric_data[metric_name][coord] = value

        return self.metric_data

    def write_output(self):

        if not hasattr(self, "metric_data"):
            raise ValueError(
                "calculate() must be run first."
            )

        with netCDF4.Dataset(self.save_fp, "w") as ds:

            lat, lon = self._create_output(ds)

            for metric in METRICS:

                grid = np.full(
                    (len(lat), len(lon)),
                    -9999.0,
                    dtype=np.float32,
                )

                for (i, j), value in (
                    self.metric_data[metric]
                    .items()
                ):
                    grid[i, j] = value

                ds.variables[metric][:] = grid

    def read_cached_metric(self, metric_name):
        """Read one metric from a lake/window NetCDF cache."""
        with netCDF4.Dataset(self.save_fp, "r") as ds:
            if metric_name not in ds.variables:
                raise KeyError(
                    f"Metric {metric_name!r} is not present in {self.save_fp}"
                )

            ds.set_auto_mask(False)
            grid = np.asarray(ds.variables[metric_name][:, :])

        computed = np.argwhere(grid != -9999.0)

        return {(int(i), int(j)): grid[i, j] for i, j in computed}

    def _create_output(self, ds, block_size):

        with netCDF4.Dataset(self.e_path) as src:

            lat = np.asarray(
                src.variables["lat"][:]
            )

            lon = np.asarray(
                src.variables["lon"][:]
            )

        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))

        ds.createVariable(
            "lat",
            "f8",
            ("lat",),
        )[:] = lat

        ds.createVariable(
            "lon",
            "f8",
            ("lon",),
        )[:] = lon

        for metric in METRICS:

            ds.createVariable(
                metric,
                "f4",
                ("lat", "lon"),
                fill_value=-9999.0,
                zlib=True,
                complevel=4,
                chunksizes=(
                    min(block_size, len(lat)),
                    min(block_size, len(lon)),
                ),
            )

        ds.lake_id = str(self.lakeID)
        ds.version = str(self.version)
        ds.variable = str(self.variable)
        ds.metric_name = "fit_metrics"
        ds.start_year = self.start_year
        ds.end_year = self.end_year

        return lat, lon
    
    def calculate_and_write_chunked(
        self,
        block_size=1000,
    ):
        """
        Calculate all fit metrics and write directly to disk in pixel chunks.

        Parameters
        ----------
        block_size : int, optional
            Number of pixels processed per chunk.
        """

        with netCDF4.Dataset(self.save_fp, "w") as ds:

            self._create_output(ds,block_size = block_size)

            with multiprocessing.Pool(
                initializer=_init_worker,
                initargs=(self.p_path, self.e_path, self.start_year,self.end_year),
                processes=min(
                    10,
                    os.cpu_count() or 4,
                ),
            ) as pool:

                for start in range(
                    0,
                    len(self.valid_coords),
                    block_size,
                ):

                    stop = min(
                        start + block_size,
                        len(self.valid_coords),
                    )

                    coords = self.valid_coords[start:stop]

                    results = pool.map(
                        _calculate_pixel_metrics,
                        coords,
                    )

                    for coord, metrics in results:

                        i, j = coord

                        for metric_name, value in metrics.items():

                            ds.variables[
                                metric_name
                            ][i, j] = value

        return self.save_fp