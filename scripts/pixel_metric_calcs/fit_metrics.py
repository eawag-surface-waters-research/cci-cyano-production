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
from scipy.stats import pearsonr, gaussian_kde
import warnings
from contextlib import contextmanager
import functions as f
from pixel_metric_calcs.base import PixelCalcBase

_GLOBALS = {}


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

    def build_path(self, metric_name, start=0, end=9999):
        """Return the cache path for one lake and one time window."""
        start_label = self.start_year if start <= self.start_year else start
        end_label = self.end_year if end >= self.end_year else end

        out_dir = os.path.join(
            self.data_folder,
            "calculated_values",
            "metrics",
            self.variable,
        )
        os.makedirs(out_dir, exist_ok=True)
        filename = f"ID{self.lakeID}_{start_label}_{end_label}.nc"
        self.save_fp = os.path.join(out_dir, filename)

    def read_input(self):
        pass

    def calculate(coord, start=0, end=9999, metrics_to_compute= None):
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
        if metrics_to_compute is None:
            metrics_to_compute = ["values_per_pixel", "r2", "MAD", "RMSE", "correlation"]

        if metrics_to_compute==["values_per_pixel"]:
            i,j = coord
            values_all = _GLOBALS["values_all"]
            qa_all = _GLOBALS["qa_all"]
            years_all = _GLOBALS["years_all"]

            values = values_all[:, i, j]
            qa_values = qa_all[:, i, j]

            mask = (values != -9999) & (qa_values == 0)
            values_m = values[mask]

            years_m = years_all[mask]

            function_start, function_end = f.define_year_range(start= start, end= end, years= years_m)

            mask_sub = (years_m>= function_start) & (years_m <=function_end)

            final_values = values_m[mask_sub]

            return (i,j), len(final_values)
        else:
            i,j = coord

            smoothing_all = _GLOBALS["smoothing_all"]
            values_all = _GLOBALS["values_all"]
            qa_all = _GLOBALS["qa_all"]
            t_all = _GLOBALS["t_all"]
            years_all = _GLOBALS["years_all"]


            smoothing = float(smoothing_all[i, j])
            values = values_all[:, i, j]
            qa_values = qa_all[:, i, j]

            mask = (values != -9999) & (qa_values== 0)
            values_m = values[mask]
            time_m = t_all[mask]
            years_m = years_all[mask]
            
            if len(values_m)>1:

                function_start, function_end = f.define_year_range(start= start, end = end, years = years_m)

                y_pred =csaps(time_m, values_m, time_m, smooth=smoothing)
                y_true = values_m

                mask_sub = (years_m>= function_start) & (years_m <=function_end)
                valid = np.isfinite(y_true) & np.isfinite(y_pred)

                if valid.sum()<1:
                        warnings.warn(f"Check data for lat, lon indices:{(i, j)}, perhaps smoothing parameter is nan or duplicates in time axis.")

                combined_mask = valid & mask_sub

                if combined_mask.sum() > 1:
                    if metrics_to_compute == ["r2"]:
                        metric = r2_score(y_true[combined_mask], y_pred[combined_mask])
                    elif metrics_to_compute == ["MAD"]:
                        metric = np.median(np.abs(y_true[combined_mask]-y_pred[combined_mask]))
                    elif metrics_to_compute == ["RMSE"]:
                        metric = np.sqrt(mean_squared_error(y_true[combined_mask], y_pred[combined_mask]))
                    elif metrics_to_compute == ["correlation"]:
                        metric, _ = pearsonr(y_true[combined_mask], y_pred[combined_mask])
                    else:
                        raise ValueError("please enter a valid metric")

                else:
                    warnings.warn(f"Not enough valid data in selected date range for indices {(i,j)}")
                    metric = np.nan
            else:
                metric = np.nan
            return (i,j), metric

    def write_output(self, file_path, col_name, data):
        """Append one metric to a lake/window NetCDF cache."""
        with netCDF4.Dataset(self.e_path) as src:
            lat = np.asarray(src.variables["lat"][:])
            lon = np.asarray(src.variables["lon"][:])

        mode = "a" if os.path.isfile(file_path) else "w"

        with netCDF4.Dataset(file_path, mode) as ds:
            if mode == "w":
                ds.createDimension("lat", len(lat))
                ds.createDimension("lon", len(lon))

                ds.createVariable("lat", "f8", ("lat",))[:] = lat
                ds.createVariable("lon", "f8", ("lon",))[:] = lon

                ds.lake_id = str(self.lakeID)
                ds.version = str(self.version)
                ds.variable = str(self.variable)

            grid = np.full(
                (len(lat), len(lon)),
                -9999.0,
                dtype=np.float32,
            )

            for (i, j), value in data.items():
                grid[i, j] = value

            if col_name in ds.variables:
                metric_var = ds.variables[col_name]
            else:
                metric_var = ds.createVariable(
                    col_name,
                    "f4",
                    ("lat", "lon"),
                    fill_value=-9999.0,
                    zlib=True,
                    complevel=4,
                )

            metric_var[:, :] = grid

    def read_cached_metric(self, file_path, col_name):
        """Read one metric from a lake/window NetCDF cache."""
        with netCDF4.Dataset(file_path, "r") as ds:
            if col_name not in ds.variables:
                raise KeyError(
                    f"Metric {col_name!r} is not present in {file_path}"
                )

            ds.set_auto_mask(False)
            grid = np.asarray(ds.variables[col_name][:, :])

        computed = np.argwhere(grid != -9999.0)

        return {(int(i), int(j)): grid[i, j] for i, j in computed}

    def compute_and_cache_metric(self, metric_name, col_name, compute_fn,
                                 start=0, end=9999):
            """ Copied over from visualization. Check for utility!
            Compute or load one metric for a lake and time window."""
            is_netcdf = self.save_format == "netcdf"
            dir_path, file_path = self.build_metric_path(metric_name, start=start, end=end)
            os.makedirs(dir_path, exist_ok=True)
    
            data = None
    
            if os.path.isfile(file_path):
                if is_netcdf:
                    with netCDF4.Dataset(file_path, "r") as ds:
                        metric_exists = col_name in ds.variables
    
                    if metric_exists:
                        data = self._read_metric_netcdf(file_path, col_name)
                else:
                    df = pd.read_csv(file_path)
                    data = dict(zip(zip(df["i"], df["j"]), df[col_name]))
    
                if data is not None:
                    missing = [coord for coord in self.valid_coords if coord not in data]
    
                    if not missing:
                        return data
    
                    warnings.warn(
                        f"Cached {metric_name} is missing {len(missing)} valid pixels; "
                        "recomputing."
                    )
            else:
                warnings.warn(
                    f"{metric_name} needs to be calculated. "
                    "Depending on lake size, this may take a while."
                )
    
            workers = partial(
                compute_fn,
                start=start,
                end=end,
                metrics_to_compute=[metric_name],
            )
    
            with multiprocessing.Pool(
                initializer=_init_worker,
                initargs=(self.p_path, self.e_path),
                processes=3,
            ) as pool:
                result = pool.map(workers, self.valid_coords)
    
            data = dict(result)
    
            if is_netcdf:
                self._write_metric_netcdf(file_path, col_name, data)
            else:
                metric_df = pd.DataFrame(
                    [
                        (i, j, value)
                        for (i, j), value in data.items()
                    ],
                    columns=["i", "j", col_name],
                )
                metric_df.to_csv(file_path, index=False)
    
            return data

