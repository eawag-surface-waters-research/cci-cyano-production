import os
from pathlib import Path
import logging
import netCDF4
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
import warnings

from base import PixelCalcBase
import functions as f

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


def _init_kde_worker(p_path, var_names):
    """Initialise per-process globals for KDE pixel extraction.

    Called once per worker process by multiprocessing.Pool. Preloads all
    required phenology arrays as numpy arrays so per-pixel work is pure
    in-memory indexing with no NetCDF I/O in the hot path.

    Parameters
    ----------
    p_path : str
        Path to the phenology NetCDF file.
    var_names : list of str
        NetCDF variable names to preload (determined by assemble_kde_data).
    """
    with netCDF4.Dataset(p_path) as nc:
        for var in var_names:
            v = nc.variables[var]
            nlat, nlon, nrec = v.shape
            arr = np.empty((nlat, nlon, nrec), dtype=v.dtype)
            # read pixel by pixel (v[i, j, :]) rather than a bulk v[:] read - netCDF4
            # 1.7.4 silently misattributes data between pixels when read this way for
            # files with an unlimited 'record' dimension, verified against the trusted
            # per-pixel access pattern used elsewhere in this class (e.g. _load_pixel_data)
            for i in range(nlat):
                for j in range(nlon):
                    arr[i, j, :] = v[i, j, :]
            _GLOBALS[f"kde_{var}"] = arr


def _init_bloom_probability_worker(kde):
    _GLOBALS["bloom_kde"] = kde


def _evaluate_bloom_probability_rows(args):
    """Evaluate KDE density for a chunk of y-grid rows."""
    y_rows, xi = args
    kde = _GLOBALS["bloom_kde"]

    Xi, Yi = np.meshgrid(xi, y_rows)
    density = kde(
        np.vstack([Xi.ravel(), Yi.ravel()])
    ).reshape(Xi.shape)

    return density


class SpatialAgg(PixelCalcBase):
    # self.__init__(metric_name = "background_median")
    def build_path(self):
        out_dir = os.path.join(
                self.out_folder,
                "calculated_values",
                "spatial_aggregation_values",
                self.variable,
            )
        os.makedirs(out_dir, exist_ok=True)
        self.save_fp = os.path.join(out_dir, f"ID{self.lakeID}_background_agg.nc")

    def read_input(self):
        self.data_path = os.path.join(self.out_folder,'extract',self.variable,f"{self.lakeID}.nc")
    
    def calculate(self):
        with netCDF4.Dataset(self.data_path) as nc:
            lat = np.asarray(nc.variables["lat"][:])
            lon = np.asarray(nc.variables["lon"][:])
            t_all = f.unix_to_datenum(nc.variables["time"][:])

            variable_name = getattr(nc, "variable")
            data_var = nc.variables[variable_name]
            qa_var = nc.variables[getattr(nc, "qa")]

            ntime = len(nc.dimensions["time"])
            nlat = len(nc.dimensions["lat"])
            nlon = len(nc.dimensions["lon"])

            coords = np.asarray(self.valid_coords, dtype=int)

            # Remove border cells once
            interior_mask = (
            (coords[:, 0] >= 1) & (coords[:, 0] < nlat - 1) &
            (coords[:, 1] >= 1) & (coords[:, 1] < nlon - 1)
            )
            coords = coords[interior_mask]

            if coords.size == 0:
                # self._write_aggregation_netcdf(t_all,
                #     np.empty(0, dtype=int), np.empty(0, dtype=int),
                #     np.empty(0), np.empty(0),
                #     np.empty((ntime, 0), dtype=np.float32),
                # )
                return

            i_idx = coords[:, 0]
            j_idx = coords[:, 1]

            # indices for median_grid, which is smaller by 1 border cell each side
            ii = i_idx - 1
            jj = j_idx - 1

            lat_vals = lat[i_idx]
            lon_vals = lon[j_idx]

            n_pixels = len(coords)
            # netcdf path accumulates a (ntime, n_pixels) array and writes it in one
            # bulk call: writing timestep-by-timestep into chunks that span the full
            # time dimension would force a decompress/recompress of every chunk on
            # every iteration.
            values = np.full((ntime, n_pixels), np.nan, dtype=np.float32)

            for n in range(ntime):
                data_n = np.asarray(data_var[n], dtype=np.float32)
                qa_n = np.asarray(qa_var[n])

                # shape: (nlat-2, nlon-2, 3, 3)
                data_windows = np.lib.stride_tricks.sliding_window_view(data_n, (3, 3))
                qa_windows = np.lib.stride_tricks.sliding_window_view(qa_n, (3, 3))

                valid_mask = (data_windows != -9999) & (qa_windows == 0)

                masked = data_windows.astype(np.float32, copy=True)
                masked[~valid_mask] = np.nan

                # shape: (nlat-2, nlon-2)
                median_grid = np.nanmedian(masked, axis=(-2, -1))

                ma_values = median_grid[ii, jj]

                values[n, :] = ma_values

    def write_output(self, t_all, i_idx, j_idx, lat_vals, lon_vals, values):
        """Write spatial_aggregation() results to a compressed (time, pixel) NetCDF cache.

        lat/lon are stored once per pixel rather than once per row (the CSV's main
        source of bloat), and MA_value is chunked as (ntime, 1) so a per-pixel read
        - the only access pattern plot_background_pts uses - pulls exactly one
        contiguous chunk instead of scanning the whole file.
        """
        n_pixels = len(i_idx)
        with netCDF4.Dataset(self.save_fp, "w") as ds:
            ds.createDimension("time", len(t_all))
            ds.createDimension("pixel", n_pixels)

            ds.createVariable("time", "f8", ("time",))[:] = t_all
            ds.createVariable("pixel_i", "i4", ("pixel",))[:] = i_idx
            ds.createVariable("pixel_j", "i4", ("pixel",))[:] = j_idx
            ds.createVariable("lat", "f8", ("pixel",))[:] = lat_vals
            ds.createVariable("lon", "f8", ("pixel",))[:] = lon_vals

            chunksizes = (len(t_all), 1) if n_pixels > 0 else None
            ma_var = ds.createVariable(
                "MA_value", "f4", ("time", "pixel"),
                fill_value=np.nan, zlib=True, complevel=4,
                chunksizes=chunksizes,
            )
            ma_var[:, :] = values

    def read_cached_metric(self):
        """Open a cached aggregation NetCDF file and rebuild the pixel lookup index."""
        self.aggregation_ds = netCDF4.Dataset(self.save_fp, "r")
        pixel_i = np.asarray(self.aggregation_ds.variables["pixel_i"][:]).tolist()
        pixel_j = np.asarray(self.aggregation_ds.variables["pixel_j"][:]).tolist()
        self._aggregation_pixel_index = dict(zip(zip(pixel_i, pixel_j), range(len(pixel_i))))


class KDE(PixelCalcBase):

    def build_path(self):
        """Return the output directory and file path for the cached KDE events CSV.

        Returns
        -------
        base : str
            Directory path where the CSV will be written.
        file_path : str
            Full path to the CSV file.
        """

        start_label = self.start_year
        end_label = self.end_year
        # lake_name = f.sanitize_filename(self.ID_to_name(int(self.lakeID)).replace(" ", ""))
        out_dir = os.path.join(
            self.data_folder,
            "calculated_values",
            "kde_data",
            self.variable,
        )
        filename = f"ID{self.lakeID}_{start_label}_{end_label}.nc"
        os.makedirs(out_dir, exist_ok=True)
        self.save_fp = os.path.join(out_dir, filename)

    def read_input(self):
        pass

    def calculate(self):
        pass

    def write_output(self):
        """Write cached KDE event data to NetCDF."""
        columns = ["primary", "qa_column", "secondary"]

        with netCDF4.Dataset(self.save_fp, "w") as ds:
            ds.createDimension("event", len(self.kde_df))

            for column in columns:
                values = self.kde_df[column].to_numpy()

                if column == "qa_column":
                    variable = ds.createVariable(column, "f4", ("event",))
                    variable[:] = values
                else:
                    variable = ds.createVariable(column, "f8", ("event",))
                    variable[:] = values

            ds.lake_id = str(self.lakeID)
            ds.version = str(self.version)
            ds.variable = str(self.variable)

    def read_cached_metric(self):
        """Read cached KDE event data from NetCDF."""
        with netCDF4.Dataset(self.save_fp, "r") as ds:
            required = {"primary", "qa_column", "secondary"}
            missing = required.difference(ds.variables)

            if missing:
                raise ValueError(
                    f"KDE cache {self.save_fp} is missing variables: {sorted(missing)}"
                )

            return pd.DataFrame({
                "primary": np.asarray(ds.variables["primary"][:]),
                "qa_column": np.asarray(ds.variables["qa_column"][:]),
                "secondary": np.asarray(ds.variables["secondary"][:]),
            })

        
    def _extract_pixel_kde_events(self, nc, i, j,
                                   primary_vars=None, secondary_vars=None,
                                   qa_var="pks", arrays=None):
        """Extract bracketing and peak events for a single pixel.

        Parameters
        ----------
        nc : netCDF4.Dataset or None
            Open phenology dataset. Ignored when arrays is provided.
        i, j : int
            Pixel row and column indices.
        primary_vars : list of str, optional
            NetCDF variable names whose events mark the start of a bloom bracket
            (e.g. green-up). Defaults to ['green_up_advanced'].
        secondary_vars : list of str, optional
            NetCDF variable names whose events mark the end of a bloom bracket
            (e.g. green-down onset). Defaults to ['green_down_onset'].
        qa_var : str, optional
            Short name for the peak QA variable (passed through parse_qa_var_from_str).
            Defaults to 'pks'.
        arrays : dict of str -> np.ndarray, optional
            Preloaded full arrays keyed by NetCDF variable name. When provided,
            pixel slices are taken from these in-memory arrays instead of reading
            from nc, which is the fast path used by _kde_pixel_worker.

        Returns
        -------
        pandas.DataFrame
            Columns: year.DOY, i, j, green_up_advanced, peaks, green_down_onset.
            Empty DataFrame if the pixel has no events.
        """
        if primary_vars is None:
            primary_vars = ["green_up_advanced"]
        if secondary_vars is None:
            secondary_vars = ["green_down_onset"]

        def _get(var_name):
            if arrays is not None:
                return arrays[var_name][i, j, :]
            return nc.variables[var_name][i, j, :]

        frames = []

        for var in primary_vars:
            var_x = f.coerce_varname_to_var_x(var)
            raw = f.remove_nan(_get(var_x))
            if len(raw) == 0:
                continue
            dt = pd.to_datetime(raw, unit="s", utc=True)
            frames.append(pd.DataFrame({
                "year.DOY": dt.year + dt.day_of_year / 1000,
                "i": i, "j": j,
                "primary": True,
                "qa_column": np.nan,
                "secondary": False,
            }))

        qa_var_parsed = f.parse_qa_var_from_str(qa_var)
        if qa_var_parsed is not None:
            qa_x = f.coerce_varname_to_var_x(qa_var_parsed)
            qa_x_raw = np.array(_get(qa_x))
            pk_mask = ~np.isnan(qa_x_raw)
            if pk_mask.any():
                pks_dt = pd.to_datetime(qa_x_raw[pk_mask], unit="s", utc=True)
                frames.append(pd.DataFrame({
                    "year.DOY": pks_dt.year + pks_dt.day_of_year / 1000,
                    "i": i, "j": j,
                    "primary": False,
                    "qa_column": np.array(_get(qa_var_parsed))[pk_mask].astype(int),
                    "secondary": False,
                }))

        for var in secondary_vars:
            var_x = f.coerce_varname_to_var_x(var)
            raw = f.remove_nan(_get(var_x))
            if len(raw) == 0:
                continue
            dt = pd.to_datetime(raw, unit="s", utc=True)
            frames.append(pd.DataFrame({
                "year.DOY": dt.year + dt.day_of_year / 1000,
                "i": i, "j": j,
                "primary": False,
                "qa_column": np.nan,
                "secondary": True,
            }))

        if not frames:
            return pd.DataFrame(columns=["year.DOY", "i", "j", "primary", "qa_column", "secondary"])
        return pd.concat(frames, ignore_index=True)

    @staticmethod
    def _kde_pixel_worker(coord, primary_vars, secondary_vars, qa_var):
        """Multiprocessing worker for a single pixel's KDE event extraction.

        Uses preloaded numpy arrays from _GLOBALS (populated by _init_kde_worker)
        for pure in-memory indexing — no NetCDF I/O in the hot path.

        Parameters
        ----------
        coord : tuple of int
            (i, j) grid index pair.
        primary_vars, secondary_vars : list of str
            Passed through to the extraction logic.
        qa_var : str
            Short name for the peak QA variable.

        Returns
        -------
        pandas.DataFrame
            Same schema as _extract_pixel_kde_events; empty if no events found.
        """
        i, j = coord
        arrays = {k[4:]: v for k, v in _GLOBALS.items() if k.startswith("kde_")}
        return KDE._extract_pixel_kde_events(
            None, None, i, j, primary_vars, secondary_vars, qa_var, arrays=arrays
        )


    # MOVE TO PixelCalc
    def assemble_kde_data(self, primary_vars=None, secondary_vars=None, qa_var="pks"):
        """Collect bracketing and peak events lake-wide into a DataFrame.

        Iterates all valid pixels inside the 1 km-inset lake boundary and gathers
        events for each variable type. Each row represents one event occurrence.

        Parameters
        ----------
        primary_vars : list of str, optional
            Variables marking the start of a bloom bracket. Defaults to ['green_up_advanced'].
        secondary_vars : list of str, optional
            Variables marking the end of a bloom bracket. Defaults to ['green_down_onset'].
        qa_var : str, optional
            Short name for the peak QA variable. Defaults to 'pks'.

        Returns
        -------
        pandas.DataFrame
            Index named 'year.DOY' — a float of the form year + DOY/1000
            (e.g. 2005.150 = year 2005, day-of-year 150).
            Columns: i, j, primary (bool), qa_column (float), secondary (bool).
            Row count equals the total number of events across all variables and pixels.
        """

        print(f"assemble kde started at: {datetime.datetime.now()}")
        g = self._load_extracted_globals()
        lats = g["lat"]
        lons = g["lon"]

        inset_coords = [
            (i, j) for (i, j) in self.valid_coords
            if self.prepped_geom.contains(Point(lons[j], lats[i]))
        ]

        # Compute the full set of NetCDF variable names needed so the initializer
        # can preload them as numpy arrays — eliminates per-pixel disk reads.
        var_names_set = set()
        _pv = primary_vars if primary_vars is not None else ["green_up_advanced"]
        _sv = secondary_vars if secondary_vars is not None else ["green_down_onset"]
        for var in _pv:
            var_names_set.add(f.coerce_varname_to_var_x(var))
        qa_var_parsed = f.parse_qa_var_from_str(qa_var)
        if qa_var_parsed is not None:
            var_names_set.add(f.coerce_varname_to_var_x(qa_var_parsed))
            var_names_set.add(qa_var_parsed)
        for var in _sv:
            var_names_set.add(f.coerce_varname_to_var_x(var))
        var_names = list(var_names_set)

        worker = partial(
            PhenologyVisualization._kde_pixel_worker,
            primary_vars=primary_vars, secondary_vars=secondary_vars, qa_var=qa_var,
        )
        n_workers = min(10, os.cpu_count() or 4)
        chunksize = max(1, len(inset_coords) // (n_workers * 4))
        with multiprocessing.Pool(
            initializer=_init_kde_worker, initargs=(self.p_path, var_names), processes=n_workers
        ) as pool:
            results = list(pool.imap_unordered(worker, inset_coords, chunksize=chunksize))

        frames = [df for df in results if not df.empty]
        if not frames:
            return pd.DataFrame(columns=["i", "j", "primary", "qa_column", "secondary"])
        print(f"assemble kde finished at: {datetime.datetime.now()}")
        return pd.concat(frames, ignore_index=True).set_index("year.DOY")


    # MOVE TO PixelCalc
    def _fit_bloom_kde(self, qa_value = None, start_year = 0, end_year = 9999):
        """
        Load/filter cached bloom events and fit a 2D gaussian_kde on
        (green-up advance DOY, green-down onset DOY).

        Returns
        -------
        tuple or None
            (kde, start, end, qa_filtered_set), or None if there isn't
            enough data to fit a KDE.
        """
        dir_path, file_path = self.build_kde_path()
        # the cached CSV doesn't retain per-pixel (i, j) identity (only pooled
        # primary/qa_column/secondary event columns), so unlike compute_and_cache_metric
        # we can't detect staleness by checking pixel coverage - fall back to comparing
        # against the source NetCDFs' modification time instead
        source_mtime = max(os.path.getmtime(self.p_path), os.path.getmtime(self.e_path))
        cache_is_stale = os.path.isfile(file_path) and os.path.getmtime(file_path) < source_mtime

        if os.path.isfile(file_path) and not cache_is_stale:
            if self.save_format == "netcdf":
                compressed_df = self._read_kde_nc(file_path)
            else:
                compressed_df = pd.read_csv(file_path)
        else:
            if cache_is_stale:
                warnings.warn(f"Cached KDE events for lake ID {self.lakeID} predate the source data; recomputing.")
            else:
                warnings.warn("KDE events need to be calculated. Depending on the lake size this may take a while.")
            os.makedirs(dir_path, exist_ok=True)
            df = self.assemble_kde_data()
            compressed_df = self.prep_kde_data(df)
            if self.save_format == "netcdf":
                self._write_kde_netcdf(file_path, compressed_df)
            else:
                compressed_df.to_csv(file_path, index=False)

        qa_filtered_set = None
        if qa_value is not None:
            if type(qa_value) != set:
                warnings.warn("qa_value needs to be a set")
            else:
                compressed_df = compressed_df[compressed_df['qa_column'].isin(qa_value)]
                qa_filtered_set = qa_value

        if len(compressed_df) < 2:
            warnings.warn("Not enough data to plot kde")
            return None

        years_all = np.unique(list(compressed_df["primary"].astype(int) ) + list(compressed_df["secondary"].astype(int) ))
        start, end = f.define_year_range(start_year, end_year, years_all)
        plot_df = f.sort_by_year(compressed_df, start_year=start, end_year=end)

        x = np.round((plot_df["primary"].values % 1) * 1000).astype(int)
        y = np.round((plot_df["secondary"].values % 1) * 1000).astype(int)
        y[y< x] += 365

        try:
            kde = gaussian_kde(np.vstack([x, y]))
        except np.linalg.LinAlgError:
            # too few / too degenerate (collinear or duplicate) points for a 2D KDE -
            # e.g. can happen with a restrictive qa_value filter that leaves very few events
            warnings.warn(f"Not enough distinct {self.variable} events to plot KDE for lake ID {self.lakeID}.")
            return None

        return kde, start, end, qa_filtered_set


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


class BloomProb(PixelCalcBase):

    def build_path(self,
        qa_value=None,
        start_year=0,
        end_year=9999,
        interval=21,
        resolution=1,
        x_max=400,
        y_max=730,
    ):
        """Return the output directory and file path for the cached bloom probability events netcdf.

        Returns
        -------
        base : str
            Directory path where the netcdf will be written.
        file_path : str
            Full path to the netcdf file.
        """
        ext = "nc" if self.save_format == "netcdf" else "csv"

        start_label = self.start_year if start_year <= self.start_year else start_year
        end_label = self.end_year if end_year >= self.end_year else end_year

        if qa_value is None:
            qa_label = "all"
        else:
            qa_label = "qa" + "".join(map(str, sorted(qa_value)))

        base = os.path.join(
            self.data_folder,
            "calculated_values",
            "bloom_prob",
            self.variable,
        )

        filename = (
            f"ID{self.lakeID}_{start_label}_{end_label}_"
            f"{qa_label}_i{interval}_r{resolution}_x{x_max}_y{y_max}.{ext}"
        )

        return base, os.path.join(base, filename)

    def read_input(self):
        pass

    def calculate(self, qa_value = None, start_year = 0, end_year = 9999,
                                     interval = 21, x_max = 400, y_max = 730, resolution = 1):
        """
        Probability of a `interval`-day bloom window under the fitted KDE, for
        a window starting at every `resolution`-day offset spanning
        [-interval, x_max+interval] x [-interval, y_max+interval] (green-up
        advance DOY x green-down onset DOY). Windows overlap (a 1-day step
        with a 21-day window means neighbors share 20 days), so `probability`
        values do NOT sum to 1 - each is its own "P(bloom falls in *this*
        21-day window)", not a share of a disjoint partition.

        Computed as a Riemann-sum approximation: the KDE density is evaluated
        once on the full `resolution`-day grid (a single vectorized call),
        then every window's integral is read off a 2D summed-area table
        (prefix sums), rather than calling the exact but expensive
        kde.integrate_box() per window - that brute-force approach means
        ~316k individual calls at interval=21/resolution=1, on the order of
        an hour or more; this is seconds, at the cost of a small
        discretization error from the finite `resolution` (~0.1% in testing
        at resolution=1).

        Returns
        -------
        pandas.DataFrame or None
            One row per window with x_low/x_high/y_low/y_high/probability
            columns. None if there isn't enough data to fit a KDE.
        """

        start_label = (
        self.start_year if start_year <= self.start_year else start_year
        )
        end_label = (
        self.end_year if end_year >= self.end_year else end_year
        )

        dir_path, cache_path = self.build_bloom_prob_path(
            qa_value=qa_value,
            start_year=start_label,
            end_year=end_label,
            interval=interval,
            resolution=resolution,
            x_max=x_max,
            y_max=y_max,
        )
        os.makedirs(dir_path, exist_ok=True)

        source_mtime = max(
            os.path.getmtime(self.p_path),
            os.path.getmtime(self.e_path),
        )

        if (
            os.path.isfile(cache_path)
            and os.path.getmtime(cache_path) >= source_mtime
        ):
            cached, metadata = self._read_bloom_prob_cache(cache_path)

            return (
                cached,
                metadata["start_year"],
                metadata["end_year"],
                metadata["qa_value"],
            )

        fit = self._fit_bloom_kde(qa_value=qa_value, start_year=start_year, end_year=end_year)
        if fit is None:
            return None
        kde, start, end, qa_filtered_set = fit

        warnings.warn("Bloom probabilities need to be calculated. Depending on the lake size this may take a while.")

        xi = np.arange(-interval, x_max + interval, resolution)
        yi = np.arange(-interval, y_max + interval, resolution)

        # Evaluate the KDE in parallel by chunks of y-grid rows.
        n_workers = min(os.cpu_count() or 4, len(yi))
        row_chunks = np.array_split(yi, n_workers)

        with multiprocessing.Pool(
            processes=n_workers,
            initializer=_init_bloom_probability_worker,
            initargs=(kde,),
        ) as pool:
            density_parts = pool.map(
                _evaluate_bloom_probability_rows,
                [(rows, xi) for rows in row_chunks],
            )

        Zi = np.vstack(density_parts) * resolution ** 2

        cumsum = np.cumsum(np.cumsum(Zi, axis=0), axis=1)
        cumsum = np.pad(cumsum, ((1, 0), (1, 0)))

        n_cells = interval // resolution
        window_sum = (
            cumsum[n_cells:, n_cells:]
            - cumsum[:-n_cells, n_cells:]
            - cumsum[n_cells:, :-n_cells]
            + cumsum[:-n_cells, :-n_cells]
        )
        window_sum = np.clip(window_sum, 0, None)

        n_y, n_x = window_sum.shape
        x_low, y_low = np.meshgrid(xi[:n_x], yi[:n_y])

        probability_df = pd.DataFrame({
            "x_low": x_low.ravel(),
            "x_high": x_low.ravel() + interval,
            "y_low": y_low.ravel(),
            "y_high": y_low.ravel() + interval,
            "probability": window_sum.ravel(),
        })

        if self.save_format == "netcdf":
            self._write_bloom_prob_netcdf(cache_path, probability_df)
        else:
            probability_df.to_csv(cache_path, index=False)

        return probability_df, start, end, qa_filtered_set
    

    def write_output(self, file_path, bloom_prob_df):
        """Write cached bloom probability data to NetCDF."""
        columns = ['x_low', 'x_high','y_low','y_high','probability']

        with netCDF4.Dataset(file_path, "w") as ds:
            ds.createDimension("event", len(bloom_prob_df))

            for column in columns:
                values = bloom_prob_df[column].to_numpy()

                if column == "qa_column":
                    variable = ds.createVariable(column, "f4", ("event",))
                    variable[:] = values
                else:
                    variable = ds.createVariable(column, "f8", ("event",))
                    variable[:] = values

            ds.lake_id = str(self.lakeID)
            ds.version = str(self.version)
            ds.variable = str(self.variable)

    def read_metadata_from_path(self, file_path):
        """Extract bloom-probability metadata encoded in the cache filename."""
        filename = Path(file_path).name

        pattern = (
            r"^ID(?P<lake_id>[^_]+)_"
            r"(?P<start>[^_]+)_"
            r"(?P<end>[^_]+)_"
            r"(?P<qa>all|qa[0-9-]+)_"
            r"i(?P<interval>[^_]+)_"
            r"r(?P<resolution>[^_]+)_"
            r"x(?P<x_max>[^_]+)_"
            r"y(?P<y_max>[^.]+)\."
            r"(?P<extension>csv|nc)$"
        )

        match = re.match(pattern, filename)
        if match is None:
            raise ValueError(f"Invalid bloom probability cache filename: {filename}")

        values = match.groupdict()
        qa_label = values["qa"]

        return {
            "lake_id": values["lake_id"],
            "start_year": int(values["start"]),
            "end_year": int(values["end"]),
            "qa_value": (
                None
                if qa_label == "all"
                else {int(value) for value in qa_label.removeprefix("qa")}
            ),
            "interval": int(values["interval"]),
            "resolution": float(values["resolution"]),
            "x_max": float(values["x_max"]),
            "y_max": float(values["y_max"]),
        }

    def read_cached_metric(self, file_path):
        """Read cached bloom proabability data from NetCDF."""
        with netCDF4.Dataset(file_path, "r") as ds:
            req_cols = ['x_low', 'x_high','y_low','y_high','probability'] # used to maintain order
            required = set(req_cols)
            missing = required.difference(ds.variables)

            if missing:
                raise ValueError(
                    f"Bloom prob cache {file_path} is missing variables: {sorted(missing)}"
                )

            bloom_df = pd.DataFrame({req_var: np.asarray(ds.variables[req_var][:]) for req_var in required})
            return bloom_df[req_cols]

    def read_full_cache(self, file_path):
        """Read a bloom-probability cache and return its data and filename metadata."""
        if Path(file_path).suffix == ".nc":
            bloom_df = self._read_bloom_prob_nc(file_path)
        else:
            bloom_df = pd.read_csv(file_path)

        metadata = self._bloom_prob_metadata_from_path(file_path)
        return bloom_df, metadata