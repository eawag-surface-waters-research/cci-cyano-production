import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime
import os
import csv
import statistics
import warnings
from matplotlib.colors import ListedColormap, BoundaryNorm
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from csaps import csaps
from functions import unix_to_datetime, unix_to_datenum, datenum_to_datetime, remove_nan, define_year_range, plot_lake_outline, grab_metrics, grab_time_data, plot_map_data, set_labels, prep_dimark_data
import multiprocessing
from functools import partial
import shapely.ops as ops
from pyproj import CRS, Transformer
import geopandas
from shapely.prepared import prep
from shapely.geometry import Point
from numpy.lib.stride_tricks import sliding_window_view



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
    t_all = unix_to_datenum(time_raw)

    # Convert once, not per pixel
    time_dt = np.array(datenum_to_datetime(t_all))
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


class PhenologyVisualization:
        shapefile_path = None

        def __init__(self, extract_path, phenology_path):
                """Initialise a PhenologyVisualization instance for a single lake and variable.

                Parameters
                ----------
                extract_path : str
                    Path to the extract NetCDF file containing satellite observations.
                phenology_path : str
                    Path to the phenology parameter NetCDF file produced by the spline
                    fitting step. The lake ID, product version, and variable name are
                    derived automatically from the directory structure of this path.

                Raises
                ------
                Warning
                    If set_shapefile_path has not been called before instantiation.
                """
                if self.shapefile_path is None:
                        raise Warning("Please define your path to the lake CCI shapefile. This can be done at a class level using PhenologyVisualization.set_shapefile_path(your_path)")
                else:
                        self.geom = geopandas.read_file(self.shapefile_path)
                self.p_path = phenology_path
                self.e_path = extract_path
                self.info = (f"Version: {os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(phenology_path))))[1:]} \n",
                        f" Variable: {os.path.basename(os.path.dirname(os.path.dirname(phenology_path)))} \n" ,
                        f"Lake ID: {os.path.basename(phenology_path)[:-3]}")
                self.version = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(phenology_path))))[1:]
                self.variable = os.path.basename(os.path.dirname(phenology_path))
                self.methods = [method for method in dir(PhenologyVisualization) if callable(getattr(PhenologyVisualization, method)) and not method.startswith("__")]
                self.valid_coords = self.valid_index_pairs()
                self.out_folder =  os.path.dirname(os.path.dirname(os.path.dirname(self.p_path)))
                self.aggregation_df = None
                self.geom_shrunk = None
                self._extracted_globals = None
                self._pixel_cache = {}


        def shortname_to_name(self, short_name:str):
                """Return the full lake name for a given short_name identifier.

                Parameters
                ----------
                short_name : str
                    The short_name value to look up in the shapefile attribute table.

                Returns
                -------
                str
                    Formatted string containing the lake name.
                """
                df = self.geom
                return f"ID: {df.loc[df["short_name"]==short_name, "name"][1]}"

        def ID_to_name(self, id:int):
                """Return the lake name for a given numeric lake ID.

                Parameters
                ----------
                id : int
                    The numeric lake ID to look up in the shapefile attribute table.

                Returns
                -------
                str
                    The lake name string.
                """
                df = self.geom
                return list(df.loc[df["id"]==id, "name"])[0]

        def name_to_ID(self, name:str):
                """Return the lake ID for a given lake name.

                Parameters
                ----------
                name : str
                    The lake name to look up in the shapefile attribute table.

                Returns
                -------
                str
                    Formatted string containing the lake ID.
                """
                df = self.geom
                return f"ID: {list(df.loc[df["name"]==name, "id"])[0]}"

        def name_to_shortname(self, name:str):
                """Return the short_name for a given lake name.

                Parameters
                ----------
                name : str
                    The lake name to look up in the shapefile attribute table.

                Returns
                -------
                str
                    Formatted string containing the short_name.
                """
                df = self.geom
                return f"short_name: {list(df.loc[df["name"]==name, "short_name"])[0]}"

        def index_to_lat_lon(self, lat_index, lon_index):
                """Return the geographic coordinates for a grid index pair.

                Parameters
                ----------
                lat_index : int
                    Row index in the extract grid.
                lon_index : int
                    Column index in the extract grid.

                Returns
                -------
                str
                    Formatted string with the latitude and longitude values.
                """
                with netCDF4.Dataset(self.e_path) as nc:
                        lats = nc.variables["lat"][:]
                        lons = nc.variables["lon"][:]
                        lat = lats[lat_index]
                        lon = lons[lon_index]
                return f"Lat, Lon: {lat}, {lon}"
        
        
        @classmethod
        def set_shapefile_path(cls, path: str):
                """Set the CCI lake shapefile path used by all instances.

                Must be called once before any PhenologyVisualization instance is created.

                Parameters
                ----------
                path : str
                    Absolute or relative path to the CCI lake shapefile (.shp).
                """
                cls.shapefile_path = path



        def valid_index_pairs(self):
                """Return grid index pairs that have more than one valid QA-passing observation.

                Reads the extract NetCDF and identifies pixels where the observation count
                (non-fill, QA==0) exceeds one — the minimum required for spline fitting.

                Returns
                -------
                list of tuple of int
                    List of (row, col) index pairs with sufficient valid observations.
                """
                with netCDF4.Dataset(self.e_path) as nc:
                        variable = getattr(nc, "variable")
                        qa_variable = getattr(nc, "qa")

                        values = np.asarray(nc.variables[variable][:])
                        qa = np.asarray(nc.variables[qa_variable][:])

                        valid_mask = (values != -9999) & (qa == 0)
                        valid_counts = np.sum(valid_mask, axis=0)

                        return [tuple(int(x) for x in idx) for idx in np.argwhere(valid_counts > 1)]
                

        def create_DataFrame(self, latitude, longitude):
                """Build a combined long-format DataFrame of all phenology variables for a single pixel.

                Reads all _x (time) and _y (value) variable pairs from the phenology NetCDF
                and concatenates them into a single DataFrame.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel in the grid.
                longitude : int
                    Column (lon) index of the pixel in the grid.

                Returns
                -------
                pandas.DataFrame
                    Long-format DataFrame with columns Value, Variable, latitude, longitude
                    and a datetime index named Time.
                """
                p = netCDF4.Dataset(self.p_path)
                exclude = ['lat','lon','smoothing_parameter','trgs_qa','data_gap_start','data_gap_end']
                l = list(set(list(p.variables))-set(exclude))
                variables_x = sorted([i for i in l if i[-1]== "x"])
                variables_y = sorted([i for i in l if i[-1]== "y"])
                lat = np.array(p.variables["lat"])
                lon = np.array(p.variables["lon"])

                result = {}

                for x,y in zip(variables_x, variables_y):
                        var_x = unix_to_datetime(remove_nan(p[x][latitude,longitude,:]))
                        var_y = remove_nan(p[y][latitude,longitude,:])
                        var_label = [x[:-2]]*len(var_y)
                        df = pd.DataFrame({"Value":var_y,
                                        "Variable": var_label,
                                        "latitude": lat[latitude],
                                        "longitude": lon[longitude]},
                                        index = var_x)
                        df.index.names = ["Time"]
                        result[y[:-2]] = df
                combined_df = pd.concat(result.values())

                return combined_df


        def shrink_geometry_by_1km(self, geom):
                """Return a lake geometry eroded inward by 1 km.

                Projects the geometry to a local Azimuthal Equidistant CRS centred on the
                lake centroid, applies a -1000 m buffer, then reprojects back to WGS84.
                The result is also stored in self.geom_shrunk.

                Parameters
                ----------
                geom : shapely.geometry.base.BaseGeometry
                    WGS84 lake geometry (Polygon or MultiPolygon) to shrink.

                Returns
                -------
                shapely.geometry.base.BaseGeometry
                    Inward-buffered WGS84 geometry.
                """
                # centroid in lon/lat
                lon0, lat0 = geom.centroid.x, geom.centroid.y

                # local Azimuthal Equidistant projection centered on the lake
                local_crs = CRS.from_proj4(
                        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m +no_defs"
                )

                wgs84 = CRS.from_epsg(4326)

                to_local = Transformer.from_crs(wgs84, local_crs, always_xy=True).transform
                to_wgs84 = Transformer.from_crs(local_crs, wgs84, always_xy=True).transform

                geom_local = ops.transform(to_local, geom)
                geom_shrunk_local = geom_local.buffer(-1000)   # minus 1000 m = inward 1 km
                geom_shrunk = ops.transform(to_wgs84, geom_shrunk_local)

                self.geom_shrunk = geom_shrunk
                return geom_shrunk
        

       
        @staticmethod
        def compute_metric_score(coord, start, end, metric_to_compute= ["values_per_pixel", "r2", "MAD", "RMSE", "correlation"]):
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
                metric_to_compute : list of str
                    Single-element list naming the metric. One of:
                    ['values_per_pixel'], ['r2'], ['MAD'], ['RMSE'], ['correlation'].

                Returns
                -------
                tuple
                    ((i, j), metric_value) where metric_value is an int for
                    'values_per_pixel' and a float (or np.nan) for all others.
                """
                if metric_to_compute==["values_per_pixel"]:
                        i,j = coord
                        values_all = _GLOBALS["values_all"]
                        qa_all = _GLOBALS["qa_all"]
                        years_all = _GLOBALS["years_all"]

                        values = values_all[:, i, j]
                        qa_values = qa_all[:, i, j]

                        mask = (values != -9999) & (qa_values == 0)
                        values_m = values[mask]

                        years_m = years_all[mask]

                        function_start, function_end = define_year_range(start= start, end= end, years= years_m)

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

                                function_start, function_end = define_year_range(start= start, end = end, years = years_m)

                                y_pred =csaps(time_m, values_m, time_m, smooth=smoothing)
                                y_true = values_m

                                mask_sub = (years_m>= function_start) & (years_m <=function_end)
                                valid = np.isfinite(y_true) & np.isfinite(y_pred)

                                if valid.sum()<1:
                                        warnings.warn(f"Check data for lat, lon indices:{(i, j)}, perhaps smoothing parameter is nan or duplicates in time axis.")

                                combined_mask = valid & mask_sub

                                if combined_mask.sum() > 1:
                                        if metric_to_compute == ["r2"]:
                                                metric = r2_score(y_true[combined_mask], y_pred[combined_mask])
                                        elif metric_to_compute == ["MAD"]:
                                                metric = np.median(np.abs(y_true[combined_mask]-y_pred[combined_mask]))
                                        elif metric_to_compute == ["RMSE"]:
                                                metric = np.sqrt(mean_squared_error(y_true[combined_mask], y_pred[combined_mask]))
                                        elif metric_to_compute == ["correlation"]:
                                                metric, _ = pearsonr(y_true[combined_mask], y_pred[combined_mask])
                                        else:
                                                raise ValueError("please enter a valid metric")

                                else:
                                        warnings.warn(f"Not enough valid data in selected date range for indices {(i,j)}")
                                        metric = np.nan
                        else:
                                metric = np.nan
                        return (i,j), metric




        def build_metric_path(self, metric_name, start, end):
                """Return the output directory and file path for a cached metric CSV.

                The filename encodes the year range: full_ts.csv for the complete series,
                ts_{start}_to_{end}.csv otherwise. Sentinel values 0 and 9999 are replaced
                by 2002 and 2024 respectively in the filename.

                Parameters
                ----------
                metric_name : str
                    Name of the metric (e.g. 'r2', 'MAD', 'RMSE', 'correlation',
                    'values_per_pixel').
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                base : str
                    Directory path where the CSV will be written.
                file_path : str
                    Full path to the CSV file.
                """
                base = os.path.join(self.out_folder, "calculated_values", "metrics", metric_name, f"v{self.version}", self.variable)
                if start == 0 and end == 9999:
                        fname = "full_ts.csv"
                elif start == 0:
                        fname = f"ts_{2002}_to_{end}.csv"
                elif end == 9999:
                        fname = f"ts_{start}_to_{2024}.csv"
                else:
                        fname = f"ts_{start}_to_{end}.csv"
                
                return base, os.path.join(base, fname)
                
        

        def compute_and_cache_metric(self, metric_name, col_name, compute_fn, start, end):
                """Compute a metric for all valid pixels in parallel and cache to CSV.

                On the first call the metric is computed using a multiprocessing pool
                (3 workers) and written to CSV. Subsequent calls with the same metric
                and year range load directly from the cached file.

                Parameters
                ----------
                metric_name : str
                    Name of the metric, passed to build_metric_path and compute_fn.
                col_name : str
                    Column name used when reading or writing the CSV cache.
                compute_fn : callable
                    Worker function (typically compute_metric_score) executed by each
                    pool worker.
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                dict
                    Mapping of (i, j) tuples to metric values for all valid pixels.
                """
                dir_path, file_path = self.build_metric_path(metric_name, start, end)
                os.makedirs(dir_path, exist_ok=True)
                if os.path.isfile(file_path):
                        df = pd.read_csv(file_path)
                        return dict(zip(zip(df["i"], df["j"]), df[col_name]))
                warnings.warn(f"{metric_name} need to be calculated. Depending on the lake size this may take a while.")
                workers = partial(compute_fn, start=start, end=end, metric_to_compute = [metric_name])
                with multiprocessing.Pool(initializer=_init_worker, initargs=(self.p_path, self.e_path), processes=3) as pool:
                        result = pool.map(workers, self.valid_coords)
                data = dict(result)
                t_df = pd.DataFrame([(i, j, v) for (i, j), v in data.items()], columns=["i", "j", col_name])
                t_df.to_csv(file_path, index=False)
                return data


        def r2_scores(self, time_split):
                """Return R² scores for all valid pixels over a single time window.

                Parameters
                ----------
                time_split : list of tuple of int
                    Single-element list containing one (start, end) year tuple.
                    Use [(0, 9999)] for the full time series.

                Returns
                -------
                dict
                    Mapping of (i, j) grid index tuples to R² float values.
                """
                for start, end in time_split:
                        return self.compute_and_cache_metric(metric_name="r2", col_name="r2_scores", compute_fn=PhenologyVisualization.compute_metric_score, start=start, end=end)

        def MAD_scores(self, time_split):
                """Return Median Absolute Deviation scores for all valid pixels over a single time window.

                Parameters
                ----------
                time_split : list of tuple of int
                    Single-element list containing one (start, end) year tuple.
                    Use [(0, 9999)] for the full time series.

                Returns
                -------
                dict
                    Mapping of (i, j) grid index tuples to MAD float values.
                """
                for start, end in time_split:
                        return self.compute_and_cache_metric(metric_name="MAD", col_name="mad_scores",compute_fn= PhenologyVisualization.compute_metric_score,start= start,end= end)

        def RMSE_scores(self, time_split):
                """Return Root Mean Squared Error scores for all valid pixels over a single time window.

                Parameters
                ----------
                time_split : list of tuple of int
                    Single-element list containing one (start, end) year tuple.
                    Use [(0, 9999)] for the full time series.

                Returns
                -------
                dict
                    Mapping of (i, j) grid index tuples to RMSE float values.
                """
                for start, end in time_split:
                        return self.compute_and_cache_metric(metric_name="RMSE", col_name="rmse_scores", compute_fn=PhenologyVisualization.compute_metric_score, start=start, end=end)

        def correlation_scores(self, time_split):
                """Return Pearson correlation scores for all valid pixels over a single time window.

                Parameters
                ----------
                time_split : list of tuple of int
                    Single-element list containing one (start, end) year tuple.
                    Use [(0, 9999)] for the full time series.

                Returns
                -------
                dict
                    Mapping of (i, j) grid index tuples to Pearson r float values.
                """
                for start, end in time_split:
                        return self.compute_and_cache_metric(metric_name="correlation", col_name="correlation_scores", compute_fn=PhenologyVisualization.compute_metric_score, start=start, end=end)

        def values_per_pixel(self, time_split):
                """Return valid observation counts for all pixels over a single time window.

                Parameters
                ----------
                time_split : list of tuple of int
                    Single-element list containing one (start, end) year tuple.
                    Use [(0, 9999)] for the full time series.

                Returns
                -------
                dict
                    Mapping of (i, j) grid index tuples to integer observation counts.
                """
                for start, end in time_split:
                        return self.compute_and_cache_metric(metric_name="values_per_pixel", col_name="number_of_values",compute_fn= PhenologyVisualization.compute_metric_score, start=start,end= end)



        def spatial_aggregation(self):
                """Compute or load per-pixel 3×3 neighbourhood median values for all timesteps.

                For each valid interior pixel and each timestep, computes the median of all
                non-fill, QA==0 observations within the 3×3 pixel neighbourhood using
                stride-trick windowing. Border pixels are excluded as they cannot form a
                complete 3×3 window.

                The result is stored in self.aggregation_df and written to a CSV for reuse.
                If the CSV already exists, it is loaded directly without recomputation.

                Returns
                -------
                None
                    Result is stored in self.aggregation_df as a pandas.DataFrame with
                    columns: time, i, j, lat, lon, MA_value.
                """

                out_dir = os.path.join(
                        self.out_folder,
                        "calculated_values", "spatial_aggregation_values",
                        f"v{self.version}", self.variable,
                )
                os.makedirs(out_dir, exist_ok=True)

                file_path = os.path.join(out_dir, "aggregation_background_values.csv")

                if os.path.isfile(file_path):
                        self.aggregation_df = pd.read_csv(file_path)
                        return

                warnings.warn("spatial aggregation needs to be calculated. Depending on lake size this could take a while.")

                with netCDF4.Dataset(self.e_path) as nc:
                        lat = np.asarray(nc.variables["lat"][:])
                        lon = np.asarray(nc.variables["lon"][:])
                        t_all = unix_to_datenum(nc.variables["time"][:])

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
                                aggregation_df = pd.DataFrame(
                                        columns=["time", "i", "j", "lat", "lon", "MA_value"]
                                )
                                aggregation_df.to_csv(file_path, index=False)
                                self.aggregation_df = aggregation_df
                                return

                        i_idx = coords[:, 0]
                        j_idx = coords[:, 1]

                        # indices for median_grid, which is smaller by 1 border cell each side
                        ii = i_idx - 1
                        jj = j_idx - 1

                        lat_vals = lat[i_idx]
                        lon_vals = lon[j_idx]

                        frames = []   # <- this must exist before the loop

                        for n in range(ntime):
                                data_n = np.asarray(data_var[n], dtype=np.float32)
                                qa_n = np.asarray(qa_var[n])

                                # shape: (nlat-2, nlon-2, 3, 3)
                                data_windows = sliding_window_view(data_n, (3, 3))
                                qa_windows = sliding_window_view(qa_n, (3, 3))

                                valid_mask = (data_windows != -9999) & (qa_windows == 0)

                                masked = data_windows.astype(np.float32, copy=True)
                                masked[~valid_mask] = np.nan

                                # shape: (nlat-2, nlon-2)
                                median_grid = np.nanmedian(masked, axis=(-2, -1))

                                ma_values = median_grid[ii, jj]

                                frames.append(
                                        pd.DataFrame({
                                        "time": np.full(len(coords), t_all[n]),
                                        "i": i_idx,
                                        "j": j_idx,
                                        "lat": lat_vals,
                                        "lon": lon_vals,
                                        "MA_value": ma_values,
                                        })
                                )

                aggregation_df = pd.concat(frames, ignore_index=True)
                aggregation_df.to_csv(file_path, index=False)
                self.aggregation_df = aggregation_df



        def pixel_map(self, latitude, longitude, ax):
                """Plot a grayscale coverage map with the selected pixel marked.

                Reads the summary grid from the extract NetCDF and displays it in
                grayscale, masking cells with fewer than 2 valid observations. The
                selected pixel is overlaid as a red star.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel to highlight.
                longitude : int
                    Column (lon) index of the pixel to highlight.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the map.

                Returns
                -------
                matplotlib.image.AxesImage
                    The imshow image object.
                """
                with netCDF4.Dataset(self.e_path) as nc:
                        summary = np.array(nc.variables["summary"][:, :])


                # mask invalid cells
                masked_summary = np.ma.masked_where(summary <= 2, summary)

                cmap = ListedColormap("gray")

                # plot as geographic map
                im = ax.imshow(masked_summary, cmap=cmap, aspect="auto", origin="lower")

                # plot selected pixel
                ax.plot(longitude, latitude, "r*", markersize=14, zorder=5, label="Pixel")

                ax.set_xlabel("Lon index")
                ax.set_ylabel("Lat index")
                text_str = f"Pixel Location\n Lake ID:{os.path.basename(self.p_path)[:-3]}"
                ax.set_title(text_str)
                ax.legend()

                return im

        def interactive_pixel_map(self, ax):
                """Plot a clickable coverage map for pixel selection.

                Displays the summary coverage grid in grayscale. Clicking a valid cell
                prints its (lat_idx, lon_idx) and marks it with a dot and text label.
                Requires an interactive matplotlib backend (e.g. %matplotlib widget).

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Axes on which to draw the map.

                Returns
                -------
                int
                    Matplotlib canvas connection ID, which can be passed to
                    ax.figure.canvas.mpl_disconnect() to remove the click handler.
                """
                with netCDF4.Dataset(self.e_path) as nc:
                        summary = np.array(nc.variables["summary"][:, :])

                masked_summary = np.ma.masked_where(summary <= 2, summary)

                cmap = ListedColormap("gray")

                ax.imshow(masked_summary, cmap=cmap, aspect="auto", origin="lower")

                ax.set_xlabel("Lon index")
                ax.set_ylabel("Lat index")
                ax.set_title(f"Pixel Location\n Lake ID: {os.path.basename(self.p_path)[:-3]}")

                def on_click(event):
                        """Label the clicked valid pixel with its lat/lon indices."""
                        print("click detected")

                        if event.inaxes is not ax or event.xdata is None or event.ydata is None:
                                return

                        lon_idx = int(round(event.xdata))
                        lat_idx = int(round(event.ydata))

                        lat_idx = max(0, min(lat_idx, summary.shape[0] - 1))
                        lon_idx = max(0, min(lon_idx, summary.shape[1] - 1))

                        if summary[lat_idx, lon_idx] <= 2:
                                print("invalid cell")
                                return

        
                        ax.text(
                        lon_idx, lat_idx,
                        f"{lat_idx},{lon_idx}",
                        color="black",
                        fontsize=10,
                        ha="center", va="center"
                )

                        ax.plot(lon_idx, lat_idx, "ro", markersize=6)
                        ax.figure.canvas.draw_idle()

                cid = ax.figure.canvas.mpl_connect("button_press_event", on_click)

                return cid

        def metric_map(self, metric_scores:dict, metric_str:str, fig, ax, colormap = None, colorbar_extent= [0,1]):
                """Plot a spatial heatmap of a precomputed per-pixel metric.

                Pixels outside the 1 km-inset lake boundary are masked. The lake outline
                is overlaid from the shapefile geometry.

                Parameters
                ----------
                metric_scores : dict
                    Mapping of (i, j) grid index tuples to metric float values, as
                    returned by r2_scores, MAD_scores, RMSE_scores, etc.
                metric_str : str
                    Label used for the colorbar and plot title (e.g. 'R$^2$', 'RMSE').
                fig : matplotlib.figure.Figure
                    Figure used to attach the colorbar.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the heatmap.
                colormap : str or matplotlib.colors.Colormap, optional
                    Colormap passed to imshow. Defaults to 'RdYlBu'.
                colorbar_extent : list of float, optional
                    [vmin, vmax] for the colorbar. Defaults to [0, 1].

                Returns
                -------
                matplotlib.image.AxesImage
                    The imshow image object.

                Raises
                ------
                ValueError
                    If the lake ID derived from p_path is not found in the shapefile.
                """
                lake_id = int(os.path.basename(self.p_path)[:-3])
                lake_row = self.geom[self.geom["id"] == lake_id]
                if lake_row.empty:
                        raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
                
                geom = lake_row.geometry.iloc[0]
                buffered_geom = self.shrink_geometry_by_1km(geom)
                buffered_geom_prepared = prep(buffered_geom)

                map_data, extent = grab_metrics(self.e_path, metric_scores, buffered_geom_prepared)

                im = plot_map_data(colormap, map_data, extent, ax, cmap_extent=colorbar_extent)
                plot_lake_outline(geometry=geom, ax=ax)
                set_labels(ax, fig, im,
                           title=f"{metric_str}-Scores for Lake: ID {lake_id}",
                           colorbar_label=metric_str)

                return im






        def interactive_metric_map(self, metric_scores, metric_str:str, fig, ax):
                """Plot a clickable spatial heatmap of a precomputed per-pixel metric.

                Like metric_map but with click interaction: clicking a pixel prints and
                labels its (lat_idx, lon_idx). Colorbar range is fixed to [0, 1].
                Requires an interactive matplotlib backend (e.g. %matplotlib widget).

                Parameters
                ----------
                metric_scores : dict
                    Mapping of (i, j) grid index tuples to metric float values.
                metric_str : str
                    Label used for the colorbar and plot title.
                fig : matplotlib.figure.Figure
                    Figure used to attach the colorbar.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the heatmap.

                Returns
                -------
                int
                    Matplotlib canvas connection ID for removing the click handler.

                Raises
                ------
                ValueError
                    If the lake ID derived from p_path is not found in the shapefile.
                """
                lake_id = int(os.path.basename(self.p_path)[:-3])

                lake_row = self.geom[self.geom["id"] == lake_id]
                if lake_row.empty:
                        raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
                geom = lake_row.geometry.iloc[0]
                buffered_geom = self.shrink_geometry_by_1km(geom)
                buffered_geom_prepared = prep(buffered_geom)

                map_data, extent = grab_metrics(self.e_path, metric_scores, buffered_geom_prepared)
                g = self._load_extracted_globals()
                lats = g["lat"]
                lons = g["lon"]

                im = plot_map_data(None, map_data, extent, ax, cmap_extent=[0, 1])
                plot_lake_outline(geometry=geom, ax=ax)
                set_labels(ax, fig, im,
                           title=f"{metric_str}-Scores for Lake: ID {lake_id}",
                           colorbar_label=metric_str)

                def on_click(event):
                        """Label the clicked metric-map pixel with its lat/lon indices."""
                        if event.inaxes is not ax or event.xdata is None or event.ydata is None:
                                return

                        clicked_lon = event.xdata
                        clicked_lat = event.ydata

                        lon_idx = int(np.abs(lons - clicked_lon).argmin())
                        lat_idx = int(np.abs(lats - clicked_lat).argmin())

                        ax.text(
                        lons[lon_idx], lats[lat_idx],
                        f"{lat_idx},{lon_idx}",
                        color="black",
                        fontsize=10,
                        ha="center", va="center"
                        )
                        ax.plot(lons[lon_idx], lats[lat_idx], "ro", markersize=6)
                        fig.canvas.draw_idle()

                cid = fig.canvas.mpl_connect("button_press_event", on_click)
                return cid
        

        def time_map(self, fig, ax, year, peaks=True, max = True, colorbar=True):
                """Map the day-of-year of a phenological event across all pixels for one year.

                For each valid pixel within the 1 km-inset lake boundary, extracts the
                summer peak or green-up midpoint DOY (restricted to DOY 160–250) for the
                given year and displays it as a spatial heatmap with a fixed colorbar.

                Parameters
                ----------
                fig : matplotlib.figure.Figure
                    Figure used to attach the colorbar.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the heatmap.
                year : int
                    Calendar year for which to map the phenological event.
                peaks : bool, optional
                    If True (default), map summer peaks (pks_x / pks_y).
                    If False, map green-up midpoints (green_up_mid_x / green_up_mid_y).
                max : bool, optional
                    If True (default), use the highest-amplitude peak in the DOY window.
                    If False, use the first peak in the window.

                Returns
                -------
                matplotlib.image.AxesImage
                    The imshow image object.
                """
                lake_id = int(os.path.basename(self.p_path)[:-3])

                lake_row = self.geom[self.geom["id"] == lake_id]
                geom = lake_row.geometry.iloc[0]
                buffered_geom = self.shrink_geometry_by_1km(geom)
                buffered_geom_prepared = prep(buffered_geom)

                var_x = "pks_x" if peaks else "green_up_mid_x"
                var_y = "pks_y" if peaks else "green_up_mid_y"
                extrema_label = "Peak" if peaks else "Green Mid Up"

                map_data, extent = grab_time_data(self.e_path, self.p_path, self.valid_coords,
                                                  buffered_geom_prepared, var_x, var_y, year, max)
                im = plot_map_data("rainbow", map_data, extent, ax, cmap_extent=[160, 250])
                plot_lake_outline(geometry=geom, ax=ax)
                if colorbar:
                        set_labels(ax, fig, im,
                                   title=f"{extrema_label} Day of Year\n Lake ID: {lake_id}\n Year: {year}",
                                   colorbar_label="Day of Year",
                                   colorbar_ticks=[160, 185, 215, 250])
                else:
                        ax.set_title(f"{extrema_label} Day of Year\n Lake ID: {lake_id}\n Year: {year}")
                        ax.set_xlabel("Lon index")
                        ax.set_ylabel("Lat index")
                        ax.legend()
                return im
        

        def time_map_panel(self, years, nrow, ncol, peaks = True, max = True):
                extrema_label = "Peak" if peaks else "Green Mid Up"
                fig, axs = plt.subplots(nrow, ncol, constrained_layout=True, squeeze=False, figsize=(ncol * 5, nrow * 4))
                im = None
                for year, ax in zip(years, axs.flatten()):
                        im = self.time_map(fig=fig, ax=ax, year=year, peaks=peaks, max=max, colorbar=False)
                        ax.set_title(str(year), fontsize=20)
                        ax.set_ylabel("Lat index", fontsize=15)
                        ax.set_xlabel("Lon index", fontsize=15)
                        ax.tick_params(labelsize=15)

                for ax in axs.flatten()[len(years):]:
                        ax.set_visible(False)

                if im is not None:
                        cbar = fig.colorbar(im, ax=axs.ravel().tolist(), location="right", shrink=0.8)
                        cbar.set_label("Day of Year", fontsize=20)
                        cbar.set_ticks([160, 185, 215, 250])
                        cbar.ax.tick_params(labelsize=15)

                fig.suptitle(f"{extrema_label} Day of Year", fontsize=25)
                
                plt.show()



        def single_day_map(self, date):
                """Plot a spatial map of chlorophyll values for a single observation date.

                Creates a new figure internally. Values failing the fill-value or QA filter
                are masked. The lake outline is overlaid from the shapefile geometry.

                Parameters
                ----------
                date : datetime.datetime
                    A UTC-aware timestamp matching one entry in the extract time axis.
                    Typically obtained from a pixel time series, e.g.::

                        series = vis.load_pixel_data(i, j)
                        vis.single_day_map(series.idxmax())

                Returns
                -------
                matplotlib.image.AxesImage
                    The imshow image object.

                Raises
                ------
                ValueError
                    If date is not found in the extract time axis.
                """
                date = pd.Timestamp(date)
                if date.tzinfo is None:
                        date = date.tz_localize("UTC")
                else:
                        date = date.tz_convert("UTC")
                g = self._load_extracted_globals()
                t_all = g["t_all"]
                idx = np.argwhere(datenum_to_datetime(t_all) == date)
                if len(idx) == 0:
                        raise ValueError(f"Date {date} not found in the extract time axis.")
                time_index = idx[0, 0]

                with netCDF4.Dataset(self.e_path) as nc:
                        values = np.array(nc.variables[g["variable"]][time_index, :, :])
                        qa     = np.array(nc.variables[g["qa"]][time_index, :, :])
                        lats   = nc.variables["lat"][:]
                        lons   = nc.variables["lon"][:]

                mask = (values != -9999) & (qa == 0)
                lake_id = int(os.path.basename(self.p_path)[:-3])
                lake_row = self.geom[self.geom["id"] == lake_id]
                geom = lake_row.geometry.iloc[0]
                buffered_geom = self.shrink_geometry_by_1km(geom)
                buffered_geom_prepared = prep(buffered_geom)

                map_data = np.full(values.shape, np.nan)
                for i in range(values.shape[0]):
                        for j in range(values.shape[1]):
                                if mask[i, j] and buffered_geom_prepared.contains(Point(lons[j], lats[i])):
                                        map_data[i, j] = values[i, j]
                fig, ax = plt.subplots(1, 1, figsize=(10, 5))
                im = ax.imshow(map_data, cmap="winter", aspect="auto", origin="lower",
                               extent=[lons.min(), lons.max(), lats.min(), lats.max()])

                label = False
                if geom.geom_type == "Polygon":
                        x, y = geom.exterior.xy
                        ax.plot(x, y, color="black", linewidth=1, label="Lake Outline")
                elif geom.geom_type == "MultiPolygon":
                        for poly in geom.geoms:
                                x, y = poly.exterior.xy
                                ax.plot(x, y, color="black", linewidth=1,
                                        label="Lake Outline" if not label else None)
                                label = True

                fig.colorbar(im, orientation='vertical', label=f'{self.variable} (ug/L)')
                plt.xticks([])
                plt.yticks([])
                plt.title(date.strftime("%Y-%m-%d"))
                return im




        def _load_extracted_globals(self):
                """Lazily load and cache shared arrays from the extract dataset.

                On the first call, opens the extract NetCDF and stores lat, lon, the full
                time array as datenums, and the variable and QA attribute names. Subsequent
                calls return the cached dict without reopening the file.

                Returns
                -------
                dict
                    Dictionary with keys: 'lat', 'lon', 't_all', 'variable', 'qa'.
                """
                if self._extracted_globals is None:
                        with netCDF4.Dataset(self.e_path) as nc:
                                self._extracted_globals= {
                                        "lat": np.asarray(nc.variables["lat"]),
                                        "lon": np.asarray(nc.variables["lon"]),
                                        "t_all":    unix_to_datenum(nc.variables["time"]),
                                        "variable": getattr(nc, "variable"),
                                        "qa":       getattr(nc, "qa"),
                                }
                return self._extracted_globals
        
        def _load_pixel_data(self, i,j):
                """Lazily load and cache all phenology arrays for a single pixel.

                On the first call for (i, j), reads values, QA, smoothing parameter, peaks,
                troughs, and green-up/green-down midpoints from both NetCDF files. Results
                are cached in self._pixel_cache for reuse across subsequent calls.

                Parameters
                ----------
                i : int
                    Row (lat) index of the pixel.
                j : int
                    Column (lon) index of the pixel.

                Returns
                -------
                dict
                    Dictionary with keys: 'values', 'qa', 'smoothing', 'pks_x', 'pks_y',
                    'trgs_x', 'trgs_y', 'midUP_x', 'midUP_y', 'midDOWN_x', 'midDOWN_y'.
                """
                if (i,j) not in self._pixel_cache:
                        g = self._load_extracted_globals()
                        with netCDF4.Dataset(self.e_path) as nc:
                                values = np.array(nc.variables[g["variable"]][:, i, j])
                                qa     = np.array(nc.variables[g["qa"]][:, i, j])
                        with netCDF4.Dataset(self.p_path) as nc:
                                smoothing = float(nc.variables["smoothing_parameter"][i, j])
                                pks_x = unix_to_datetime(remove_nan(nc.variables["pks_x"][i, j, :]))
                                pks_y = remove_nan(nc.variables["pks_y"][i, j, :])
                                trgs_x = unix_to_datetime(remove_nan(nc.variables["trgs_x"][i, j, :]))
                                trgs_y = remove_nan(nc.variables["trgs_y"][i, j, :])
                                midUP_x    = unix_to_datetime(remove_nan(nc.variables["green_up_mid_x"][i, j, :]))
                                midUP_y    = remove_nan(nc.variables["green_up_mid_y"][i, j, :])
                                midDOWN_x  = unix_to_datetime(remove_nan(nc.variables["green_down_mid_x"][i, j, :]))
                                midDOWN_y  = remove_nan(nc.variables["green_down_mid_y"][i, j, :])
                        self._pixel_cache[(i,j)] = {
                                "values": values, "qa": qa, "smoothing": smoothing, "pks_x":pks_x, "pks_y": pks_y,
                                "trgs_x": trgs_x, "trgs_y": trgs_y,  "midUP_x": midUP_x, "midUP_y": midUP_y,
                                 "midDOWN_x": midDOWN_x, "midDOWN_y": midDOWN_y,
                        }
                return self._pixel_cache[(i,j)]

        def load_pixel_data(self, i, j):
                """Return a datetime-indexed Series of valid observations for pixel (i, j).

                Filters the raw extract time series to keep only observations that pass
                both the fill-value check (value != -9999) and QA flag == 0.

                Parameters
                ----------
                i : int
                    Row (lat) index of the pixel.
                j : int
                    Column (lon) index of the pixel.

                Returns
                -------
                pandas.Series
                    Float values indexed by datetime, containing only valid (QA==0) observations.
                """
                g = self._load_extracted_globals()
                px = self._load_pixel_data(i, j)
                t_all = g["t_all"]
                mask     = (px["values"] != -9999) & (px["qa"] == 0)
                values_m = px["values"][mask]
                time_dt   = datenum_to_datetime(t_all[mask])
                return pd.Series(index=time_dt,data=values_m)


        def extrema_plot(self, latitude, longitude, ax,  peak = True, aggregation= False,  start = 0, end = 9999, background_pts = True, purple_chla21= False, show_legend = True):
                """Plot detected peaks or troughs as a stem plot with optional background scatter.

                Displays summer peaks or winter troughs for the pixel at (latitude, longitude)
                as vertical stems. Background observations may be shown as a raw scatter or
                3×3 spatial median (aggregation=True). Negative values are flagged with red
                crosses and trigger a warning.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the plot.
                peak : bool, optional
                    If True (default), plot summer peaks. If False, plot troughs.
                aggregation : bool, optional
                    If True, show the 3×3 neighbourhood median instead of raw scatter.
                    Requires spatial_aggregation() to have been called or available cache.
                start : int, optional
                    First year to display (inclusive). 0 = earliest in the series.
                end : int, optional
                    Last year to display (inclusive). 9999 = latest in the series.
                background_pts : bool, optional
                    If True (default), show scatter background observations.
                    Cannot be False when aggregation is True.
                purple_chla21 : bool, optional
                    Unused colour-override flag kept for API compatibility.

                Returns
                -------
                float or None
                    The upper y-axis limit set for the plot, or None if there is no data.

                Raises
                ------
                ValueError
                    If background_pts is False and aggregation is True simultaneously.
                """


                g = self._load_extracted_globals()
                px = self._load_pixel_data(latitude, longitude)
                lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
                smoothing = px["smoothing"]
                if peak:
                        pks_x, pks_y = px["pks_x"], px["pks_y"]
                        mask_pks = np.array([(d.year <= end) & (d.year >= start) for d in pks_x])
                        x_sub = pks_x[mask_pks]
                        y_sub = pks_y[mask_pks]

                else:
                        trgs_x, trgs_y = px["trgs_x"], px["trgs_y"]
                        mask_pks = np.array([(d.year <= end) & (d.year >= start) for d in trgs_x])
                        x_sub = trgs_x[mask_pks]
                        y_sub = trgs_y[mask_pks]

                mask     = (px["values"] != -9999) & (px["qa"] == 0)
                values_m = px["values"][mask]
                time_m   = t_all[mask]



                if len(values_m) > 1:

                        limits = sorted(datenum_to_datetime(time_m))
                        if start == 0:
                                function_start = min(limits).year
                        else:
                                function_start =start
                        if end == 9999:
                                function_end= max(limits).year
                        else:
                                function_end = end
                        neg_values_sub =[]
                        neg_label_before = False
                        phenology_name = os.path.basename(os.path.dirname(self.p_path))
                        if purple_chla21:
                                if phenology_name == "phycocyanin":
                                        label = "phyco"
                                        color = "blue"
                                elif phenology_name == "chla_mean":
                                        label = "chla v2.1"
                                        color = "purple"
                                else:
                                        label = "chla v3.1"
                                        color = "green"

                        else:
                                if phenology_name == "phycocyanin":
                                        label = "phyco"
                                        color = "blue"
                                elif phenology_name == "chla_mean":
                                        label = "chla v2.1"
                                        color = "lightgreen"
                                else:
                                        label = "chla v3.1"
                                        color = "green"

                        if not background_pts and aggregation:
                                raise ValueError("Either aggreagte background points or not plot them at all.")

                        if aggregation:
                                if self.aggregation_df is None:
                                        self.spatial_aggregation()

                                background_sub = self.aggregation_df[(self.aggregation_df["i"]==latitude) & (self.aggregation_df["j"]==longitude)]
                                background_time = background_sub["time"]
                                background_time = background_time.to_numpy()

                                background_values = background_sub["MA_value"]

                                ax.scatter(datenum_to_datetime(background_time), background_values, color=color, alpha=0.3, s=10, label=f"{label} Data")
                        elif background_pts and not aggregation:
                                ax.scatter(datenum_to_datetime(time_m), values_m, color=color, alpha=0.3, s=10, label=f"{label} Data")
                        else:
                                pass

                        extrema_label = "Peaks" if peak else "Troughs"
                        ax.stem(x_sub, y_sub, linefmt=color, label=f"{label} {extrema_label}", basefmt = " ")
                        if (y_sub < 0).any():
                                mask =  y_sub<0
                                pks_x_neg_before = x_sub[mask]
                                pks_y_neg_before = y_sub[mask]
                                label = "Negative Value" if not neg_label_before else None
                                ax.scatter(pks_x_neg_before, pks_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                neg_values_sub.append(len(pks_x_neg_before))
                                neg_label_before = True
                                warnings.warn(f"Negative Peak(s) in time period {start}-{end}", Warning)

                        if show_legend:
                                ax.legend(loc="upper left", ncol= 2)
                        if peak:
                                textstr = f"Peak Comparison\n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                        else:
                                textstr = f"Trough Comparison\n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                        ax.set_title(textstr)
                        ax.xaxis.set_minor_locator(mdates.YearLocator())
                        ax.grid(axis="x", which="minor", linewidth=0.5)
                        ax.grid(axis="x", which="major", linewidth=0.5)
                        ax.grid(axis="y", linewidth=0.5)
                        ax.set_ylabel("[ug/L]")

                        ax.set_xlim(pd.to_datetime('01-01-' + str(function_start), format='%d-%m-%Y') , pd.to_datetime('31-12-' + str(function_end), format='%d-%m-%Y'))
                        pks_lim_sub = sorted(y_sub)
                        if max(pks_lim_sub)> 10:
                                ymax = pks_lim_sub[-2]+0.5
                                ax.set_ylim(-0.5, ymax)
                        else:
                                ymax = pks_lim_sub[-1]+0.5
                                ax.set_ylim(-0.5, ymax)
                        return ymax
                else:
                        warnings.warn("No data to plot (check valid indices)")
                        return None
                        




        def extrema_comparison(self, other1,  latitude, longitude, ax,  peak = True, aggregation= False, start = 0, end = 9999, background_pts = False, other2= None, purple_chla21= False, show_legend= False):
                """Overlay extrema plots from two or three PhenologyVisualization instances on one axis.

                Calls extrema_plot for self and other1 (and optionally other2), sharing the
                same axes so that peaks or troughs from different products (e.g. chla v2.1
                vs v3.1) can be compared directly. All instances must reference the same lake.

                Parameters
                ----------
                other1 : PhenologyVisualization
                    Second instance to overlay (must share the same lake ID as self).
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                ax : matplotlib.axes.Axes
                    Axes on which to draw all overlaid plots.
                peak : bool, optional
                    If True (default), compare peaks. If False, compare troughs.
                aggregation : bool, optional
                    If True, use the 3×3 neighbourhood median as background scatter.
                start : int, optional
                    First year to display (inclusive). 0 = earliest.
                end : int, optional
                    Last year to display (inclusive). 9999 = latest.
                background_pts : bool, optional
                    If True, show background scatter for each overlay. Default False.
                other2 : PhenologyVisualization, optional
                    Optional third instance to overlay. Must share the same lake ID.
                purple_chla21 : bool, optional
                    Colour-override flag forwarded to each extrema_plot call.

                Returns
                -------
                None

                Raises
                ------
                Warning
                    If self, other1, or other2 do not all reference the same lake ID.
                """


                g = self._load_extracted_globals()
                lat = g["lat"]
                lon = g["lon"]
                lakeID1 = os.path.basename(self.p_path)[:-3]
                lakeID2 = os.path.basename(other1.p_path)[:-3]
                if other2:
                        lakeID3 = os.path.basename(other2.p_path)[:-3]

                if lakeID1 != lakeID2:
                        raise Warning("Comparison must be made on the same lake!")
                if other2:
                        if lakeID2!= lakeID3:
                                raise Warning("Comparison must be made on the same lake!")
                        ymax1 = self.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
                        ymax2 = other1.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
                        ymax3 = other2.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
                        y_lims = [ymax1, ymax2, ymax3]
                        phenology_name1 = os.path.basename(os.path.dirname(self.p_path))
                        phenology_name2 = os.path.basename(os.path.dirname(other1.p_path))
                        phenology_name3 = os.path.basename(os.path.dirname(other2.p_path))
                        
                        if phenology_name1 == "phycocyanin":
                                phenology_name1 = "phyco"
                        elif phenology_name1 == "chla_mean":
                                phenology_name1 = "chla v2.1"
                        else:
                                phenology_name1 = "chla v3.1"

                        if phenology_name2 == "phycocyanin":
                                phenology_name2 = "phyco"
                        elif phenology_name2 == "chla_mean":
                                phenology_name2 = "chla v2.1"
                        else:
                                phenology_name2 = "chla v3.1"

                        if phenology_name3 == "phycocyanin":
                                phenology_name3 = "phyco"
                        elif phenology_name3 == "chla_mean":
                                phenology_name3 = "chla v2.1"
                        else:
                                phenology_name3 = "chla v3.1"


                        if peak:
                                textr =  f"{phenology_name1}, {phenology_name2} vs {phenology_name3} Peaks \n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                        else:
                                textr =  f"{phenology_name1}, {phenology_name2} vs {phenology_name3} Troughs \n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                        ax.set_title(textr)
                        ax.set_ylim(top = max(y_lims))
                        ax.xaxis.set_minor_locator(mdates.YearLocator())
                        ax.grid(axis="x", which="minor", linewidth=0.5)
                        ax.grid(axis="x", which="major", linewidth=0.5)
                        ax.grid(axis="y", linewidth=0.5)
                        ax.set_ylabel("[ug/L]")
                        ax.legend(loc="upper left", ncol= 2)

                else:

                        ymax1 = self.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak= peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
                        ymax2 = other1.extrema_plot(latitude=latitude, longitude=longitude, ax = ax,  peak= peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
                        y_lims = [ymax1, ymax2]
                        phenology_name1 = os.path.basename(os.path.dirname(self.p_path))
                        phenology_name2 = os.path.basename(os.path.dirname(other1.p_path))

                        if phenology_name1 == "phycocyanin":
                                phenology_name1 = "phyco"
                        elif phenology_name1 == "chla_mean":
                                phenology_name1 = "chla v2.1"
                        else:
                                phenology_name1 = "chla v3.1"

                        if phenology_name2 == "phycocyanin":
                                phenology_name2 = "phyco"
                        elif phenology_name2 == "chla_mean":
                                phenology_name2 = "chla v2.1"
                        else:
                                phenology_name2 = "chla v3.1"
                        if peak: 
                                textr =  f"{phenology_name1} vs {phenology_name2} Peaks \n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                        else:
                                textr =  f"{phenology_name1} vs {phenology_name2} Troughs \n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                        
                        ax.set_title(textr)
                        ax.set_ylim(top = max(y_lims))
                        ax.xaxis.set_minor_locator(mdates.YearLocator())
                        ax.grid(axis="x", which="minor", linewidth=0.5)
                        ax.grid(axis="x", which="major", linewidth=0.5)
                        ax.grid(axis="y", linewidth=0.5)
                        ax.set_ylabel("[ug/L]")
                        ax.legend(loc="upper left", ncol= 2)





        def single_plot_background(self, latitude, longitude, ax, fig, aggregation = False, start= 0, end= 9999):
                """Plot a pixel time series with QA-coloured scatter, spline, and phenological events.

                Like single_plot, but colours each background scatter point by its QA flag
                using a discrete colormap and adds a QA colorbar to the figure. Peaks,
                troughs, green-up and green-down midpoints are overlaid as scatter markers.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the plot.
                fig : matplotlib.figure.Figure
                    Figure used to attach the QA colorbar.
                aggregation : bool, optional
                    If True, replace raw scatter with the 3×3 neighbourhood median.
                    Requires spatial_aggregation() to have been called or available cache.
                start : int, optional
                    First year to display (inclusive). 0 = earliest in the series.
                end : int, optional
                    Last year to display (inclusive). 9999 = latest in the series.

                Returns
                -------
                None
                """

                with netCDF4.Dataset(self.e_path) as nc:
                        # change np.array to np.asarry to avoid getting DeprecationError
                        lat = np.asarray(nc.variables["lat"])
                        lon = np.asarray(nc.variables["lon"])
                        t_all = unix_to_datenum(nc.variables["time"])


                with netCDF4.Dataset(self.p_path) as nc:
                        smoothing = float(nc.variables["smoothing_parameter"][latitude, longitude])
                        pks_x = unix_to_datetime(remove_nan(nc.variables["pks_x"][latitude, longitude, :]))
                        pks_y = remove_nan(nc.variables["pks_y"][latitude, longitude, :])
                        mask_pks = np.array([(d.year <= end) & (d.year >= start) for d in pks_x])
                        pks_x_sub = pks_x[mask_pks]
                        pks_y_sub = pks_y[mask_pks]

                        trgs_x = unix_to_datetime(remove_nan(nc.variables["trgs_x"][latitude, longitude, :]))
                        trgs_y = remove_nan(nc.variables["trgs_y"][latitude, longitude, :])
                        mask_trgs = np.array([(d.year <= end) & (d.year >= start) for d in trgs_x])
                        trgs_x_sub = trgs_x[mask_trgs]
                        trgs_y_sub = trgs_y[mask_trgs]

                        midUP_x = unix_to_datetime(remove_nan(nc.variables["green_up_mid_x"][latitude, longitude, :]))
                        midUP_y = remove_nan(nc.variables["green_up_mid_y"][latitude, longitude, :])
                        mask_midUP = np.array([(d.year <= end) & (d.year >= start) for d in midUP_x])
                        midUP_x_sub = midUP_x[mask_midUP]
                        midUP_y_sub = midUP_y[mask_midUP]

                        midDOWN_x = unix_to_datetime(remove_nan(nc.variables["green_down_mid_x"][latitude, longitude, :]))
                        midDOWN_y = remove_nan(nc.variables["green_down_mid_y"][latitude, longitude, :])
                        mask_midDOWN= np.array([(d.year <= end) & (d.year >= start) for d in midDOWN_x])
                        midDOWN_x_sub = midDOWN_x[mask_midDOWN]
                        midDOWN_y_sub = midDOWN_y[mask_midDOWN]

                with netCDF4.Dataset(self.e_path) as nc:
                        variable = getattr(nc, "variable")
                        variable_qa  = getattr(nc, "qa")
                        values = np.array(nc.variables[variable][:, latitude, longitude])
                        values_qa = np.array(nc.variables[variable_qa][:, latitude, longitude])
                        mask = (values != -9999)
                        values_m = values[mask]
                        time_m = t_all[mask]
                        qa_mask = values_qa[mask]
                        qa_unique = sorted(np.unique(qa_mask))

                if len(values_m) > 1:

                        limits = sorted(datenum_to_datetime(time_m))
                        if start == 0:
                                function_start = min(limits).year
                        else:
                                function_start = start
                        if end == 9999:
                                function_end= max(limits).year
                        else:
                                function_end = end
                        neg_values_sub =[]
                        smooth_x = np.arange(t_all.min(), t_all.max() + 1, 1)
                        smooth_y = csaps(time_m, values_m, smooth_x, smooth=smoothing)

                        y_pred =csaps(time_m, values_m, time_m, smooth=smoothing)
                        y_true = values_m

                        time_slice = np.array(datenum_to_datetime(time_m))

                        mask_sub = np.array([(d.year <=function_end) & (d.year >= function_start) for d in time_slice])
                        if mask_sub.sum()>2:
                                rmse_sub = np.sqrt(mean_squared_error(y_true[mask_sub], y_pred[mask_sub]))
                                r2_sub = r2_score(y_true[mask_sub], y_pred[mask_sub])

                                rmse_tot = np.sqrt(mean_squared_error(y_true, y_pred))
                                r2_tot = r2_score(y_true, y_pred)

                                # data1 = {"time": time_m,"data": values_m}
                                # average_df1 = pd.DataFrame(data1)

                                # data2 = {"data": values_m}
                                # average_df2 = pd.DataFrame(data2)

                                # average_df2["SMA15"] = average_df2.rolling(window = 15).mean()

                                # average_background = pd.merge(average_df1, average_df2, on='data', how='inner')

                                neg_label_before = False

                                qa_unique = sorted(np.unique(qa_mask))

                                # one discrete color per unique QA value
                                cmap = plt.cm.get_cmap("tab10", len(qa_unique))
                                cmap_new = ListedColormap(cmap(np.arange(len(qa_unique))))

                                # map actual QA values to category indices 0..N-1
                                qa_to_idx = {qa: i for i, qa in enumerate(qa_unique)}
                                qa_idx = np.array([qa_to_idx[q] for q in qa_mask])

                                # boundaries around category indices
                                bounds = np.arange(-0.5, len(qa_unique) + 0.5, 1)
                                norm = BoundaryNorm(bounds, cmap_new.N)

                                if aggregation:
                                        if self.aggregation_df is None:
                                                self.spatial_aggregation()

                                        background_sub = self.aggregation_df[(self.aggregation_df["i"]==latitude) & (self.aggregation_df["j"]==longitude)]
                                        background_time = background_sub["time"]
                                        background_time = background_time.to_numpy()

                                        background_values = background_sub["MA_value"]

                                        sc = ax.scatter(datenum_to_datetime(background_time), background_values, c=qa_idx, cmap = cmap_new, norm = norm, alpha=1, s=10, label="Data")
                                else:
                                        sc = ax.scatter(datenum_to_datetime(time_m), values_m, c=qa_idx, cmap = cmap_new, norm = norm, alpha=1, s=10, label="Data")
                                ax.plot(datenum_to_datetime(smooth_x), smooth_y, color="black", linewidth=1, label="Spline")
                                ax.scatter(pks_x_sub, pks_y_sub, color="orange", s=50, marker="o", zorder=4, label="Peaks")
                                if (pks_y_sub < 0).any():
                                        mask =  pks_y_sub<0
                                        pks_x_neg_before = pks_x_sub[mask]
                                        pks_y_neg_before = pks_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(pks_x_neg_before, pks_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(pks_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Peak(s) in time period {start}-{end}", Warning)

                                ax.scatter(trgs_x_sub, trgs_y_sub, color="blue", s=50, marker="o", zorder=4, label="Troughs")
                                if (trgs_y_sub < 0).any():
                                        mask =  trgs_y_sub<0
                                        trgs_x_neg_before = trgs_x_sub[mask]
                                        trgs_y_neg_before = trgs_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(trgs_x_neg_before, trgs_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(trgs_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Troughs(s) in time period {start}-{end}", Warning)
                                ax.scatter(midUP_x_sub, midUP_y_sub, color="mediumseagreen", s=30, marker="^", zorder=4, label="Mid Up")
                                if (midUP_y_sub < 0).any():
                                        mask =  midUP_y_sub<0
                                        midUP_x_neg_before = midUP_x_sub[mask]
                                        midUP_y_neg_before = midUP_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(midUP_x_neg_before, midUP_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(midUP_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Mid Up(s) in time period {start}-{end}", Warning)
                                ax.scatter(midDOWN_x_sub, midDOWN_y_sub, color="darkgreen", s=30, marker="v", zorder=4, label="Mid Down")
                                if (midDOWN_y_sub < 0).any():
                                        mask =  midDOWN_y_sub<0
                                        midDOWN_x_neg_before = midDOWN_x_sub[mask]
                                        midDOWN_y_neg_before = midDOWN_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(midDOWN_x_neg_before, midDOWN_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(midDOWN_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Mid Down(s) in time period {start}-{end}", Warning)

                                ax.legend(loc="upper left", ncol= 2)
                                textstr = f"{os.path.basename(os.path.dirname(self.p_path))}\n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}\n Total RMSE, R$^2$: {round(rmse_tot,4)}, {round(r2_tot, 4)}"
                                ax.set_title(textstr)
                                ax.grid()
                                ax.set_ylabel("[ug/L]")
                                ax.set_xlim(pd.to_datetime('01-01-' + str(function_start), format='%d-%m-%Y') , pd.to_datetime('31-12-' + str(function_end), format='%d-%m-%Y'))
                                cbar = fig.colorbar(sc, ax=ax, boundaries=bounds)
                                cbar.set_label("QA indicators")

                                # ticks at category centers
                                cbar.set_ticks(np.arange(len(qa_unique)))
                                cbar.set_ticklabels([str(q) for q in qa_unique])


                                pks_lim_sub = sorted(pks_y_sub)
                                if max(pks_lim_sub)> 10:
                                        ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-2]+0.5)
                                else:
                                        ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-1]+0.5)

                                if neg_values_sub:
                                        ax.text(0.99,0.99,f"# Neg.values: {sum(neg_values_sub)} \n RMSE: {round(rmse_sub,4)}\n R$^2$: {round(r2_sub,4)}", transform = ax.transAxes,   ha= "right", va= "top")
                                else:
                                        ax.text(0.99,0.99,f"RMSE:{round(rmse_sub,4)} \n R$^2$: {round(r2_sub,4)}", transform = ax.transAxes,   ha= "right", va= "top")
                        else:
                                warnings.warn("Not enough data to plot or compute metrics for chosen time interval")
                else:
                        warnings.warn("No data to plot")



        def count_peaks(self, latitude, longitude, start= 0, end= 9999):
                """Return the number of detected peaks for a pixel within a year range.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                start : int, optional
                    First year to include (inclusive). 0 = earliest in the series.
                end : int, optional
                    Last year to include (inclusive). 9999 = latest in the series.

                Returns
                -------
                int
                    Number of peaks falling within the specified year range.
                """
                with netCDF4.Dataset(self.p_path) as nc:
                        pks_x = unix_to_datetime(remove_nan(nc.variables["pks_x"][latitude, longitude, :]))
                        mask_pks = np.array([(d.year <= end) & (d.year >= start) for d in pks_x])
                        pks_x_sub = pks_x[mask_pks]

                        return len(pks_x_sub)
                
        def pixel_r2(self, latitude, longitude, start, end):
                """Return the R² score for a single pixel within a year range.

                Delegates to r2_scores, which may trigger full-lake parallel computation
                and CSV caching on the first call for this year range.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                float
                    R² score for the pixel, or np.nan if insufficient data.
                """
                scores = self.r2_scores([(start, end)])
                return scores[(latitude, longitude)]

        def pixel_rmse(self, latitude, longitude, start, end):
                """Return the RMSE for a single pixel within a year range.

                Delegates to RMSE_scores, which may trigger full-lake parallel computation
                and CSV caching on the first call for this year range.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                float
                    RMSE value for the pixel, or np.nan if insufficient data.
                """
                scores = self.RMSE_scores([(start, end)])
                return scores[(latitude, longitude)]

        def pixel_mad(self, latitude, longitude, start, end):
                """Return the Median Absolute Deviation for a single pixel within a year range.

                Delegates to MAD_scores, which may trigger full-lake parallel computation
                and CSV caching on the first call for this year range.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                float
                    MAD value for the pixel, or np.nan if insufficient data.
                """
                scores = self.MAD_scores([(start, end)])
                return scores[(latitude, longitude)]

        def pixel_correlation(self, latitude, longitude, start, end):
                """Return the Pearson correlation coefficient for a single pixel within a year range.

                Delegates to correlation_scores, which may trigger full-lake parallel
                computation and CSV caching on the first call for this year range.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                float
                    Pearson r for the pixel, or np.nan if insufficient data.
                """
                scores = self.correlation_scores([(start, end)])
                return scores[(latitude, longitude)]

        def pixel_values(self, latitude, longitude, start, end):
                """Return the valid observation count for a single pixel within a year range.

                Delegates to values_per_pixel, which may trigger full-lake parallel
                computation and CSV caching on the first call for this year range.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                start : int
                    First year of the evaluation window (0 = full series start).
                end : int
                    Last year of the evaluation window (9999 = full series end).

                Returns
                -------
                int
                    Number of valid (QA==0, value!=-9999) observations in the window.
                """
                scores = self.values_per_pixel([(start, end)])
                return scores[(latitude, longitude)]



        def single_plot(self, latitude, longitude, ax, aggregation = False, start= 0, end= 9999, annotation = None):
                """Plot raw observations, the smoothed spline, and all phenological events for a pixel.

                Displays a scatter of valid (QA==0) observations or 3×3 aggregated values,
                overlaid with the csaps spline and scatter markers for peaks, troughs,
                green-up midpoints, and green-down midpoints. Negative values are flagged
                with red crosses and trigger a warning. Fit metrics (RMSE, R², MAD) are
                annotated on the plot.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the plot.
                aggregation : bool, optional
                    If True, replace raw scatter with the 3×3 neighbourhood spatial median.
                    Requires spatial_aggregation() to have been called or available cache.
                start : int, optional
                    First year to display (inclusive). 0 = earliest in the series.
                end : int, optional
                    Last year to display (inclusive). 9999 = latest in the series.

                Returns
                -------
                None
                """


                g  = self._load_extracted_globals()
                px = self._load_pixel_data(latitude, longitude)
                lat, lon, t_all = g["lat"], g["lon"], g["t_all"]

                smoothing = px["smoothing"]

                mask_pks     = np.array([(d.year <= end) & (d.year >= start) for d in px["pks_x"]])
                pks_x_sub    = px["pks_x"][mask_pks]
                pks_y_sub    = px["pks_y"][mask_pks]

                mask_trgs    = np.array([(d.year <= end) & (d.year >= start) for d in px["trgs_x"]])
                trgs_x_sub   = px["trgs_x"][mask_trgs]
                trgs_y_sub   = px["trgs_y"][mask_trgs]

                mask_midUP   = np.array([(d.year <= end) & (d.year >= start) for d in px["midUP_x"]])
                midUP_x_sub  = px["midUP_x"][mask_midUP]
                midUP_y_sub  = px["midUP_y"][mask_midUP]

                mask_midDOWN  = np.array([(d.year <= end) & (d.year >= start) for d in px["midDOWN_x"]])
                midDOWN_x_sub = px["midDOWN_x"][mask_midDOWN]
                midDOWN_y_sub = px["midDOWN_y"][mask_midDOWN]

                mask     = (px["values"] != -9999) & (px["qa"] == 0)
                values_m = px["values"][mask]
                time_m   = t_all[mask]

                if len(values_m) > 1:

                        limits = sorted(datenum_to_datetime(time_m))
                        if start == 0:
                                function_start = min(limits).year
                        else:
                                function_start = start
                        if end == 9999:
                                function_end= max(limits).year
                        else:
                                function_end = end
                        neg_values_sub =[]
                        smooth_x = np.arange(t_all.min(), t_all.max() + 1, 1)
                        smooth_y = csaps(time_m, values_m, smooth_x, smooth=smoothing)

                        y_pred =csaps(time_m, values_m, time_m, smooth=smoothing)
                        y_true = values_m

                        time_slice = np.array(datenum_to_datetime(time_m))

                        mask_sub = np.array([(d.year <=function_end) & (d.year >= function_start) for d in time_slice])
                        if mask_sub.sum()>2:
                                rmse_sub = np.sqrt(mean_squared_error(y_true[mask_sub], y_pred[mask_sub]))
                                r2_sub = r2_score(y_true[mask_sub], y_pred[mask_sub])
                                mad_sub = np.median(np.abs(y_true[mask_sub]-y_pred[mask_sub]))

                                rmse_tot = np.sqrt(mean_squared_error(y_true, y_pred))
                                r2_tot = r2_score(y_true, y_pred)
                                mad_tot = np.median(np.abs(y_true-y_pred))



                                neg_label_before = False

                                if aggregation:
                                        if self.aggregation_df is None:
                                                self.spatial_aggregation()

                                        background_sub = self.aggregation_df[(self.aggregation_df["i"]==latitude) & (self.aggregation_df["j"]==longitude)]

                                        background_time = background_sub["time"]
                                        background_time = background_time.to_numpy()

                                        background_values = background_sub["MA_value"]

                                        ax.scatter(datenum_to_datetime(background_time), background_values, color="grey", alpha=0.3, s=10, label="Data")
                                else:
                                        ax.scatter(datenum_to_datetime(time_m), values_m, color="grey", alpha=0.3, s=10, label="Data")
                                ax.plot(datenum_to_datetime(smooth_x), smooth_y, color="black", linewidth=1, label="Spline")
                                ax.scatter(pks_x_sub, pks_y_sub, color="orange", s=50, marker="o", zorder=4, label="Peaks")
                                if (pks_y_sub < 0).any():
                                        mask =  pks_y_sub<0
                                        pks_x_neg_before = pks_x_sub[mask]
                                        pks_y_neg_before = pks_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(pks_x_neg_before, pks_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(pks_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Peak(s) in time period {start}-{end}", Warning)

                                ax.scatter(trgs_x_sub, trgs_y_sub, color="blue", s=50, marker="o", zorder=4, label="Troughs")
                                if (trgs_y_sub < 0).any():
                                        mask =  trgs_y_sub<0
                                        trgs_x_neg_before = trgs_x_sub[mask]
                                        trgs_y_neg_before = trgs_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(trgs_x_neg_before, trgs_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(trgs_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Troughs(s) in time period {start}-{end}", Warning)

                                ax.scatter(midUP_x_sub, midUP_y_sub, color="mediumseagreen", s=30, marker="^", zorder=4, label="Mid Up")
                                if (midUP_y_sub < 0).any():
                                        mask =  midUP_y_sub<0
                                        midUP_x_neg_before = midUP_x_sub[mask]
                                        midUP_y_neg_before = midUP_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(midUP_x_neg_before, midUP_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(midUP_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Mid Up(s) in time period {start}-{end}", Warning)

                                ax.scatter(midDOWN_x_sub, midDOWN_y_sub, color="darkgreen", s=30, marker="v", zorder=4, label="Mid Down")
                                if (midDOWN_y_sub < 0).any():
                                        mask =  midDOWN_y_sub<0
                                        midDOWN_x_neg_before = midDOWN_x_sub[mask]
                                        midDOWN_y_neg_before = midDOWN_y_sub[mask]
                                        label = "Negative Value" if not neg_label_before else None
                                        ax.scatter(midDOWN_x_neg_before, midDOWN_y_neg_before, color="red", s=50, marker="x", zorder=6, label=label)
                                        neg_values_sub.append(len(midDOWN_x_neg_before))
                                        neg_label_before = True
                                        warnings.warn(f"Negative Mid Down(s) in time period {start}-{end}", Warning)

                                ax.legend(loc="upper left", ncol= 2)
                                textstr = f"{os.path.basename(os.path.dirname(self.p_path))}\n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}\n Total RMSE, R$^2$, MAD: {round(rmse_tot,4)}, {round(r2_tot, 4)}, {round(mad_tot,4)}"
                                ax.set_title(textstr)
                                ax.xaxis.set_minor_locator(mdates.YearLocator())
                                ax.grid(axis="x", which="minor", linewidth=0.5)
                                ax.grid(axis="x", which="major", linewidth=0.5)
                                ax.grid(axis="y")
                                ax.set_ylabel("[ug/L]")
                                ax.set_xlim(pd.to_datetime('01-01-' + str(function_start), format='%d-%m-%Y') , pd.to_datetime('31-12-' + str(function_end), format='%d-%m-%Y'))


                                pks_lim_sub = sorted(pks_y_sub)
                                # if max(pks_lim_sub)> 10:
                                #         ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-2]+0.5)
                                # else:
                                #         ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-1]+0.5)



                                trgs_lim_sub = sorted(trgs_y_sub)

                                if len(pks_lim_sub) > 0 and len(trgs_lim_sub) > 0:

                                        ymax = pks_lim_sub[-1]

                                        if ymax > 10 and len(pks_lim_sub) > 1:
                                                ymax = pks_lim_sub[-2]

                                        ymin = trgs_lim_sub[0]

                                        ax.set_ylim(ymin - 0.5, ymax + 0.5)

                                else:
                                        warnings.warn(
                                                f"No peaks/troughs available for {start}-{end}; using automatic y-limits."
                                        )
                                if not annotation:
                                        if neg_values_sub:
                                                ax.text(0.99,0.99,f"# Neg.values: {sum(neg_values_sub)} \n RMSE: {round(rmse_sub,3)}\n R$^2$: {round(r2_sub,3)}\n MAD: {round(mad_sub, 3)}", transform = ax.transAxes,   ha= "right", va= "top", zorder = 10)
                                        else:
                                                ax.text(0.99,0.99,f"RMSE:{round(rmse_sub,3)} \n R$^2$: {round(r2_sub,3)}\n MAD: {round(mad_sub, 3)}", transform = ax.transAxes,   ha= "right", va= "top", zorder = 10)
                                else:
                             
                                        lines = []
                                        if "R2" in annotation:
                                                lines.append(f"R$^2$: {round(r2_sub, 3)}")
                                        if "RMSE" in annotation:
                                                lines.append(f"RMSE: {round(rmse_sub, 3)}")
                                        if "MAD" in annotation:
                                                lines.append(f"MAD: {round(mad_sub, 3)}")
                                        if "neg" in annotation and neg_values_sub:
                                                lines.append(f"# Neg.values: {sum(neg_values_sub)}")
                                        
                                        if lines:
                                                ax.text(0.99, 0.99, "\n".join(lines),
                                                        transform=ax.transAxes, ha="right", va="top", zorder = 10)
                                        
                        else:
                                warnings.warn("Not enough data to plot or compute metrics for chosen time interval")
                else:
                        warnings.warn("No data to plot")


        def split_plot(self, latitude, longitude, ax0, ax1, aggregation = False, start0= 0, end0= 9999, start1= 0, end1=9999):
                """Plot two year-windowed single_plots side by side for the same pixel.

                Calls single_plot twice — once on ax0 with [start0, end0] and once on ax1
                with [start1, end1]. Intended for comparing two non-overlapping time periods
                at the same pixel.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                ax0 : matplotlib.axes.Axes
                    Axes for the first time window.
                ax1 : matplotlib.axes.Axes
                    Axes for the second time window.
                aggregation : bool, optional
                    If True, use 3×3 spatial median scatter on both panels.
                start0 : int, optional
                    First year of the first window. 0 = earliest in the series.
                end0 : int, optional
                    Last year of the first window. 9999 = latest in the series.
                start1 : int, optional
                    First year of the second window. 0 = earliest in the series.
                end1 : int, optional
                    Last year of the second window. 9999 = latest in the series.

                Returns
                -------
                None
                """
                if start0 and start1 == 0 and end0 and end1 == 9999:
                        warnings.warn("split_plot needs a least end0 and start1 parameter, otherwise use full_plot")

                self.single_plot(latitude = latitude, longitude= longitude, ax=ax0, aggregation = aggregation, start= start0, end = end0)
                self.single_plot(latitude = latitude, longitude= longitude, ax=ax1, aggregation = aggregation, start=start1, end=end1)

        def full_plot(self, latitude, longitude, ax, aggregation = False):
                """Plot the complete valid time series for a pixel.

                Auto-detects the first and last years with valid observations and
                delegates to single_plot with those bounds.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                ax : matplotlib.axes.Axes
                    Axes on which to draw the plot.
                aggregation : bool, optional
                    If True, use 3×3 neighbourhood spatial median scatter.

                Returns
                -------
                None
                """
                with netCDF4.Dataset(self.e_path) as nc:
                        t_all = unix_to_datenum(nc.variables["time"])
                        variable = getattr(nc, "variable")
                        values = np.array(nc.variables[variable][:, latitude, longitude])
                        mask = (values != -9999) & (np.array(nc.variables[getattr(nc, 'qa')][:, latitude, longitude]) == 0)
                        values_m = values[mask]
                        time_m = t_all[mask]

                if len(values_m) > 1:
                        limits = sorted(datenum_to_datetime(time_m))
                        full_plot_start = min(limits).year
                        full_plot_end = max(limits).year
                        self.single_plot(latitude= latitude, longitude= longitude, ax = ax, aggregation = aggregation, start=full_plot_start, end= full_plot_end)

                else:
                        warnings.warn("No data to plot")


        def single_years_plot(self, latitude, longitude, years, ncol, nrow, annotation, ylim=None):
                """Plot one panel per year in a grid, each showing phenology for a single pixel.

                Creates a figure with ``nrow × ncol`` subplots. Each subplot calls
                :meth:`single_plot` for one year in *years*, with enlarged markers and
                month-number x-axis labels. Panels beyond ``len(years)`` are hidden.

                Parameters
                ----------
                latitude : int
                    Row (lat) index of the pixel.
                longitude : int
                    Column (lon) index of the pixel.
                years : list of int
                    Calendar years to display, one per panel.
                ncol : int
                    Number of subplot columns.
                nrow : int
                    Number of subplot rows. Must satisfy ``ncol * nrow >= len(years)``.
                annotation : list of str or None
                    Passed to :meth:`single_plot`. Controls which fit metrics are shown
                    (e.g. ``["R2", "RMSE", "MAD"]``). ``None`` shows all metrics.
                ylim : tuple of (float, float) or None, optional
                    If provided, sets the y-axis limits as ``(bottom, top)`` for every
                    panel, overriding the automatic limits set by :meth:`single_plot`.

                Returns
                -------
                matplotlib.figure.Figure
                """
                _MARKER_SIZES = {"Data": 50, "Peaks": 150, "Troughs": 150, "Mid Up": 150, "Mid Down": 150}

                fig, axs = plt.subplots(nrow, ncol, constrained_layout=True, squeeze=False, figsize=(ncol * 5, nrow * 4))
                for year, ax in zip(years, axs.flatten()):
                        self.single_plot(latitude, longitude, ax, start=year, end=year, annotation=annotation)

                        for col in ax.collections:
                                if col.get_label() in _MARKER_SIZES:
                                        col.set_sizes([_MARKER_SIZES[col.get_label()]])

                        if ax.texts:
                                ax.texts[-1].set_fontsize(15)

                        legend = ax.get_legend()
                        if legend is not None:
                                legend.remove()
                        if ylim is not None:
                                ax.set_ylim(ylim)
                        else:
                                ax.set_ylim(bottom=-0.5,top = ax.get_ylim()[1]*1.5)
                        ax.set_title(str(year), fontsize=20)
                        ax.set_ylabel("[ug/L]", fontsize=15)
                        ax.xaxis.set_major_locator(mdates.MonthLocator())
                        ax.xaxis.set_major_formatter(mdates.DateFormatter('%#m'))
                        ax.tick_params(labelsize=15)

                for ax in axs.flatten()[len(years):]:
                        ax.set_visible(False)

        def single_plot_insitu(self, latitude,longitude,ax,insitu_df,aggregation=False, start = 0, end = 9999, insitu_date_col="datetime", insitu_value_col="chlorophyll_a", insitu_station_col=None, station_id=None, max_depth = 5):
                """
                Plot satellite observations + spline + phenology + in situ overlay.
                """

                # -------------------------------------------------
                # FIRST: draw the original plot
                # -------------------------------------------------

                self.single_plot(
                        latitude=latitude,
                        longitude=longitude,
                        ax=ax,
                        aggregation=aggregation,
                        start=start,
                        end=end)

                # -------------------------------------------------
                # PREPARE IN SITU DATA
                # -------------------------------------------------

                insitu_mean = prep_dimark_data(insitu_df = insitu_df,start=start, end=end,  insitu_date_col=insitu_date_col, insitu_value_col=insitu_value_col, insitu_station_col=insitu_station_col, station_id=station_id, max_depth = max_depth)
                
                # -------------------------------------------------
                # OVERLAY IN SITU DATA
                # -------------------------------------------------

                ax.scatter(
                        insitu_mean[insitu_date_col],
                        insitu_mean[insitu_value_col],
                        color="red",
                        marker="D",
                        s=20,
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=5,
                        alpha= 0.5,
                        label=f"In Situ (<{max_depth}m)"
                )

                # optional connecting line
                # ax.plot(
                #         insitu[insitu_date_col],
                #         insitu[insitu_value_col],
                #         color="red",
                #         alpha=0.5,
                #         linewidth=1
                # )

                # -------------------------------------------------
                # UPDATE LEGEND
                # -------------------------------------------------

                handles, labels = ax.get_legend_handles_labels()

                by_label = dict(zip(labels, handles))

                ax.legend(
                        by_label.values(),
                        by_label.keys(),
                        loc="upper left",
                        ncol=2
                )


















