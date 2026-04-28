import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
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
from functions import unix_to_datetime, unix_to_datenum, datenum_to_datetime, remove_nan, define_year_range
import multiprocessing
from functools import partial
import shapely.ops as ops
from pyproj import CRS, Transformer
import geopandas
from shapely.geometry import Point
from shapely.prepared import prep
from numpy.lib.stride_tricks import sliding_window_view



_GLOBALS = {}


def _init_worker(p_path, e_path):
    """
    Initialize global data for a worker process by loading and preprocessing
    NetCDF datasets.

    This function is intended to be used in a multiprocessing context. It loads
    two NetCDF files (typically parameter and evaluation datasets), extracts
    relevant variables, performs one-time preprocessing (such as time conversion),
    and stores everything in a shared global dictionary (`_GLOBALS`) for fast
    access by worker tasks.

    Parameters
    ----------
    p_path : str
        File path to the NetCDF dataset containing parameter data
        (e.g., smoothing parameters).
    e_path : str
        File path to the NetCDF dataset containing evaluation data
        (e.g., variable values, QA flags, time, lat/lon).

    Notes
    -----
    - Time values are converted from Unix timestamps to datenums and then to
      datetime objects. The corresponding year for each timestamp is extracted
      and cached.
    - All large arrays (values, QA flags, smoothing parameters, etc.) are loaded
      into memory once per worker to avoid repeated I/O overhead.
    - The function relies on attributes `variable` and `qa` being present in
      the evaluation dataset (`nc_e`).
    - Results are stored in the module-level `_GLOBALS` dictionary for reuse
      across worker calls.

    Globals Set
    -----------
    _GLOBALS : dict
        Dictionary populated with the following keys:
        - "nc_p": NetCDF dataset for parameter data
        - "nc_e": NetCDF dataset for evaluation data
        - "variable": Name of the main variable
        - "qa": Name of the QA variable
        - "t_all": Time values in datenumber format
        - "years_all": Extracted year for each time step
        - "smoothing_all": Smoothing parameter array
        - "values_all": Main variable data array
        - "qa_all": QA flag array
        - "lats": Latitude array
        - "lons": Longitude array
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


class PhenologyEDA:
        shapefile_path = None

        def __init__(self, extract_path, phenology_path):
                if self.shapefile_path is None:
                        raise Warning("Please define your path to the lake CCI shapefile. This can be done at a class level using PhenologyEDA.set_shapefile_path(your_path)")
                else:
                        self.geom = geopandas.read_file(self.shapefile_path)
                self.p_path = phenology_path
                self.e_path = extract_path
                self.info = (f"Version: {os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(phenology_path))))[1:]} \n",
                        f" Variable: {os.path.basename(os.path.dirname(os.path.dirname(phenology_path)))} \n" ,
                        f"Lake ID: {os.path.basename(phenology_path)[:-3]}")
                self.version = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(phenology_path))))[1:]
                self.variable = os.path.basename(os.path.dirname(phenology_path))
                self.methods = [method for method in dir(PhenologyEDA) if callable(getattr(PhenologyEDA, method)) and not method.startswith("__")]
                self.valid_coords = self.valid_index_pairs()
                self.out_folder =  os.path.dirname(os.path.dirname(os.path.dirname(self.p_path)))
                self.aggregation_df = None
                self.geom_shrunk = None


        def shortname_to_name(self, short_name:str):
                df = self.geom
                return f"ID: {df.loc[df["short_name"]==short_name, "name"][1]}"

        def ID_to_name(self, id:int):
                df = self.geom
                return list(df.loc[df["id"]==id, "name"])[0]

        def name_to_ID(self, name:str):
                df = self.geom
                return f"ID: {list(df.loc[df["name"]==name, "id"])[0]}"

        def name_to_shortname(self, name:str):
                df = self.geom
                return f"short_name: {list(df.loc[df["name"]==name, "short_name"])[0]}"

        def index_to_lat_lon(self, lat_index, lon_index):
                with netCDF4.Dataset(self.e_path) as nc:
                        lats = nc.variables["lat"][:]
                        lons = nc.variables["lon"][:]
                        lat = lats[lat_index]
                        lon = lons[lon_index]
                return f"Lat, Lon: {lat}, {lon}"
        
        
        @classmethod
        def set_shapefile_path(cls, path: str):
                """
    Set the class-level path to the shapefile.

    This method assigns a file path to a class attribute so it can be
    accessed globally by all instances of the class.

    Parameters
    ----------
    path : str
        Path to the shapefile (e.g., .shp file) to be used.

    Notes
    -----
    - This updates a class attribute (`shapefile_path`), not an instance attribute.
    - All existing and future instances of the class will use this updated path.
    """
                cls.shapefile_path = path



        def valid_index_pairs(self):
                """
    Identify spatial indices with sufficient valid data points.

    This method reads the evaluation NetCDF dataset, extracts the target
    variable and its associated quality assurance (QA) variable, and
    determines which spatial locations contain more than one valid
    observation over time.

    A data point is considered valid if:
    - Its value is not equal to the missing value flag (-9999), and
    - Its QA flag is equal to 0.

    Returns
    -------
    list of tuple of int
        A list of (row, column) index pairs corresponding to spatial
        locations where the number of valid observations exceeds one.

    Notes
    -----
    - The method assumes the NetCDF dataset has attributes `variable`
      and `qa` that specify the variable names.
    - Validation is performed across the time dimension (axis=0).
    - The dataset is opened and closed within the method using a context manager.
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
                """
    Create a pandas DataFrame from phenology NetCDF data for a specific location.

    This method reads the phenology dataset stored at `self.p_path`, extracts
    paired variables ending in 'x' (time) and 'y' (values), and constructs
    individual DataFrames for each variable. The resulting DataFrames are then
    concatenated into a single DataFrame indexed by time.

    Processing steps:
    - Excludes non-relevant variables (e.g., 'lat', 'lon', quality flags).
    - Identifies variable pairs where:
        - 'x' variables represent time (converted from Unix timestamps).
        - 'y' variables represent observed values.
    - Removes NaN values from both time and value arrays.
    - Converts time values to datetime format.
    - Creates a labeled DataFrame for each variable.
    - Concatenates all variable DataFrames into a single DataFrame.

    Returns:
        pandas.DataFrame: A combined DataFrame containing all variables with:
            - Index: DatetimeIndex labeled "Time"
            - Columns:
                - "Value": Observed values
                - "Variable": Name of the variable

    Notes:
        - The method assumes that variables ending in 'x' and 'y' are paired
          and aligned in order.
        - Latitude and longitude are used as indices to extract point data
          from the NetCDF dataset.
        - Helper functions `remove_nan` and `unix_to_datetime` must be defined.

    Raises:
        KeyError: If expected variables are missing in the dataset.
        IndexError: If latitude/longitude indices are out of bounds.
        IOError: If the NetCDF file cannot be opened.
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
                """
    Shrink a geometry inward by 1 kilometer using a local projection.

    This method reduces the size of an input geometry by applying a negative
    buffer of 1 km. To ensure distance accuracy, the geometry is first
    projected from geographic coordinates (WGS84) into a local Azimuthal
    Equidistant (AEQD) projection centered on the geometry's centroid.
    The buffer operation is performed in meters, and the result is then
    transformed back to WGS84 coordinates.

    Parameters
    ----------
    geom : shapely.geometry.base.BaseGeometry
        Input geometry in WGS84 coordinates (longitude/latitude).

    Returns
    -------
    shapely.geometry.base.BaseGeometry
        The geometry shrunk inward by 1 km, in WGS84 coordinates.

    Notes
    -----
    - A local AEQD projection is used to preserve distance accuracy for
      the buffering operation.
    - A negative buffer value (-1000 meters) is applied to shrink the geometry.
    - If the geometry is smaller than 1 km in any dimension, the result
      may be empty or invalid.
    - The resulting geometry is also stored in `self.geom_shrunk`.
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




        def build_metric_path(self, metric_name, start= 0, end= 9999):
                if metric_name == "values_per_pixel":
                        base = os.path.join(
                        self.out_folder,
                        "calculated_values", "metrics", metric_name, os.path.basename(os.path.dirname(self.p_path)), os.path.basename(self.p_path)[:-3])
                        if start == 0 and end == 9999:
                                fname = "full_ts.csv"
                        elif start == 0:
                                fname = f"ts_end{end}.csv"
                        elif end == 9999:
                                fname = f"ts_start{start}.csv"
                        else:
                                fname = f"ts_start{start}_end{end}.csv"
                        return base, os.path.join(base, fname)

                else:
                        base = os.path.join(
                        self.out_folder,
                        "calculated_values", "metrics", metric_name,
                        os.path.basename(os.path.dirname(self.p_path)),
                        os.path.basename(self.p_path)[:-3]
                        )
                        if start == 0 and end == 9999:
                                fname = "full_ts.csv"
                        elif start == 0:
                                fname = f"ts_end_{end}.csv"
                        elif end == 9999:
                                fname = f"ts_start_{start}.csv"
                        else:
                                fname = f"ts_start_{start}_end_{end}.csv"
                        return base, os.path.join(base, fname)
                
        

        def compute_and_cache_metric(self, metric_name, col_name, compute_fn, start, end):
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


        def r2_scores(self, start=0, end=9999):
                return self.compute_and_cache_metric(metric_name="r2", col_name="r2_scores", compute_fn=PhenologyEDA.compute_metric_score, start=start, end=end)

        def MAD_scores(self, start=0, end=9999):
                return self.compute_and_cache_metric(metric_name="MAD", col_name="mad_scores",compute_fn= PhenologyEDA.compute_metric_score,start= start,end= end)

        def RMSE_scores(self, start=0, end=9999):
                return self.compute_and_cache_metric(metric_name="RMSE", col_name="rmse_scores", compute_fn=PhenologyEDA.compute_metric_score, start=start, end=end)

        def correlation_scores(self, start=0, end=9999):
                return self.compute_and_cache_metric(metric_name="correlation", col_name="correlation_scores", compute_fn=PhenologyEDA.compute_metric_score, start=start, end=end)
        
        def values_per_pixel(self, start=0, end=9999):
                return self.compute_and_cache_metric(metric_name="values_per_pixel", col_name="number_of_values",compute_fn= PhenologyEDA.compute_metric_score, start=start,end= end)





        def spatial_aggregation(self):
                """
                Load or compute spatially aggregated background values for valid pixels.
                """

                out_dir = os.path.join(
                        self.out_folder,
                        "calculated_values/spatial_aggregation_values",
                        os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(self.p_path)))),
                        os.path.basename(os.path.dirname(self.p_path)),
                        os.path.basename(self.p_path)[:-3],
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
                """
    Plot a map of valid pixels and highlight a selected pixel location.

    This method visualizes a summary layer from the evaluation NetCDF dataset
    and overlays a user-specified pixel location. Cells with low or invalid
    values are masked, and the remaining valid cells are displayed in a
    grayscale map.

    Parameters
    ----------
    latitude : int or float
        Latitude index (row index) of the pixel to highlight.
    longitude : int or float
        Longitude index (column index) of the pixel to highlight.
    ax : matplotlib.axes.Axes
        Matplotlib axes object on which the map will be plotted.

    Returns
    -------
    matplotlib.image.AxesImage
        The image object returned by `imshow`, useful for adding a colorbar.

    Notes
    -----
    - The method reads the `"summary"` variable from the NetCDF dataset.
    - Values less than or equal to 2 are treated as invalid and masked.
    - The map is displayed using a grayscale colormap.
    - The highlighted pixel is plotted as a red star marker.
    - Axis labels refer to index space (not geographic coordinates).
    - The plot title includes a lake identifier derived from `self.p_path`.
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
                """
    Display an interactive pixel map and allow users to mark valid cells
    by clicking on them.

    This method plots the NetCDF `summary` layer on the given axes, masking
    cells with values less than or equal to 2. A mouse click callback is
    then attached to the figure so that clicking on a valid pixel adds a
    marker and its `(lat_idx, lon_idx)` indices to the plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes object on which the interactive map is drawn.

    Returns
    -------
    int
        The Matplotlib callback connection id returned by
        `mpl_connect("button_press_event", ...)`. This can be used later
        to disconnect the callback if needed.

    Notes
    -----
    - This function requires an interactive Matplotlib backend. In Jupyter
      environments, you must enable it using: %matplotlib widget
    - To change it back to the original backend, use: %matplotlib inline
    - The method reads the `"summary"` variable from `self.e_path`.
    - Cells with `summary <= 2` are treated as invalid and masked from the
      displayed image.
    - Clicked coordinates are rounded to the nearest integer grid indices.
    - Clicks outside the axes area are ignored.
    - Clicked indices are clipped to remain within array bounds.
    - If the clicked cell is invalid, no marker or label is added.
    - Valid clicks add:
        * a red circular marker at the selected pixel
        * a text label showing `"lat_idx,lon_idx"`
    - The plot title includes a lake identifier derived from `self.p_path`.
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
                lake_id = int(os.path.basename(self.p_path)[:-3])

                lake_row = self.geom[self.geom["id"] == lake_id]
                if lake_row.empty:
                        raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
                geom = lake_row.geometry.iloc[0]
                buffered_geom = self.shrink_geometry_by_1km(geom)
                buffered_geom_prepared = prep(buffered_geom)

                with netCDF4.Dataset(self.e_path) as nc:
                        summary = np.array(nc.variables["summary"][:, :])
                        lats = nc.variables["lat"][:]
                        lons = nc.variables["lon"][:]
                        map_data = np.full(summary.shape, np.nan)
                        for (i, j), r2 in metric_scores.items():
                                lon = lons[j]
                                lat = lats[i]
                                if buffered_geom_prepared.contains(Point(lon, lat)):
                                        map_data[i, j] = r2

                        if colormap:
                                cmap = colormap
                        else:
                                cmap = "RdYlBu"

                        # plot as geographic map
                        if colorbar_extent == [0,1]:
                                im = ax.imshow(map_data, cmap=cmap, aspect="auto", origin="lower", vmin= colorbar_extent[0], vmax = colorbar_extent[1], extent=[lons.min(), lons.max(), lats.min(), lats.max()])
                        else:
                                im = ax.imshow(map_data, cmap=cmap, aspect="auto", origin="lower", extent=[lons.min(), lons.max(), lats.min(), lats.max()])



                        label = False
                        if geom.geom_type == "Polygon":
                                x, y = geom.exterior.xy
                                label = "Lake Outline" if not label else None
                                ax.plot(x, y, color="black", linewidth=1, label = label)
                                label = True
                        elif geom.geom_type == "MultiPolygon":
                                for poly in geom.geoms:
                                        x, y = poly.exterior.xy
                                        label = "Lake Outline" if not label else None
                                        ax.plot(x, y, color="black", linewidth=1, label = label)
                                        label = True

                        ax.set_xlabel("Lon index")
                        ax.set_ylabel("Lat index")
                        text_str = f"{metric_str}-Scores for Lake: ID {os.path.basename(self.p_path)[:-3]}"
                        ax.set_title(text_str)
                        fig.colorbar(im, ax= ax, label= metric_str)
                        ax.legend()
                return im


        def interactive_metric_map(self, metric_scores, metric_str:str, fig, ax):
                lake_id = int(os.path.basename(self.p_path)[:-3])

                lake_row = self.geom[self.geom["id"] == lake_id]
                if lake_row.empty:
                        raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
                geom = lake_row.geometry.iloc[0]
                buffered_geom = self.shrink_geometry_by_1km(geom)
                buffered_geom_prepared = prep(buffered_geom)

                with netCDF4.Dataset(self.e_path) as nc:
                        summary = np.array(nc.variables["summary"][:, :])
                        lats = nc.variables["lat"][:]
                        lons = nc.variables["lon"][:]
                        map_data = np.full(summary.shape, np.nan)
                        for (i, j), r2 in metric_scores.items():
                                lon = lons[j]
                                lat = lats[i]
                                if buffered_geom_prepared.contains(Point(lon, lat)):
                                        map_data[i, j] = r2

                        # plot as geographic map
                        im = ax.imshow(map_data, cmap="RdYlBu", aspect="auto", origin="lower", vmin= 0, vmax = 1, extent=[lons.min(), lons.max(), lats.min(), lats.max()])

                        label = False
                        if geom.geom_type == "Polygon":
                                x, y = geom.exterior.xy
                                label = "Lake Outline" if not label else None
                                ax.plot(x, y, color="black", linewidth=1, label = label)
                                label = True
                        elif geom.geom_type == "MultiPolygon":
                                for poly in geom.geoms:
                                        x, y = poly.exterior.xy
                                        label = "Lake Outline" if not label else None
                                        ax.plot(x, y, color="black", linewidth=1, label = label)
                                        label = True

                        ax.set_xlabel("Lon index")
                        ax.set_ylabel("Lat index")
                        text_str = f"{metric_str} for Lake: ID {os.path.basename(self.p_path)[:-3]}"
                        ax.set_title(text_str)
                        fig.colorbar(im, ax= ax, label= metric_str)
                        ax.legend()


                def on_click(event):
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




        def extrema_plot(self, latitude, longitude, ax,  peak = True, aggregation= False,  start = 0, end = 9999, background_pts = True, purple_chla21= False):
                """
    Plot the time series and detected peak values for a selected pixel.

    This method visualizes observed data for a single pixel together with
    detected peak magnitudes and dates from the phenology parameter file.
    Optionally, it can plot spatially aggregated background values instead
    of the raw pixel time series. The displayed peaks and x-axis range can
    be restricted to a selected year interval.

    Parameters
    ----------
    latitude : int
        Row index of the selected pixel.
    longitude : int
        Column index of the selected pixel.
    ax : matplotlib.axes.Axes
        Matplotlib axes object on which the plot is drawn.
    aggregation : bool, optional
        If True, plot spatially aggregated background values from
        `self.aggregation_df` instead of the raw pixel time series.
        If False, plot the raw valid observations for the selected pixel.
        Default is False.
    start : int, optional
        Start year for the displayed period. If set to 0, the earliest
        available year in the valid time series is used. Default is 0.
    end : int, optional
        End year for the displayed period. If set to 9999, the latest
        available year in the valid time series is used. Default is 9999.

    Returns
    -------
    None
        This method plots directly onto `ax` and does not return a value.

    Notes
    -----
    - The method reads latitude, longitude, and time information from
      `self.e_path`, and smoothing/peak information from `self.p_path`.
    - Peak dates are read from `pks_x` and converted from Unix time to
      datetime objects.
    - Peak magnitudes are read from `pks_y`.
    - Only peaks whose years fall within the selected `[start, end]`
      interval are displayed.
    - Raw observations are filtered so that only values not equal to
      `-9999` and with QA flag equal to `0` are plotted.
    - If `aggregation=True` and `self.aggregation_df` has not yet been
      created, `self.spatial_aggregation()` is called automatically.
    - The plotted data color depends on the phenology variable inferred
      from the parent directory name of `self.p_path`:
        * `"phycocyanin"` -> red
        * `"chla_mean"` -> purple
        * otherwise -> green
    - Peaks are drawn with `ax.stem`, while observations/background values
      are drawn with `ax.scatter`.
    - Negative peak values are highlighted with red “x” markers and trigger
      a warning.
    - The title includes the lake ID and the geographic coordinates of the
      selected pixel.
    - The y-axis upper limit is adjusted based on the largest peak values.

    Warns
    -----
    UserWarning
        If spatial aggregation must be calculated before plotting.
    Warning
        If negative peaks are present in the selected time period.
    UserWarning
        If there is not enough valid data to produce the plot.
    """


                with netCDF4.Dataset(self.e_path) as nc:
                        # change np.array to np.asarry to avoid getting DeprecationError
                        lat = np.asarray(nc.variables["lat"])
                        lon = np.asarray(nc.variables["lon"])
                        t_all = unix_to_datenum(nc.variables["time"])


                with netCDF4.Dataset(self.p_path) as nc:
                        smoothing = float(nc.variables["smoothing_parameter"][latitude, longitude])
                        if peak:

                                pks_x = unix_to_datetime(remove_nan(nc.variables["pks_x"][latitude, longitude, :]))
                                pks_y = remove_nan(nc.variables["pks_y"][latitude, longitude, :])
                                mask_pks = np.array([(d.year <= end) & (d.year >= start) for d in pks_x])
                                x_sub = pks_x[mask_pks]
                                y_sub = pks_y[mask_pks]
                        else:
                                trgs_x = unix_to_datetime(remove_nan(nc.variables["trgs_x"][latitude, longitude, :]))
                                trgs_y = remove_nan(nc.variables["trgs_y"][latitude, longitude, :])
                                mask_trgs = np.array([(d.year <= end) & (d.year >= start) for d in trgs_x])
                                x_sub = trgs_x[mask_trgs]
                                y_sub = trgs_y[mask_trgs]

                        with netCDF4.Dataset(self.e_path) as nc:
                                variable = getattr(nc, "variable")
                                values = np.array(nc.variables[variable][:, latitude, longitude])
                                mask = (values != -9999) & (np.array(nc.variables[getattr(nc, 'qa')][:, latitude, longitude]) == 0)
                                values_m = values[mask]
                                time_m = t_all[mask]



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
                                if background_pts or purple_chla21:
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
                                                warnings.warn("Aggregation needs to be calculated, this may take a few minutes")
                                                self.spatial_aggregation()

                                        background_sub = self.aggregation_df[(self.aggregation_df["i"]==latitude) & (self.aggregation_df["j"]==longitude)]
                                        background_time = background_sub["time"]
                                        background_time = background_time.to_numpy()

                                        background_values = background_sub["MA_value"]

                                        ax.scatter(datenum_to_datetime(background_time), background_values, color=color, alpha=0.1, s=10, label=f"{label} Data")
                                elif background_pts and not aggregation:
                                        ax.scatter(datenum_to_datetime(time_m), values_m, color=color, alpha=0.1, s=10, label=f"{label} Data")
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


                                ax.legend(loc="upper left", ncol= 2)
                                if peak:
                                        textstr = f"Peak Comparison\n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                                else:
                                        textstr = f"Trough Comparison\n Lake ID:{os.path.basename(self.p_path)[:-3]}\n lat, lon: {round(float(lat[latitude]), 4)}, {round(float(lon[longitude]),4)}"
                                ax.set_title(textstr)
                                ax.set_ylabel("[ug/L]")
                                ax.grid()

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
                        



        def extrema_comparison(self, other1,  latitude, longitude, ax,  peak = True, aggregation= False, start = 0, end = 9999, background_pts = False, other2= None, purple_chla21= False):
                """
    Compare peak time series from two phenology datasets for the same lake
    and pixel location.

    This method overlays peak plots from `self` and another PhenologyEDA-like
    object on the same axes, allowing visual comparison of peak timing and
    magnitude for the same lake and grid cell. Optionally, spatially
    aggregated background values can be plotted instead of the raw pixel
    time series.

    Parameters
    ----------
    other : object
        Another object with the same interface as `self`, providing at least
        `p_path`, `peak_plot`, and compatible phenology data for comparison.
    latitude : int
        Row index of the selected pixel.
    longitude : int
        Column index of the selected pixel.
    ax : matplotlib.axes.Axes
        Matplotlib axes object on which both peak plots are drawn.
    aggregation : bool, optional
        If True, compare peaks against spatially aggregated background values
        instead of raw pixel time series. Default is False.
    start : int, optional
        Start year for the displayed comparison period. If set to 0, the
        earliest available year in each valid time series is used.
        Default is 0.
    end : int, optional
        End year for the displayed comparison period. If set to 9999, the
        latest available year in each valid time series is used.
        Default is 9999.

    Returns
    -------
    None
        This method plots directly onto `ax` and does not return a value.

    Raises
    ------
    Warning
        If `self` and `other` do not refer to the same lake ID.

    Notes
    -----
    - The lake ID is derived from the filename portion of `self.p_path` and
      `other.p_path`. Comparison is only allowed when both lake IDs match.
    - The actual plotting is delegated to `self.peak_plot(...)` and
      `other.peak_plot(...)`, both of which are drawn on the same axes.
    - The title is replaced after plotting to reflect the comparison between
      the two phenology variable names.
    - Phenology names are inferred from the parent directory names of
      `self.p_path` and `other.p_path`.
    - The title includes the lake ID and the geographic coordinates of the
      selected pixel, taken from `self.e_path`.
    - Any warnings or plotting behavior from `peak_plot` also apply here.
    """


                with netCDF4.Dataset(self.e_path) as nc:
                        # change np.array to np.asarry to avoid getting DeprecationError
                        lat = np.asarray(nc.variables["lat"])
                        lon = np.asarray(nc.variables["lon"])
                with netCDF4.Dataset(self.p_path) as nc:
                        lakeID1 = os.path.basename(self.p_path)[:-3]

                with netCDF4.Dataset(other1.p_path) as nc:
                        lakeID2 = os.path.basename(other1.p_path)[:-3]
                if other2:
                        with netCDF4.Dataset(other2.p_path) as nc:
                                lakeID3 = os.path.basename(other2.p_path)[:-3]

                if lakeID1 != lakeID2:
                        raise Warning("Comparison must be made on the same lake!")
                if other2:
                        if lakeID2!= lakeID3:
                                raise Warning("Comparison must be made on the same lake!")
                        ymax1 = self.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21)
                        ymax2 = other1.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21)
                        ymax3 = other2.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21)
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
                        ax.grid(True)

                else:

                        ymax1 = self.extrema_plot(latitude=latitude, longitude=longitude, ax = ax, peak= peak, aggregation = aggregation, start = start, end = end)
                        ymax2 = other1.extrema_plot(latitude=latitude, longitude=longitude, ax = ax,  peak= peak, aggregation = aggregation, start = start, end = end)
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
                        ax.set_ylabel("[ug/L]")
                        ax.grid()




        def single_plot_background(self, latitude, longitude, ax, fig, aggregation = False, start= 0, end= 9999):
                """
    Generate a time series plot for a single spatial point including raw data,
    smoothed spline, and detected phenological features (peaks, troughs, and
    midpoints), with optional spatial aggregation.

    This function reads environmental and processed data from NetCDF files,
    applies smoothing using cubic splines, computes performance metrics
    (RMSE and R²), and visualizes the results on a given matplotlib axis.

    Parameters
    ----------
    latitude : int
        Index of the latitude dimension in the dataset.
    longitude : int
        Index of the longitude dimension in the dataset.
    ax : matplotlib.axes.Axes
        Matplotlib axes object where the plot will be drawn.
    fig : matplotlib.figure.Figure
        Matplotlib figure object (used for colorbar placement).
    aggregation : bool, optional (default=False)
        If True, use spatially aggregated data instead of raw observations.
        Requires `self.aggregation_df` to be precomputed.
    start : int, optional (default=0)
        Start year (inclusive) for filtering the time series.
        If 0, the earliest available year is used.
    end : int, optional (default=9999)
        End year (inclusive) for filtering the time series.
        If 9999, the latest available year is used.

    Behavior
    --------
    - Loads latitude, longitude, and time information from the environmental dataset.
    - Extracts variable values and quality assurance (QA) flags for the selected grid point.
    - Filters out invalid values (-9999) and applies temporal subsetting.
    - Applies cubic smoothing spline (csaps) to generate a continuous signal.
    - Computes RMSE and R² for both the full dataset and the selected time subset.
    - Plots:
        * Raw or aggregated data (colored by QA flags)
        * Smoothed spline curve
        * Peaks, troughs, and midpoints (up/down)
    - Highlights negative feature values and issues warnings if present.
    - Adds legend, title, axis labels, and a QA-based colorbar.

    Notes
    -----
    - Requires external helper functions: `unix_to_datenum`, `unix_to_datetime`,
      `datenum_to_datetime`, and `remove_nan`.
    - Uses `csaps` for smoothing and `sklearn.metrics` for RMSE and R².
    - If aggregation is enabled and not yet computed, `self.spatial_aggregation()`
      will be called automatically.
    - Warnings are issued when:
        * No data is available
        * Not enough data points exist for metric computation
        * Negative peaks/troughs/midpoints are detected

    Returns
    -------
    None
        The function modifies the provided matplotlib axis in-place.
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
                                                warnings.warn("Aggregation needs to be calculated, this may take a few minutes")
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
                with netCDF4.Dataset(self.p_path) as nc:
                        pks_x = unix_to_datetime(remove_nan(nc.variables["pks_x"][latitude, longitude, :]))
                        mask_pks = np.array([(d.year <= end) & (d.year >= start) for d in pks_x])
                        pks_x_sub = pks_x[mask_pks]

                        return len(pks_x_sub)
                
        def pixel_r2(self, latitude, longitude, start= 0, end= 9999):
                scores = self.r2_scores(start=start, end = end)
                return scores[(latitude, longitude)]
        
        def pixel_rmse(self, latitude, longitude, start=0, end=9999):
                scores = self.RMSE_scores(start=start, end=end)
                return scores[(latitude, longitude)]

        def pixel_mad(self, latitude, longitude, start=0, end=9999):
                scores = self.MAD_scores(start=start, end=end)
                return scores[(latitude, longitude)]

        def pixel_correlation(self, latitude, longitude, start=0, end=9999):
                scores = self.correlation_scores(start=start, end=end)
                return scores[(latitude, longitude)]

        def pixel_values(self, latitude, longitude, start=0, end=9999):
                scores = self.values_per_pixel(start=start, end=end)
                return scores[(latitude, longitude)]



        def single_plot(self, latitude, longitude, ax, aggregation = False, start= 0, end= 9999):
                """
    Plot a smoothed phenology time series and detected phenological events for a
    single lake pixel.

    This method reads extract and phenology data from NetCDF files, selects the
    time series at the instance's latitude and longitude indices, and plots the
    observed values together with a cubic smoothing spline on the provided
    matplotlib axis. It also overlays detected phenological markers, including
    peaks, troughs, mid green-up, and mid green-down events, optionally filtered
    to a given year range.

    The function calculates RMSE and R² for both the full time series and the
    selected year subset, and displays these metrics on the plot. Negative event
    values are highlighted and reported with warnings.

    Args:
        ax (matplotlib.axes.Axes): Axis on which to draw the plot.
        start (int, optional): First year to include in the displayed subset.
            Defaults to 0, which uses the earliest year in the data.
        end (int, optional): Last year to include in the displayed subset.
            Defaults to 9999, which uses the latest year in the data.

    Returns:
        None: The function modifies `ax` in place and does not return a value.

    Notes:
        - Observations with value `-9999` or nonzero QA flags are excluded.
        - Spline smoothing uses the `smoothing_parameter` stored in the
          phenology dataset.
        - Negative peak, trough, mid-up, or mid-down values are marked with red
          "x" symbols and raise warnings.
        - Plot limits and annotations are adjusted automatically based on the
          selected data.

    Raises:
        Warning: If negative phenology values are detected in the selected time
            period.
        Warning: If no valid data are available to plot.
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
                        values = np.array(nc.variables[variable][:, latitude, longitude])
                        mask = (values != -9999) & (np.array(nc.variables[getattr(nc, 'qa')][:, latitude, longitude]) == 0)
                        values_m = values[mask]
                        time_m = t_all[mask]

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
                                                warnings.warn("Aggregation needs to be calculated, this may take a few minutes")
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
                                ax.grid()
                                ax.set_ylabel("[ug/L]")
                                ax.set_xlim(pd.to_datetime('01-01-' + str(function_start), format='%d-%m-%Y') , pd.to_datetime('31-12-' + str(function_end), format='%d-%m-%Y'))


                                pks_lim_sub = sorted(pks_y_sub)
                                if max(pks_lim_sub)> 10:
                                        ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-2]+0.5)
                                else:
                                        ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-1]+0.5)

                                if neg_values_sub:
                                        ax.text(0.99,0.99,f"# Neg.values: {sum(neg_values_sub)} \n RMSE: {round(rmse_sub,4)}\n R$^2$: {round(r2_sub,4)}\n MAD: {round(mad_sub, 4)}", transform = ax.transAxes,   ha= "right", va= "top")
                                else:
                                        ax.text(0.99,0.99,f"RMSE:{round(rmse_sub,4)} \n R$^2$: {round(r2_sub,4)}\n MAD: {round(mad_sub, 4)}", transform = ax.transAxes,   ha= "right", va= "top")
                        else:
                                warnings.warn("Not enough data to plot or compute metrics for chosen time interval")
                else:
                        warnings.warn("No data to plot")


        def split_plot(self, latitude, longitude, ax0, ax1, aggregation = False, start0= 0, end0= 9999, start1= 0, end1=9999):
                """
    Plot two separate time ranges of the phenology data on two axes.

    This method creates a split visualization by calling `single_plot` twice,
    allowing comparison of two different time periods side by side (or stacked),
    each rendered on its own matplotlib axis.

    Args:
        ax0 (matplotlib.axes.Axes): Axis for the first time range.
        ax1 (matplotlib.axes.Axes): Axis for the second time range.
        start0 (int, optional): Start year for the first plot. Defaults to 0
            (uses earliest available year).
        end0 (int, optional): End year for the first plot. Defaults to 9999
            (uses latest available year).
        start1 (int, optional): Start year for the second plot. Defaults to 0.
        end1 (int, optional): End year for the second plot. Defaults to 9999.

    Returns:
        None: The function modifies the provided axes in place.

    Notes:
        - Internally calls `self.single_plot` for each axis.
        - A warning is issued if both plots use the full time range, as this
          defeats the purpose of splitting the visualization.
    """
                if start0 and start1 == 0 and end0 and end1 == 9999:
                        warnings.warn("split_plot needs a least end0 and start1 parameter, otherwise use full_plot")

                self.single_plot(latitude = latitude, longitude= longitude, ax=ax0, aggregation = aggregation, start= start0, end = end0)
                self.single_plot(latitude = latitude, longitude= longitude, ax=ax1, aggregation = aggregation, start=start1, end=end1)

        def full_plot(self, latitude, longitude, ax, aggregation = False):
                """
    Plot the full available phenology time series on a single axis.

    This method extracts the complete valid time series for the selected
    latitude and longitude from the extract dataset (`self.e_path`), determines
    the full temporal extent of the data, and delegates plotting to
    `self.single_plot`.

    Invalid values (marked as -9999) and observations with nonzero QA flags
    are removed before determining the time range.

    Args:
        ax (matplotlib.axes.Axes): Axis on which to draw the plot.

    Returns:
        None: The function modifies `ax` in place and does not return a value.

    Notes:
        - The full time range is automatically determined from the filtered data.
        - Internally calls `self.single_plot` with the computed start and end years.
        - Requires at least two valid data points to produce a plot.
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

















