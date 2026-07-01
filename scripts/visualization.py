import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import numpy as np
import datetime
import os
from pathlib import Path
import csv
import statistics
import warnings
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, Patch
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr, gaussian_kde
from csaps import csaps
import functions as f
import multiprocessing
from functools import partial
import shapely.ops as ops
from pyproj import CRS, Transformer
import geopandas
from shapely.prepared import prep
from shapely.geometry import Point
from numpy.lib.stride_tricks import sliding_window_view
import colorcet as cc
import time



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
            _GLOBALS[f"kde_{var}"] = np.asarray(nc.variables[var][:])


# color_sets_4x4: bivariate color palettes for the 4×4 heatmap legend.
# Color palettes are adapted from:
# https://medium.com/@leodpereda/mastering-bivariate-maps-with-plotly-a-step-by-step-guide-ad9cae150d8a
color_sets_4x4 = {
    'pink-blue': [
        '#e8e8e8', '#cfdfe6', '#95d3d8', '#5ac8c8',
        '#e0c0d5', '#c4bfd7', '#8fb9d1', '#5698b9',
        '#d398c1', '#b29fcb', '#7f93be', '#4a72a7',
        '#be64ac', '#9762ae', '#6760a5', '#3b4994'
    ],
    'teal-red': [
    '#e8e8e8', '#e7c0c0', '#dc8b8b', '#c85a5a',
    '#c1dde4', '#d2b7d5', '#c08cb9', '#985356',
    '#8fd0da', '#b69dcf', '#9d77b4', '#7f5a9e',
    '#64acbe', '#8c7ec0', '#6c5cad', '#4b3f98'
    ],
    'teal-red1': [
        '#e7e7e7', '#cfd9dd', '#9abec9', '#64abbd',
        '#ddc4c4', '#c1b4b8', '#949ca3', '#5d7783',
        '#d08e8e', '#b18488', '#8c62aa', '#7058a8',
        '#c75a5a', '#a15355', '#5a4ba5', '#5a4ba5'
    ],
    'blue-orange': [
        '#fef1e4', '#fcd1b7', '#f8a474', '#f3742d',
        '#c8e1e6', '#c6c0bc', '#bb917c', '#ab5f37',
        '#7fc8e6', '#7ea2bb', '#7a7879', '#6f5a47',
        '#18aee5', '#2b96c7', '#3a7ea9', '#5c473d'
    ]
}


class PhenologyVisualization:
    shapefile_path = None
    QA_LEVELS = (0, 1, 2)
    
    QA_CONFIG = {
        0: {
            "label": "Good",
            "style":{"marker": "o", "color": "green"},
            "style_alt":{"marker": "o"}
        },
        1: {
            "label": "Fair",
            "style":{"marker": "o", "color": "orange"},
             "style_alt":{"marker": "s"}
        },
        2: {
            "label": "Poor",
            "style":{"marker": "o", "color": "red"},
             "style_alt":{"marker": "x"}
        },
    }


    VAR_CONFIG = {
        "phycocyanin": {
            "label": "phyco",
            "style": {"color": "blue"},
        },
        "chla_mean": {
            "label": "chla v2.1",
            "style": {"color": "lightgreen"},
            "style_alt": {"color_alt": "purple"},
        },
        "chla": {
            "label": "chla v3.0",
            "style": {"color": "green"},
        },
    }

    METRIC_CONFIG = {
        'pks': {
            'label':'peak',
            'style':{'marker':'o','color':'blue','label':'peak'}
        },
        'trg': {
            'label':'trough',
            'style':{'marker':'o','color':'orange','label':'trough'}
        },
        'midUP': {
            'label':'green up',
            'style':{'marker':'^','color':'mediumseagreen'},
            'meta':{}
        },
        'midDOWN': {
            'label':'green down',
            'style':{'marker':'v','color':'darkgreen'},
            'meta':{}
        },
        'onsetUP': {
            'label':'green up',
            'style':{'marker':'^','color':'mediumseagreen'},
            'meta':{}
        },
        'onsetDOWN': {
            'label':'green down',
            'style':{'marker':'v','color':'darkgreen'},
            'meta':{}
        },
        'advUP': {
            'label':'green up',
            'style':{'marker':'^','color':'mediumseagreen'},
            'meta':{}
        },
        'advDOWN': {
            'label':'green down',
            'style':{'marker':'v','color':'darkgreen'},
            'meta':{}
        },
    }


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
            self.gdf = geopandas.read_file(self.shapefile_path)
        self.p_path = phenology_path
        self.e_path = extract_path
        self.version = Path(self.p_path).parents[2].name.removeprefix('v')
        self.variable = Path(self.p_path).parents[0].stem
        self.lakeID = Path(self.p_path).stem
        # self.lakename = lakeID_to_name(self.gdf,self.lakeID)
        self.info = (f"Version: {self.version} \n",
                    f"Variable: {self.variable} \n" ,
                    f"Lake ID: {self.lakeID}")
        self.methods = [method for method in dir(PhenologyVisualization) if callable(getattr(PhenologyVisualization, method)) and not method.startswith("__")]
        self.valid_coords = self.valid_index_pairs()
        self.out_folder = Path(self.p_path).parents[2]
        self.aggregation_df = None
        self.geom_shrunk = None
        self._extracted_globals = None
        self._pixel_cache = {}
        self._lake_cache = {}
        self.prep_geometry_from_shapefile()


    def ID_to_name(self, id):
        """Return the lake name for a given lake ID from the shapefile.

        Parameters
        ----------
        id : int
            Lake ID matching the 'id' column in the CCI shapefile.

        Returns
        -------
        str
            Lake name from the 'name' column.

        Raises
        ------
        ValueError
            If the ID is not found in the shapefile.
        """
        matches = self.gdf.loc[self.gdf["id"] == int(id), "name"]
        if matches.empty:
            raise ValueError(f"Lake ID {id} not found in shapefile.")
        return matches.iloc[0]


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
            

    def create_DataFrame(self, latitude_idx, longitude_idx):
        """Build a combined long-format DataFrame of all phenology variables for a single pixel.

        Reads all _x (time) and _y (value) variable pairs from the phenology NetCDF
        and concatenates them into a single DataFrame.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel in the grid.
        longitude_idx : int
            Column (lon) index of the pixel in the grid.

        Returns
        -------
        pandas.DataFrame
            Long-format DataFrame with columns Value, Variable, latitude_idx, longitude_idx
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
            var_x =f.unix_to_datetime(f.remove_nan(p[x][latitude_idx,longitude_idx,:]))
            var_y = f.remove_nan(p[y][latitude_idx,longitude_idx,:])
            var_label = [x[:-2]]*len(var_y)
            df = pd.DataFrame({"Value":var_y,
                            "Variable": var_label,
                            "latitude": lat[latitude_idx],
                            "longitude": lon[longitude_idx], 
                            "lake_ID": self.lakeID},
                            index = var_x)
            df.index.names = ["Time"]
            result[y[:-2]] = df
        combined_df = pd.concat(result.values())

        return combined_df


    def get_plot_config(self, config_name, key, **kwargs):
        """
        Generic accessor for QA, variable, and metric configs.
        """

        config_map = {
            "qa": self.QA_CONFIG,
            "var": self.VAR_CONFIG,
            "metric": self.METRIC_CONFIG,
        }

        cfg = config_map[config_name][key]

        # Base style
        style = cfg.get("style", {}).copy()

        # Optional override (e.g. purple_chla21)
        if kwargs.get("use_alt") and "style_alt" in cfg:
            alt = {("color" if k == "color_alt" else k): v for k, v in cfg["style_alt"].items()}
            style.update(alt)

        return {
            "label": cfg.get("label"),
            "style": style,
            "meta": cfg.get("meta", {}),
        }
    

    def shrink_geometry(self, geom=None,distance_m= 1000):
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
        if geom is None:
            geom = self.geometry
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
        geom_shrunk_local = geom_local.buffer(-distance_m)   # minus 1000 m = inward 1 km
        geom_shrunk = ops.transform(to_wgs84, geom_shrunk_local)

        if geom_shrunk_local.is_empty:
            raise ValueError("Geometry disappeared after shrinking (buffer too large).")

        # self.geom_shrunk = geom_shrunk
        return geom_shrunk


    def extract_geometry_from_shapefile(self):
        lake_id = int(self.lakeID)
        lake_row = self.gdf[self.gdf["id"] == lake_id]
        if lake_row.empty:
            raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
        
        # self.geometry = lake_row.geometry.iloc[0]
        return lake_row.geometry.iloc[0]
    

    def prep_geometry_from_shapefile(self):
        geom = self.extract_geometry_from_shapefile()
        geom_shrunk = self.shrink_geometry(geom)

        self.geometry = geom
        self.geom_shrunk = geom_shrunk
        self.prepped_geom = prep(geom_shrunk)
    

    @staticmethod
    def compute_metric_score(coord, start=0, end=9999, metrics_to_compute= None):
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


    def build_metric_path(self, metric_name, start=0, end=9999):
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
        base = os.path.join(self.out_folder,self.lakeID, "calculated_values", "metrics", metric_name, f"v{self.version}", self.variable)
        if start == 0 and end == 9999:
            fname = "full_ts.csv"
        elif start == 0:
            fname = f"ts_{2002}_to_{end}.csv"
        elif end == 9999:
            fname = f"ts_{start}_to_{2024}.csv"
        else:
            fname = f"ts_{start}_to_{end}.csv"
        
        return base, os.path.join(base, fname)

    def build_kde_path(self):
        """Return the output directory and file path for the cached KDE events CSV.

        Returns
        -------
        base : str
            Directory path where the CSV will be written.
        file_path : str
            Full path to the CSV file.
        """
        lake_name = self.ID_to_name(int(self.lakeID)).replace(" ", "")
        base = os.path.join(
            self.out_folder.parents[1], "lake_analysis",
             f"ID{self.lakeID}_{lake_name}", "calculated_values", "kde_data",
            f"v{self.version}", self.variable,
        )
        base, path =  base, os.path.join(base, "kde_events.csv")
        print(path)
        return base, path

    def compute_and_cache_metric(self, metric_name, col_name, compute_fn, start=0, end=9999):
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
        workers = partial(compute_fn, start=start, end=end, metrics_to_compute = [metric_name])
        with multiprocessing.Pool(initializer=_init_worker, initargs=(self.p_path, self.e_path), processes=3) as pool:
            result = pool.map(workers, self.valid_coords)
        data = dict(result)
        t_df = pd.DataFrame([(i, j, v) for (i, j), v in data.items()], columns=["i", "j", col_name])
        t_df.to_csv(file_path, index=False)
        return data


    def r2_scores(self, time_split=None):
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
        if time_split is None:
            time_split = [(0, 9999)]
        for start, end in time_split:
            return self.compute_and_cache_metric(metric_name="r2", col_name="r2_scores", compute_fn=PhenologyVisualization.compute_metric_score, start=start, end=end)


    def MAD_scores(self, time_split=None):
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
        if time_split is None:
            time_split = [(0, 9999)]
        for start, end in time_split:
            return self.compute_and_cache_metric(metric_name="MAD", col_name="mad_scores",compute_fn= PhenologyVisualization.compute_metric_score,start= start,end= end)


    def RMSE_scores(self, time_split=None):
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
        if time_split is None:
            time_split = [(0, 9999)]
        for start, end in time_split:
            return self.compute_and_cache_metric(metric_name="RMSE", col_name="rmse_scores", compute_fn=PhenologyVisualization.compute_metric_score, start=start, end=end)


    def correlation_scores(self, time_split=None):
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
        if time_split is None:
            time_split = [(0, 9999)]
        for start, end in time_split:
            return self.compute_and_cache_metric(metric_name="correlation", col_name="correlation_scores", compute_fn=PhenologyVisualization.compute_metric_score, start=start, end=end)


    def values_per_pixel(self, time_split=None):
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
        if time_split is None:
            time_split = [(0, 9999)]
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


    def pixel_map(self, latitude_idx, longitude_idx, ax):
        """Plot a grayscale coverage map with the selected pixel marked.

        Reads the summary grid from the extract NetCDF and displays it in
        grayscale, masking cells with fewer than 2 valid observations. The
        selected pixel is overlaid as a red star.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel to highlight.
        longitude_idx : int
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
        ax.plot(longitude_idx, latitude_idx, "r*", markersize=14, zorder=5, label="Pixel")

        ax.set_xlabel("Lon index")
        ax.set_ylabel("Lat index")
        text_str = f"Pixel Location\n Lake ID:{self.lakeID}"
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
        ax.set_title(f"Pixel Location\n Lake ID: {self.lakeID}")

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


    def metric_map(self, metric_scores:dict, metric_str:str, fig, ax, colormap = None, colorbar_extent= None):
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
        if colorbar_extent is None:
            colorbar_extent = [0,1]
        lake_id = int(self.lakeID)
        lake_row = self.gdf[self.gdf["id"] == lake_id]
        if lake_row.empty:
            raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
        
        geom = lake_row.geometry.iloc[0]
        buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = prep(buffered_geom)

        map_data, extent = f.grab_metrics(self.e_path, metric_scores, buffered_geom_prepared)

        im = f.plot_map_data(colormap, map_data, extent, ax, cmap_extent=colorbar_extent)
        f.plot_lake_outline(geometry=geom, ax=ax)
        f.set_labels(ax, fig, im,
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
        lake_id = int(self.lakeID)

        lake_row = self.gdf[self.gdf["id"] == lake_id]
        if lake_row.empty:
                raise ValueError(f"Lake ID {lake_id} not found in shapefile.")
        geom = lake_row.geometry.iloc[0]
        buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = prep(buffered_geom)

        map_data, extent = f.grab_metrics(self.e_path, metric_scores, buffered_geom_prepared)
        g = self._load_extracted_globals()
        lats = g["lat"]
        lons = g["lon"]

        im = f.plot_map_data(None, map_data, extent, ax, cmap_extent=[0, 1])
        f.plot_lake_outline(geometry=geom, ax=ax)
        f.set_labels(ax, fig, im,
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
    

    # def grab_DOY_data(self):
    #     return map_data, extent

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
        lake_id = int(self.lakeID)

        # lake_row = self.gdf[self.gdf["id"] == lake_id]
        geom = self.geometry #lake_row.geometry.iloc[0]
        # buffered_geom = self.shrink_geometry(geom)
        buffered_geom_prepared = self.prepped_geom

        var_x = "pks_x" if peaks else "green_up_mid_x"
        var_y = "pks_y" if peaks else "green_up_mid_y"
        extrema_label = "Peak" if peaks else "Green Mid Up"

        map_data, extent = f.grab_time_data(self.e_path, self.p_path, self.valid_coords,
                                            buffered_geom_prepared, var_x, var_y, year, max)
        im = f.plot_map_data("rainbow", map_data, extent, ax, cmap_extent=[160, 250])
        f.plot_lake_outline(geometry=geom, ax=ax)
        if colorbar:
            f.set_labels(ax, fig, im,
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
            Typically obtained from a pixel time series, e.g.:
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
        idx = np.argwhere(f.datenum_to_datetime(t_all) == date)
        if len(idx) == 0:
            raise ValueError(f"Date {date} not found in the extract time axis.")
        time_index = idx[0, 0]

        with netCDF4.Dataset(self.e_path) as nc:
            values = np.array(nc.variables[g["variable"]][time_index, :, :])
            qa     = np.array(nc.variables[g["qa"]][time_index, :, :])
            lats   = nc.variables["lat"][:]
            lons   = nc.variables["lon"][:]

        mask = (values != -9999) & (qa == 0)
        lake_id = int(self.lakeID)
        lake_row = self.gdf[self.gdf["id"] == lake_id]
        geom = lake_row.geometry.iloc[0]
        buffered_geom = self.shrink_geometry(geom)
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
                    "t_all":    f.unix_to_datenum(nc.variables["time"]),
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
                pks_x_raw = np.array(nc.variables["pks_x"][i, j, :])
                pk_mask   = ~np.isnan(pks_x_raw)
                pks_x  =f.unix_to_datetime(pks_x_raw[pk_mask])
                pks_y  = np.array(nc.variables["pks_y"][i, j, :])[pk_mask]
                pks_qa = np.array(nc.variables["pks_qa"][i, j, :])[pk_mask]
                trgs_x_raw = np.array(nc.variables["trgs_x"][i, j, :])
                trg_mask   = ~np.isnan(trgs_x_raw)
                trgs_x  = f.unix_to_datetime(trgs_x_raw[trg_mask])
                trgs_y  = np.array(nc.variables["trgs_y"][i, j, :])[trg_mask]
                trgs_qa = np.array(nc.variables["trgs_qa"][i, j, :])[trg_mask]
                midUP_x    = f.unix_to_datetime(f.remove_nan(nc.variables["green_up_mid_x"][i, j, :]))
                midUP_y    = f.remove_nan(nc.variables["green_up_mid_y"][i, j, :])
                midDOWN_x  = f.unix_to_datetime(f.remove_nan(nc.variables["green_down_mid_x"][i, j, :]))
                midDOWN_y  = f.remove_nan(nc.variables["green_down_mid_y"][i, j, :])
                onsetUP_x    = f.unix_to_datetime(f.remove_nan(nc.variables["green_up_onset_x"][i, j, :]))
                onsetUP_y    = f.remove_nan(nc.variables["green_up_onset_y"][i, j, :])
                onsetDOWN_x  = f.unix_to_datetime(f.remove_nan(nc.variables["green_down_onset_x"][i, j, :]))
                onsetDOWN_y  = f.remove_nan(nc.variables["green_down_onset_y"][i, j, :])
                advUP_x    = f.unix_to_datetime(f.remove_nan(nc.variables["green_up_advanced_x"][i, j, :]))
                advUP_y    = f.remove_nan(nc.variables["green_up_advanced_y"][i, j, :])
                advDOWN_x  = f.unix_to_datetime(f.remove_nan(nc.variables["green_down_advanced_x"][i, j, :]))
                advDOWN_y  = f.remove_nan(nc.variables["green_down_advanced_y"][i, j, :])

                gap_starts = f.unix_to_datetime(f.remove_nan(nc.variables["data_gap_start"][i, j, :]))
                gap_ends   = f.unix_to_datetime(f.remove_nan(nc.variables["data_gap_end"][i, j, :]))
            self._pixel_cache[(i,j)] = {
                "values": values, "qa": qa, "smoothing": smoothing,
                "pks_x": pks_x, "pks_y": pks_y, "pks_qa": pks_qa,
                "trgs_x": trgs_x, "trgs_y": trgs_y, "trgs_qa": trgs_qa,
                "midUP_x": midUP_x, "midUP_y": midUP_y,
                "midDOWN_x": midDOWN_x, "midDOWN_y": midDOWN_y,
                "onsetUP_x": onsetUP_x, "onsetUP_y": onsetUP_y,
                "onsetDOWN_x": onsetDOWN_x, "onsetDOWN_y": onsetDOWN_y,
                "advUP_x": advUP_x, "advUP_y": advUP_y,
                "advDOWN_x": advDOWN_x, "advDOWN_y": advDOWN_y,
                "gap_starts": gap_starts, "gap_ends": gap_ends,
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
        time_dt   = f.datenum_to_datetime(t_all[mask])
        return pd.Series(index=time_dt,data=values_m)

    


    def extrema_plot(self, latitude_idx, longitude_idx, ax,  peak = True, aggregation= False,  start = 0, end = 9999, background_pts = True, purple_chla21= False, show_legend = True):
        """Plot detected peaks or troughs as a stem plot with optional background scatter.

        Displays summer peaks or winter troughs for the pixel at (latitude_idx, longitude_idx)
        as vertical stems, with each extremum marker shaped by its QA flag
        while keeping the product's colour. Background observations may be shown as
        a raw scatter (QA==0 only) or 3×3 spatial median (when aggregation=True). Negative
        values are flagged with red crosses and trigger a warning.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        px = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        smoothing = px["smoothing"]
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])

        var = "pks" if peak else "trgs"
        plotting_data = f.grab_plotting_variables(start=start, end=end, pixel_data=px, variables=[var])
        x_sub, y_sub, qa_sub = plotting_data[var]

        mask     = (px["values"] != -9999) & (px["qa"] == 0)
        values_m = px["values"][mask]
        time_m   = t_all[mask]

        
        if len(values_m) <= 1:
            warnings.warn("No data to plot (check valid indices)")
            return None

        limits = sorted(f.datenum_to_datetime(time_m))
        function_start = min(limits).year if start == 0 else start
        function_end= max(limits).year if end ==9999 else end

        phenology_name = self.variable
        
        var_cfg = self.get_plot_config(
            "var",
            phenology_name,
            use_alt=purple_chla21
        )
        
        var_label = var_cfg["label"]
        var_style = var_cfg["style"]
        background_style = {"alpha": 0.1, "s": 10}
    
        if background_pts:
            self.plot_background_pts(
                ax=ax,
                latitude_idx = latitude_idx,
                longitude_idx = longitude_idx,
                masked_values = values_m,
                masked_time = time_m,
                aggregation= aggregation,
                **{**var_style, **background_style}
            )
        elif not background_pts and aggregation:
            warnings.warn(f"Aggregation ignored for lake ID {self.lakeID} since backrgound_pts turned off.")

        self.plot_data_gaps(ax=ax, pixel_data=px)

        ax.stem(x_sub, y_sub, markerfmt=" ", basefmt = " ",linefmt = var_style['color'])
        
        seen_labels = set()
        for q in self.QA_LEVELS:
            qm = qa_sub == q
            if not qm.any():
                continue

            qa_cfg = self.get_plot_config("qa", q, use_alt = True)
            qa_style = qa_cfg["style"]

            qa_label = qa_cfg["label"]
            if qa_label in seen_labels:
                qa_label = None
            else:
                seen_labels.add(qa_label)

            combined_style = {
                **var_style,            # color
                **qa_style,             # marker
                "s": 50,
                "edgecolors": var_style['color'],
                "linewidths": 2,
                "zorder": 4,
            }

            ax.scatter(
                x_sub[qm],
                y_sub[qm],
                label=qa_label,
                **combined_style,
            )

        if (y_sub < 0).any():
            mask =  y_sub<0
            ax.scatter(x_sub[mask], y_sub[mask], color="red", s=50, marker="x", zorder=6, label="Negative value")
            warnings.warn(f"Negative Peak(s) in time period {start}-{end}", Warning)

        if show_legend:
            ax.legend(loc="upper left", ncol= 2)
        
        plot_var = "Peak" if peak else "Trough"
        textstr = f"{plot_var} Comparison\n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}"

        ax.set_title(textstr)
        ax.xaxis.set_minor_locator(mdates.YearLocator())
        ax.grid(axis="x", which="minor", linewidth=0.5)
        ax.grid(axis="x", which="major", linewidth=0.5)
        ax.grid(axis="y", linewidth=0.5)
        ax.set_ylabel("[ug/L]")

        ax.set_xlim(
            pd.to_datetime('01-01-' + str(function_start), format='%d-%m-%Y'),
            pd.to_datetime('31-12-' + str(function_end), format='%d-%m-%Y')
        )
        pks_lim_sub = sorted(y_sub)
        if max(pks_lim_sub)> 10:
            ymax = pks_lim_sub[-2]+0.5
            ax.set_ylim(-0.5, ymax)
        else:
            ymax = pks_lim_sub[-1]+0.5
            ax.set_ylim(-0.5, ymax)
        return ymax
                    

    def extrema_comparison(self, other1,  latitude_idx, longitude_idx, ax,  peak = True, aggregation= False, start = 0, end = 9999, background_pts = True, other2= None, purple_chla21= False, show_legend= False):
        """Overlay extrema plots from two or three PhenologyVisualization instances on one axis.

        Calls extrema_plot for self and other1 (and optionally other2), sharing the
        same axes so that peaks or troughs from different products (e.g. chla v2.1
        vs v3.0) can be compared directly. All instances must reference the same lake.

        Parameters
        ----------
        other1 : PhenologyVisualization
            Second instance to overlay (must share the same lake ID as self).
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])

        lakeID1 = self.lakeID
        lakeID2 = other1.lakeID
        if other2:
            lakeID3 = other2.lakeID

        if lakeID1 != lakeID2:
            raise Warning("Comparison must be made on the same lake!")
        if other2:
            if lakeID2!= lakeID3:
                    raise Warning("Comparison must be made on the same lake!")
            ymax1 = self.extrema_plot(latitude_idx=latitude_idx, longitude_idx=longitude_idx, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
            ymax2 = other1.extrema_plot(latitude_idx=latitude_idx, longitude_idx=longitude_idx, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
            ymax3 = other2.extrema_plot(latitude_idx=latitude_idx, longitude_idx=longitude_idx, ax = ax, peak = peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
            y_lims = [ymax1, ymax2, ymax3]
            phenology_name1 = self.variable
            phenology_name2 = other1.variable
            phenology_name3 = other2.variable
            
            label_dict = {"phycocyanin": "phyco",
                            "chla_mean": "chla v2.1",
                            "chla": "chla v3.0"
                            }
            if purple_chla21:
                color_dict = {"phycocyanin": "blue",
                                "chla_mean": "purple",
                                "chla": "green"
                                }
            else:
                color_dict = {"phycocyanin": "blue",
                                "chla_mean": "lightgreen",
                                "chla": "green"
                                }

            if peak:
                textr =  f"{label_dict[phenology_name1]}, {label_dict[phenology_name2]} vs {label_dict[phenology_name3]} Peaks \n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}"
            else:
                textr =  f"{label_dict[phenology_name1]}, {label_dict[phenology_name2]} vs {label_dict[phenology_name3]} Troughs \n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}"
            ax.set_title(textr)
            ax.set_ylim(top = max(y_lims))
            ax.xaxis.set_minor_locator(mdates.YearLocator())
            ax.grid(axis="x", which="minor", linewidth=0.5)
            ax.grid(axis="x", which="major", linewidth=0.5)
            ax.grid(axis="y", linewidth=0.5)
            ax.set_ylabel("[ug/L]")
            # QA legend (markers only, no year lines)
            qa_markers = {0: "o", 1: "s", 2: "x"}
            qa_labels  = {0: "Good", 1: "Fair", 2: "Poor"}
            qa_handles = [mlines.Line2D([], [], color="black", marker=qa_markers[qa],
                            linestyle="None", markersize=6, label=qa_labels[qa])
                    for qa in self.QA_LEVELS
            ]
            # Add data gap legend entry
            qa_handles.append(Patch(facecolor="orange",
                edgecolor="orange",
                alpha=0.15,
                label="Data gap")
            )

            # Color legend for product versions
            type_handles = [
            mlines.Line2D([], [], color=color_dict[phenology_name1], marker="o",
                        linestyle="None", markersize=8, label=label_dict[phenology_name1]),
            mlines.Line2D([], [], color=color_dict[phenology_name2], marker="o",
                        linestyle="None", markersize=8, label=label_dict[phenology_name2]),
            mlines.Line2D([], [], color=color_dict[phenology_name3], marker="o",
                        linestyle="None", markersize=8, label=label_dict[phenology_name3]),]

            leg1 = ax.legend(handles=qa_handles, loc="upper left")
            ax.add_artist(leg1)
            ax.legend(handles= type_handles, loc = "upper right")
        else:
            ymax1 = self.extrema_plot(latitude_idx=latitude_idx, longitude_idx=longitude_idx, ax = ax, peak= peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
            ymax2 = other1.extrema_plot(latitude_idx=latitude_idx, longitude_idx=longitude_idx, ax = ax,  peak= peak, aggregation = aggregation, start = start, end = end, background_pts=background_pts, purple_chla21=purple_chla21, show_legend=show_legend)
            y_lims = [ymax1, ymax2]
            phenology_name1 = self.variable
            phenology_name2 = other1.variable

            label_dict = {"phycocyanin": "phyco",
                                    "chla_mean": "chla v2.1",
                                    "chla": "chla v3.0"
                                    }
            if purple_chla21:

                    color_dict = {"phycocyanin": "blue",
                                    "chla_mean": "purple",
                                    "chla": "green"
                                    }
            else:
                    color_dict = {"phycocyanin": "blue",
                                    "chla_mean": "lightgreen",
                                    "chla": "green"
                                    }
            if peak:
                    textr =  f"{label_dict[phenology_name1]} vs {label_dict[phenology_name2]} Peaks \n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}"
            else:
                    textr =  f"{label_dict[phenology_name1]} vs {label_dict[phenology_name2]} Troughs \n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}"
            
            ax.set_title(textr)
            ax.set_ylim(top = max(y_lims))
            ax.xaxis.set_minor_locator(mdates.YearLocator())
            ax.grid(axis="x", which="minor", linewidth=0.5)
            ax.grid(axis="x", which="major", linewidth=0.5)
            ax.grid(axis="y", linewidth=0.5)
            ax.set_ylabel("[ug/L]")
            # QA legend (markers only, no year lines)
            qa_markers = {0: "o", 1: "s", 2: "x"}
            qa_labels  = {0: "Good", 1: "Fair", 2: "Poor"}
            qa_handles = [mlines.Line2D([], [], color="black", marker=qa_markers[qa],
                            linestyle="None", markersize=6, label=qa_labels[qa])
                    for qa in self.QA_LEVELS
            ]
            # Add data gap legend entry
            qa_handles.append(Patch(facecolor="orange",
                    edgecolor="orange",
                    alpha=0.15,
                    label="Data gap")
            )

            # Color legend for product versions
            type_handles = [
            mlines.Line2D([], [], color=color_dict[phenology_name1], marker="o",
                        linestyle="None", markersize=8, label=label_dict[phenology_name1]),
            mlines.Line2D([], [], color=color_dict[phenology_name2], marker="o",
                        linestyle="None", markersize=8, label=label_dict[phenology_name2]),]

            leg1 = ax.legend(handles=qa_handles, loc="upper left")
            ax.add_artist(leg1)
            ax.legend(handles= type_handles, loc = "upper right")


    def single_plot_background(self, latitude_idx, longitude_idx, ax, fig, aggregation = False, start= 0, end= 9999):
        """Plot a pixel time series with QA-coloured scatter, spline, and phenological events.

        Like single_plot, but colours each background scatter point by its QA flag
        using a discrete colormap and adds a QA colorbar to the figure. Peaks,
        troughs, green-up and green-down midpoints are overlaid as scatter markers.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        g  = self._load_extracted_globals()
        pixel_data = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        smoothing = pixel_data["smoothing"]

        plotting_data = f.grab_plotting_variables(start=start, end=end, pixel_data=pixel_data)

        # No QA==0 filter here — all non-fill values kept for QA-coloured scatter
        mask     = (pixel_data["values"] != -9999)
        values_m = pixel_data["values"][mask]
        time_m   = t_all[mask]
        qa_mask  = pixel_data["qa"][mask]

        if len(values_m) == 0:
            warnings.warn("No data to plot")
            return

        smooth_x, smooth_y = f.calculate_spline(
            whole_timeframe=t_all, masked_values=values_m,
            masked_time=time_m, smoothing_parameter=smoothing
        )
        if smooth_x is None:
            warnings.warn("No data to plot")
            return

        metrics_dict, plot_time_frame = f.calculate_metrics_to_plot(
            start=start, end=end, masked_values=values_m,
            masked_time=time_m, smoothing_parameter=smoothing
        )
        if metrics_dict is None:
            return

        # QA-coloured background scatter (unique to this method)
        qa_unique = sorted(np.unique(qa_mask))
        cmap = plt.cm.get_cmap("tab10", len(qa_unique))
        cmap_new = ListedColormap(cmap(np.arange(len(qa_unique))))
        qa_to_idx = {qa: i for i, qa in enumerate(qa_unique)}
        qa_idx = np.array([qa_to_idx[q] for q in qa_mask])
        bounds = np.arange(-0.5, len(qa_unique) + 0.5, 1)

        norm = BoundaryNorm(bounds, cmap_new.N)

        sc = self.plot_background_pts(ax, latitude_idx, longitude_idx, values_m, time_m, aggregation = aggregation)

        # if aggregation:
        #     if self.aggregation_df is None:
        #             self.spatial_aggregation()
        #     background_sub    = self.aggregation_df[(self.aggregation_df["i"] == latitude_idx) & (self.aggregation_df["j"] == longitude_idx)]
        #     background_time   = background_sub["time"].to_numpy()
        #     background_values = background_sub["MA_value"]
        #     sc = ax.scatter(datenum_to_datetime(background_time), background_values, c=qa_idx, cmap=cmap_new, norm=norm, alpha=1, s=10, label="Data")
        # else:
        #     sc = ax.scatter(datenum_to_datetime(time_m), values_m, c=qa_idx, cmap=cmap_new, norm=norm, alpha=1, s=10, label="Data")

        neg_values_sub = f.plot_variables(
            ax=ax, plotting_data=plotting_data, spline_x=smooth_x, spline_y=smooth_y,
            time_frame=plot_time_frame
        )

        self.annotations_and_limits(
            ax=ax, plotting_data=plotting_data, metrics_dict=metrics_dict,
            time_frame=plot_time_frame, lat=lat, lon=lon,
            latitude_idx=latitude_idx, longitude_idx=longitude_idx,
            neg_values_sub=neg_values_sub
        )

        cbar = fig.colorbar(sc, ax=ax, boundaries=bounds)
        cbar.set_label("QA indicators")
        cbar.set_ticks(np.arange(len(qa_unique)))
        cbar.set_ticklabels([str(q) for q in qa_unique])


    def count_extrema(self, latitude_idx, longitude_idx, start= 0, end= 9999, peaks = True):
        """Return the number of detected peaks for a pixel within a year range.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        var = "pks" if peaks else "trgs"
        pixel_data = self._load_pixel_data(latitude_idx, longitude_idx)
        plotting_data = f.grab_plotting_variables(start=start, end=end, pixel_data=pixel_data, variables=[var])
        return len(plotting_data[var][0])


    def create_heatmap_output(self, latitude_idx, longitude_idx, start_year=2002, end_year=2024, fraction = False):
        """Return peak/trough counts or lake-wide fractions per year and quarter.

        When fraction=False, counts are read from the single pixel at (latitude_idx,
        longitude_idx). When fraction=True, all pixels in the lake are aggregated and
        each quarter value is expressed as the fraction of that year's total events
        occurring in that quarter (0.0 – 1.0); years with no events return 0.0.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel. Only used when fraction=False.
        longitude_idx : int
            Column (lon) index of the pixel. Only used when fraction=False.
        start_year : int, optional
            First calendar year to include (inclusive). Default 2002.
        end_year : int, optional
            Last calendar year to include (inclusive). Default 2024.
        fraction : bool, optional
            If False (default), return per-pixel integer counts.
            If True, return lake-wide fractions aggregated across all pixels.

        Returns
        -------
        dict
            Keys are years (int). Values are lists of 4 tuples, one per quarter
            (Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec). Each tuple is
            (n_peaks, n_troughs) when fraction=False, or
            (peaks_fraction, troughs_fraction) when fraction=True.
        """
        quarters = [(1, 3), (4, 6), (7, 9), (10, 12)]
        result = {}
        if not fraction:
            with netCDF4.Dataset(self.p_path) as nc:
                pks_x = f.unix_to_datetime(f.remove_nan(nc.variables["pks_x"][latitude_idx, longitude_idx, :]))
                trgs_x = f.unix_to_datetime(f.remove_nan(nc.variables["trgs_x"][latitude_idx, longitude_idx, :]))
        
                for year in range(start_year, end_year + 1):
                    year_counts = []
                    for (q_start, q_end) in quarters:
                        n_pks = sum(
                            1 for d in pks_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        n_trgs = sum(
                            1 for d in trgs_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        year_counts.append((n_pks, n_trgs))
                    result[year] = year_counts
        else:
            with netCDF4.Dataset(self.p_path) as nc:
                pks_x = f.unix_to_datetime(f.remove_nan(nc.variables["pks_x"][:, :, :]))
                trgs_x = f.unix_to_datetime(f.remove_nan(nc.variables["trgs_x"][:, :, :]))
        
                for year in range(start_year, end_year + 1):
                    year_fractions = []
                    for (q_start, q_end) in quarters:
                        n_pks = sum(
                            1 for d in pks_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        yearly_pks = sum(1 for d in pks_x if d.year == year)
                        pks_fraction = n_pks / yearly_pks if yearly_pks > 0 else 0.0

                        n_trgs = sum(
                            1 for d in trgs_x
                            if d.year == year and q_start <= d.month <= q_end
                        )
                        yearly_trgs = sum(1 for d in trgs_x if d.year == year)
                        trgs_fraction = n_trgs / yearly_trgs if yearly_trgs > 0 else 0.0

                        year_fractions.append((pks_fraction, trgs_fraction))
                    result[year] = year_fractions
        return result
    

    @staticmethod
    def _gap_rectangles(gap_starts, gap_ends):
        """Convert data gap intervals into heatmap cell coordinates.

        For each gap, finds every year/quarter cell it overlaps and returns the
        proportional position of the gap within that cell (matching how axvspan
        renders gaps in single_plot, but mapped onto the heatmap grid).

        Parameters
        ----------
        gap_starts, gap_ends : array-like of timezone-aware datetime

        Returns
        -------
        list of (year, q_idx, x_offset, x_width)
            x_offset and x_width are in [0, 1] relative to the cell width.
        """
        import calendar as _cal
        quarters_def = [(1, 3), (4, 6), (7, 9), (10, 12)]
        rects = []
        for gs, ge in zip(gap_starts, gap_ends):
            for year in range(gs.year, ge.year + 1):
                for q_idx, (m_start, m_end) in enumerate(quarters_def):
                    q_start = datetime.datetime(year, m_start, 1, tzinfo=datetime.timezone.utc)
                    last_day = _cal.monthrange(year, m_end)[1]
                    q_end = datetime.datetime(year, m_end, last_day, 23, 59, 59,
                                             tzinfo=datetime.timezone.utc)
                    q_duration = (q_end - q_start).total_seconds()
                    clipped_start = max(gs, q_start)
                    clipped_end = min(ge, q_end)
                    if clipped_end <= clipped_start or q_duration <= 0:
                        continue
                    x_offset = (clipped_start - q_start).total_seconds() / q_duration
                    x_width = (clipped_end - clipped_start).total_seconds() / q_duration
                    rects.append((year, q_idx, x_offset, x_width))
        return rects


    def yearly_heatmap_pixel(self, latitude_idx, longitude_idx, color_scheme='pink-blue', show_gaps=False):
        """Plot a bivariate heatmap of peak and trough counts or fractions by year and quarter.

        Each cell in the heatmap represents one calendar quarter of one year. The
        cell colour encodes two variables simultaneously using a 4×4 bivariate
        colour palette from color_sets_4x4.

        When whole_lake=False, counts for the single pixel at (latitude_idx, longitude_idx)
        are binned into four levels (0, 1, 2, 3+) and the cell is coloured from the
        discrete 4×4 grid using f.bivariate_legend.

        When whole_lake=True, peak and trough events are aggregated across all
        lake pixels and each cell shows the fraction of that year's total events
        falling in that quarter. Colours are interpolated continuously across the
        4×4 grid using bivariate_continuous_legend(), and the legend is a smooth 2-D gradient
        rendered by bivariate_continuous_legend.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel. Only used when whole_lake=False.
        longitude_idx : int
            Column (lon) index of the pixel. Only used when whole_lake=False.
        color_scheme : str, optional
            Key into color_sets_4x4 selecting the bivariate palette.
            One of 'pink-blue', 'teal-red', 'teal-red1', 'blue-orange'.
            Default 'pink-blue'.
        whole_lake : bool, optional
            If False (default), plot per-pixel counts with a discrete legend.
            If True, plot lake-wide fractions with a continuous gradient legend.
        show_gaps : bool, optional
            If True, overlay a grey strip on the left of each cell whose width
            is proportional to the fraction of that quarter covered by data gaps.
            Default False.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure containing the heatmap.
        ax : matplotlib.axes.Axes
            The axes on which the heatmap is drawn.
        """
        heatmap_data = self.create_heatmap_output(latitude_idx=latitude_idx, longitude_idx=longitude_idx, fraction = False)
        fig, ax, _ = f.create_empty_heatmap()
        g  = self._load_extracted_globals()
        lat, lon = g["lat"], g["lon"]
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])
        textstr = f"Yearly Heatmap for Pixel\n lat, lon: {lat_val:.4f}, {lon_val:.4f}\n {self.variable}, Lake ID:{self.lakeID}"
        ax.set_title(textstr)

        color_set = color_sets_4x4[color_scheme]

        for year, quarters in heatmap_data.items():
            for q_idx, (n_pks, n_trgs) in enumerate(quarters):
                pk_bin  = min(n_pks,  3)
                trg_bin = min(n_trgs, 3)
                color = color_set[trg_bin * 4 + pk_bin]
                ax.add_patch(Rectangle((q_idx, year - 1), 1, 1, facecolor=color, edgecolor='none'))

        if show_gaps:
            px = self._load_pixel_data(latitude_idx, longitude_idx)
            gap_rects = self._gap_rectangles(px["gap_starts"], px["gap_ends"])
            gap_patch_added = False
            for year, q_idx, x_offset, x_width in gap_rects:
                if 2002 <= year <= 2024:
                    label = "Data gap" if not gap_patch_added else None
                    ax.add_patch(Rectangle((q_idx + x_offset, year - 1), x_width, 1,
                                           facecolor='grey', alpha=0.6, edgecolor='none', label=label))
                    gap_patch_added = True
            if gap_patch_added:
                ax.legend(handles=[Patch(facecolor='grey', alpha=0.6, label='Data gap')],
                          bbox_to_anchor=(1.5, 0.5),
                          bbox_transform=ax.transAxes, fontsize=8)

        ax_legend = ax.inset_axes([1.2, 0.7, 0.3, 0.3], transform = ax.transAxes)
        f.bivariate_legend(ax_legend, color_set)

        return fig, ax


    def yearly_heatmap_lake(self, color_scheme="pink-blue", show_gaps=False):
        """Plot a bivariate heatmap of lake-wide peak/trough fractions by year and quarter.

        Parameters
        ----------
        color_scheme : str, optional
            Key into color_sets_4x4 selecting the bivariate palette. Default 'pink-blue'.
        show_gaps : bool, optional
            If True, overlay a grey strip on the left of each cell whose width is
            proportional to the fraction of that quarter covered by data gaps,
            using the first valid pixel as a representative for the whole lake.
            Default False.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ax : matplotlib.axes.Axes
        """
        quarters_def = [(1, 3), (4, 6), (7, 9), (10, 12)]
        # lat and lon are not needed as the heatmap uses all pixels from the lake, thus they can be arbitrary
        heatmap_data = self.create_heatmap_output(latitude_idx=-1, longitude_idx=-1, fraction=True)
        fig, ax, _ = f.create_empty_heatmap()
        textstr = f"Yearly Heatmap for Lake ID: {self.lakeID}\n  {self.variable}"
        ax.set_title(textstr)

        color_set = color_sets_4x4[color_scheme]

        if show_gaps and self.valid_coords:
            rep_i, rep_j = self.valid_coords[0]
            px = self._load_pixel_data(rep_i, rep_j)
            gap_starts = px["gap_starts"]
            gap_ends = px["gap_ends"]
        else:
            gap_starts, gap_ends = [], []

        gap_patch_added = False
        for year, quarters in heatmap_data.items():
            for q_idx, (pks_frac, trgs_frac) in enumerate(quarters):
                color = f.interpolate_from_color_set(pks_frac, trgs_frac, color_set)
                ax.add_patch(Rectangle((q_idx, year - 1), 1, 1, facecolor=color, edgecolor='none'))
                if show_gaps:
                    m_start, m_end = quarters_def[q_idx]
                    gap_frac = self._gap_fraction_in_quarter(gap_starts, gap_ends, year, m_start, m_end)
                    if gap_frac > 0:
                        label = "Data gap" if not gap_patch_added else None
                        ax.add_patch(Rectangle((q_idx, year - 1), gap_frac, 1,
                                               facecolor='grey', alpha=0.6, edgecolor='none', label=label))
                        gap_patch_added = True

        ax_legend = ax.inset_axes([1.2, 0.7, 0.3, 0.3], transform=ax.transAxes)
        f.bivariate_continuous_legend(ax_legend, color_set)

        if show_gaps and gap_patch_added:
            ax.legend(handles=[Patch(facecolor='grey', alpha=0.6, label='Data gap')],
                      loc='lower left', fontsize=8)

        return fig, ax


    def pixel_r2(self, latitude_idx, longitude_idx, start=0, end=9999):
        """Return the R² score for a single pixel within a year range.

        Delegates to r2_scores, which may trigger full-lake parallel computation
        and CSV caching on the first call for this year range.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        return scores[(latitude_idx, longitude_idx)]


    def pixel_rmse(self, latitude_idx, longitude_idx, start=0, end=9999):
        """Return the RMSE for a single pixel within a year range.

        Delegates to RMSE_scores, which may trigger full-lake parallel computation
        and CSV caching on the first call for this year range.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        return scores[(latitude_idx, longitude_idx)]


    def pixel_mad(self, latitude_idx, longitude_idx, start=0, end=9999):
        """Return the Median Absolute Deviation for a single pixel within a year range.

        Delegates to MAD_scores, which may trigger full-lake parallel computation
        and CSV caching on the first call for this year range.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        return scores[(latitude_idx, longitude_idx)]
    

    def pixel_correlation(self, latitude_idx, longitude_idx, start=0, end=9999):
        """Return the Pearson correlation coefficient for a single pixel within a year range.

        Delegates to correlation_scores, which may trigger full-lake parallel
        computation and CSV caching on the first call for this year range.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        return scores[(latitude_idx, longitude_idx)]


    def pixel_values(self, latitude_idx, longitude_idx, start=0, end=9999):
        """Return the valid observation count for a single pixel within a year range.

        Delegates to values_per_pixel, which may trigger full-lake parallel
        computation and CSV caching on the first call for this year range.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        return scores[(latitude_idx, longitude_idx)]
    
    
    def plot_background_pts(self, ax, latitude_idx, longitude_idx, masked_values, masked_time, aggregation = False, **style_kwargs):
        base_style = {"alpha":0.3, "color":"grey","s":10}
        style = {**base_style,**style_kwargs}

        if aggregation:
            if self.aggregation_df is None:
                self.spatial_aggregation()

            background_sub = self.aggregation_df[(self.aggregation_df["i"]==latitude_idx) & (self.aggregation_df["j"]==longitude_idx)]
            background_time = background_sub["time"].to_numpy()
            background_values = background_sub["MA_value"]

            x = f.datenum_to_datetime(background_time)
            y = background_values
        else:
            x = f.datenum_to_datetime(masked_time)
            y = masked_values

        sc = ax.scatter(x, y, label="Data", **style)
        return sc


    def plot_data_gaps(self, ax, pixel_data):
        gap_starts = pixel_data["gap_starts"]
        gap_ends   = pixel_data["gap_ends"]
        for gs, ge in zip(gap_starts, gap_ends):
            ax.axvspan(gs, ge, color="orange", alpha=0.15, zorder=0)
        if len(gap_starts) > 0:
            ax.axvspan(gap_starts[0], gap_ends[0], color="orange", alpha=0.15, zorder=0, label="Data gap")

            
    def plot_pheno_metrics(self, ax, plotting_data, spline_x, spline_y, time_frame, variables = None):

        if variables is None:
            variables = ["pks", "trgs", "midUP", "midDOWN"]

        ax.plot(f.datenum_to_datetime(spline_x), spline_y, color="black", linewidth=1, label="Spline")

        used_labels = set()
        qa_colors = {0: "blue", 1: "orange", 2: "red"}
        qa_labels = {0: "Good", 1: "Fair", 2: "Poor"}
        qa_vars = {'pks','trgs'}
        for var in variables:
            if var in qa_vars:
                qa_vals = plotting_data[var][2]

                for qa in self.QA_LEVELS:
                    pm =  qa_vals == qa
                    if not pm.any():
                        continue
                    label = qa_labels[qa]
                    ax.scatter(plotting_data[var][0][pm], plotting_data[var][1][pm], color=qa_colors[qa], s=50,
                                    marker="o", edgecolors="black", linewidths=0.5,
                                    zorder=4, label=label if label not in used_labels else None)
                    used_labels.add(label)

            else:
                ax.scatter(plotting_data[var][0], plotting_data[var][1], s=30, zorder=4, **self.METRIC_CONFIG[var]['style'])
                used_labels.add(label)
                # labeling not working right now

        neg_values_sub = f.mark_negative_values_timeseries_plot(ax,plotting_data,time_frame)
        return neg_values_sub

    def annotations_and_limits(self, ax, plotting_data, metrics_dict, time_frame, lat_val, lon_val, neg_values_sub, annotation = None):
        start, end = time_frame[0], time_frame[1]

        ax.legend(loc="upper left", ncol= 2)
        textstr = f"{self.variable}\n Lake ID:{self.lakeID}\n lat, lon: {lat_val:.4f}, {lon_val:.4f}\n Total RMSE, R$^2$, MAD: {metrics_dict['rmse'][1]:.4f}, {metrics_dict['r2'][1]:.4f}, {metrics_dict['mad'][1]:.4f}"
        ax.set_title(textstr)
        ax.xaxis.set_minor_locator(mdates.YearLocator())
        ax.grid(axis="x", which="minor", linewidth=0.5)
        ax.grid(axis="x", which="major", linewidth=0.5)
        ax.grid(axis="y")
        ax.set_ylabel("[ug/L]")
        ax.set_xlim(pd.to_datetime('01-01-' + str(start), format='%d-%m-%Y') , pd.to_datetime('31-12-' + str(end), format='%d-%m-%Y'))
        trgs_y_sub = plotting_data["trgs"][1]

        pks_lim_sub = sorted(plotting_data["pks"][1])
        # if max(pks_lim_sub)> 10:
        #         ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-2]+0.5)
        # else:
        #         ax.set_ylim(sorted(trgs_y_sub)[0]-0.5, pks_lim_sub[-1]+0.5)

        trgs_lim_sub = sorted(plotting_data["trgs"][1])

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
                ax.text(0.99,0.99,f"# Neg.values: {sum(neg_values_sub)} \n RMSE: {metrics_dict['rmse'][0]:.3f}\n R$^2$: {metrics_dict['r2'][0]:.3f}\n MAD: {metrics_dict['mad'][0]:.0f}", transform = ax.transAxes,   ha= "right", va= "top", zorder = 10)
            else:
                ax.text(0.99,0.99,f"RMSE:{metrics_dict['rmse'][0]:.3f} \n R$^2$: {metrics_dict['r2'][0]:.3f}\n MAD: {metrics_dict['mad'][0]:.3f}", transform = ax.transAxes,   ha= "right", va= "top", zorder = 10)
        else:
            lines = []
            if "R2" in annotation:
                lines.append(f"R$^2$: {round(metrics_dict['r2'][0], 3)}")
            if "RMSE" in annotation:
                lines.append(f"RMSE: {round(metrics_dict['rmse'][0], 3)}")
            if "MAD" in annotation:
                lines.append(f"MAD: {round(metrics_dict['mad'][0], 3)}")
            if "neg" in annotation and neg_values_sub:
                lines.append(f"# Neg.values: {sum(neg_values_sub)}")
            if lines:
                ax.text(0.99, 0.99, "\n".join(lines),
                    transform=ax.transAxes, ha="right", va="top", zorder = 10)


    def single_plot(self, latitude_idx, longitude_idx, ax, aggregation = False, start= 0, end= 9999, annotation = None,variables = None):
        """Plot raw observations, the smoothed spline, and all phenological events for a pixel.

        Displays a scatter of valid (QA==0) observations or 3×3 aggregated values,
        overlaid with the csaps spline and scatter markers for peaks, troughs,
        green-up midpoints, and green-down midpoints. Negative values are flagged
        with red crosses and trigger a warning. Fit metrics (RMSE, R², MAD) are
        annotated on the plot.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        pixel_data = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        lat_val = float(lat[latitude_idx])
        lon_val = float(lon[longitude_idx])
        smoothing = pixel_data["smoothing"]

        plotting_data = f.grab_plotting_variables(start = start, end = end, pixel_data=pixel_data, variables=variables)

        mask     = (pixel_data["values"] != -9999) & (pixel_data["qa"] == 0)
        values_m = pixel_data["values"][mask]
        time_m   = t_all[mask]

        if len(values_m) == 0:
            warnings.warn("No data to plot")
            return

        smooth_x, smooth_y = f.calculate_spline(whole_timeframe= t_all, masked_values=values_m, masked_time= time_m, smoothing_parameter=smoothing)

        metrics_dict, plot_time_frame = f.calculate_metrics_to_plot(start = start, end = end, masked_values= values_m, masked_time=time_m, smoothing_parameter=smoothing)

        if metrics_dict is None:
            return

        self.plot_background_pts(ax = ax, latitude_idx= latitude_idx, longitude_idx = longitude_idx, masked_values=values_m, masked_time=time_m, aggregation=aggregation)
        self.plot_data_gaps(ax = ax, pixel_data = pixel_data)
        neg_values_sub = f.plot_variables(ax = ax, plotting_data= plotting_data, spline_x= smooth_x, spline_y= smooth_y, time_frame= plot_time_frame, variables= variables)
        self.annotations_and_limits(ax = ax, plotting_data= plotting_data, metrics_dict= metrics_dict, time_frame=plot_time_frame, lat_val = lat_val, lon_val = lon_val, neg_values_sub=neg_values_sub, annotation = annotation)


    def split_plot(self, latitude_idx, longitude_idx, ax0, ax1, aggregation = False, start0= 0, end0= 9999, start1= 0, end1=9999):
        """Plot two year-windowed single_plots side by side for the same pixel.

        Calls single_plot twice — once on ax0 with [start0, end0] and once on ax1
        with [start1, end1]. Intended for comparing two non-overlapping time periods
        at the same pixel.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
        if (start0 == 0) and (start1 == 0) and (end0 == 9999) and (end1 == 9999):
            warnings.warn("split_plot needs a least end0 and start1 parameter, otherwise use full_plot")

        self.single_plot(latitude_idx = latitude_idx, longitude_idx= longitude_idx, ax=ax0, aggregation = aggregation, start= start0, end = end0)
        self.single_plot(latitude_idx = latitude_idx, longitude_idx= longitude_idx, ax=ax1, aggregation = aggregation, start=start1, end=end1)

    def full_plot(self, latitude_idx, longitude_idx, ax, aggregation = False):
        """Plot the complete valid time series for a pixel.

        Auto-detects the first and last years with valid observations and
        delegates to single_plot with those bounds.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
            t_all = f.unix_to_datenum(nc.variables["time"])
            variable = getattr(nc, "variable")
            values = np.array(nc.variables[variable][:, latitude_idx, longitude_idx])
            mask = (values != -9999) & (np.array(nc.variables[getattr(nc, 'qa')][:, latitude_idx, longitude_idx]) == 0)
            values_m = values[mask]
            time_m = t_all[mask]

        if len(values_m) > 1:
            limits = sorted(f.datenum_to_datetime(time_m))
            full_plot_start = min(limits).year
            full_plot_end = max(limits).year
            self.single_plot(latitude_idx= latitude_idx, longitude_idx= longitude_idx, ax = ax, aggregation = aggregation, start=full_plot_start, end= full_plot_end)

        else:
            warnings.warn("No data to plot")


    def single_years_plot(self, latitude_idx, longitude_idx, years, ncol, nrow, annotation, ylim=None):
        """Plot one panel per year in a grid, each showing phenology for a single pixel.

        Creates a figure with ``nrow × ncol`` subplots. Each subplot calls
        :meth:`single_plot` for one year in *years*, with enlarged markers and
        month-number x-axis labels. Panels beyond ``len(years)`` are hidden.

        Parameters
        ----------
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
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
            self.single_plot(latitude_idx, longitude_idx, ax, start=year, end=year, annotation=annotation)

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


    def single_plot_insitu(self, latitude_idx,longitude_idx,ax,insitu_df,aggregation=False, start = 0, end = 9999, insitu_date_col="datetime", insitu_value_col="chlorophyll_a", insitu_station_col=None, station_id=None, max_depth = 5):
        """
        Plot satellite observations + spline + phenology + in situ overlay.
        """

        # -------------------------------------------------
        # FIRST: draw the original plot
        # -------------------------------------------------

        self.single_plot(
            latitude_idx=latitude_idx,
            longitude_idx=longitude_idx,
            ax=ax,
            aggregation=aggregation,
            start=start,
            end=end)

        # -------------------------------------------------
        # PREPARE IN SITU DATA
        # -------------------------------------------------

        insitu_mean = f.prep_dimark_data(insitu_df = insitu_df,start=start, end=end,  insitu_date_col=insitu_date_col, insitu_value_col=insitu_value_col, insitu_station_col=insitu_station_col, station_id=station_id, max_depth = max_depth)
        
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


    def plot_background_ratio_timeseries(self, other, ax, latitude_idx, longitude_idx, color="blue"):
        """Plot the per-observation ratio between two instances at a single pixel.

        For each timestep where both instances have a valid QA-0 observation,
        computes self / other and renders the result as a scatter plot on ax.
        Observations where other equals zero are excluded to avoid division by zero.
        Both instances must refer to the same lake; a warning is raised otherwise.

        Parameters
        ----------
        other : PhenologyVisualization
            Second instance whose values form the denominator of the ratio.
        ax : matplotlib.axes.Axes
            Axes on which to draw the scatter plot.
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        color : str, optional
            Marker colour. Default is "blue".

        Returns
        -------
        None
        """
        g_self  = self._load_extracted_globals()
        g_other = other._load_extracted_globals()

        pixel_data = self._load_pixel_data(latitude_idx, longitude_idx)

        t_self  = g_self["t_all"]
        t_other = g_other["t_all"]

        lakeID1 = self.lakeID
        lakeID2 = other.lakeID
        if lakeID1 != lakeID2:
                raise Warning("Comparison must be made on the same lake!")

        phenology_name1 = self.variable
        phenology_name2 = other.variable

        pixel_self  = self._load_pixel_data(latitude_idx, longitude_idx)
        pixel_other = other._load_pixel_data(latitude_idx, longitude_idx)

        mask_self  = (pixel_self["values"]  != -9999) & (pixel_self["qa"]  == 0)
        mask_other = (pixel_other["values"] != -9999) & (pixel_other["qa"] == 0)

        df_self  = pd.DataFrame({"time": t_self[mask_self],   "value": pixel_self["values"][mask_self]})
        df_other = pd.DataFrame({"time": t_other[mask_other], "value": pixel_other["values"][mask_other]})

        merged = df_self.merge(df_other, on="time", suffixes=("_self", "_other"), how="inner")
        merged = merged[merged["value_other"] != 0]

        ratio = merged["value_self"] / merged["value_other"]

        ax.scatter(
            f.datenum_to_datetime(merged["time"].to_numpy()),
            ratio,
            color=color,
            s=10,
            label="Ratio")
        self.plot_data_gaps(ax = ax, pixel_data = pixel_data)

        ax.legend(loc = "upper left")
        ax.set_ylabel(f"{phenology_name1}/ {phenology_name2}")
        ax.set_title(f"Background Points Ratio Lake ID: {lakeID1}")
        ax.xaxis.set_minor_locator(mdates.YearLocator())
        ax.grid(axis="x", which="minor", linewidth=0.5)
        ax.grid(axis="x", which="major", linewidth=0.5)
        ax.grid(axis="y")


    def plot_background_ratio_v_self(self, other, ax, latitude_idx, longitude_idx, color="blue"):
        """Plot the per-observation ratio between two instances at a single pixel.

        For each timestep where both instances have a valid QA-0 observation,
        computes self / other and renders the result as a scatter plot on ax.
        Observations where other equals zero are excluded to avoid division by zero.
        Both instances must refer to the same lake; a warning is raised otherwise.

        Parameters
        ----------
        other : PhenologyVisualization
            Second instance whose values form the denominator of the ratio.
        ax : matplotlib.axes.Axes
            Axes on which to draw the scatter plot.
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        color : str, optional
            Marker colour. Default is "blue".

        Returns
        -------
        None
        """
        g_self  = self._load_extracted_globals()
        g_other = other._load_extracted_globals()

        t_self  = g_self["t_all"]
        t_other = g_other["t_all"]

        lakeID1 = self.lakeID
        lakeID2 = other.lakeID
        if lakeID1 != lakeID2:
                raise Warning("Comparison must be made on the same lake!")

        phenology_name1 = self.variable
        phenology_name2 = other.variable

        pixel_self  = self._load_pixel_data(latitude_idx, longitude_idx)
        pixel_other = other._load_pixel_data(latitude_idx, longitude_idx)

        mask_self  = (pixel_self["values"]  != -9999) & (pixel_self["qa"]  == 0)
        mask_other = (pixel_other["values"] != -9999) & (pixel_other["qa"] == 0)

        df_self  = pd.DataFrame({"time": t_self[mask_self],   "value": pixel_self["values"][mask_self]})
        df_other = pd.DataFrame({"time": t_other[mask_other], "value": pixel_other["values"][mask_other]})

        merged = df_self.merge(df_other, on="time", suffixes=("_self", "_other"), how="inner")
        merged = merged[merged["value_other"] != 0]

        ratio = merged["value_self"] / merged["value_other"]
        
        ax.scatter(
            merged["value_other"].values,
            ratio,
            color=color,
            s=10,
            label="Ratio")

        ax.set_ylabel(f"{phenology_name1}/ {phenology_name2}")
        ax.set_xlabel(f"{phenology_name2}")
        ax.set_title(f"Background Points Ratio Lake ID: {lakeID1}")
        ax.grid(axis="x", which="minor", linewidth=0.5)
        ax.grid(axis="x", which="major", linewidth=0.5)
        ax.grid(axis="y")


    def yearly_cubic_spline(self, ax, latitude_idx, longitude_idx, years=None):
        """Overlay csaps splines for multiple years on a common fractional-month x-axis.

        Fits a single spline over the full valid time series for the pixel at
        (latitude_idx, longitude_idx), then slices it year by year and plots each slice
        against fractional month (1–12) using a distinct colour from the
        cc.glasbey_light palette. Detected peaks and troughs are overlaid as
        scatter markers coloured by QA level.

        A vertical colorbar maps year indices to their assigned colours, and two
        separate legends show QA marker styles (Good / Fair / Poor) and event
        types (Peak / Trough).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes on which to draw the overlaid splines and markers.
        latitude_idx : int
            Row (lat) index of the pixel.
        longitude_idx : int
            Column (lon) index of the pixel.
        years : list of str, optional
            Calendar years to include. Each element must be a string (e.g. '2005').
            Defaults to 2002–2024. Years with fewer than 3 spline points are
            skipped with a warning.

        Returns
        -------
        None
        """
        if years is None:
            years = range(2002, 2025)

        g = self._load_extracted_globals()
        px = self._load_pixel_data(latitude_idx, longitude_idx)
        lat, lon, t_all = g["lat"], g["lon"], g["t_all"]
        smoothing = px["smoothing"]
        lat_val = round(float(lat[latitude_idx]), 4)
        lon_val = round(float(lon[longitude_idx]), 4)

        # QA legend (markers only, no year lines)
        qa_markers = {0: "o", 1: ".", 2: "x"}
        qa_labels = {0: "Good", 1: "Fair", 2: "Poor"}
        qa_handles = [
            mlines.Line2D(
                [], [],
                color="black",
                marker=qa_markers[qa],
                linestyle="None",
                markersize=6,
                label=qa_labels[qa]
            )
            for qa in self.QA_LEVELS
        ]
        
        cmap = cc.cm.rainbow
        year_colors = {year: cmap(i / max(len(years) - 1, 1)) for i, year in enumerate(years)}
        cmap_discrete = ListedColormap([year_colors[y] for y in years])

        # Fit one spline over all valid data
        mask_all = (px["values"] != -9999) & (px["qa"] == 0)
        values_m_all = px["values"][mask_all]
        time_m_all = t_all[mask_all]
        smooth_x_all, smooth_y_all = f.calculate_spline(
            whole_timeframe=t_all, masked_values=values_m_all,
            masked_time=time_m_all, smoothing_parameter=smoothing
        )
        if smooth_x_all is None:
            warnings.warn("No data to plot")
            return
        smooth_dates_all = np.array(f.datenum_to_datetime(smooth_x_all))

        for year in years:
            plotting_data = f.grab_plotting_variables(
                start=year,
                end=year,
                pixel_data=px,
                variables=["pks", "trgs"]
            )
            pks_x_sub, pks_y_sub, pks_qa_sub   = plotting_data["pks"]
            trgs_x_sub, trgs_y_sub, trgs_qa_sub = plotting_data["trgs"]

            # Subset the pre-fitted spline to this year and convert to fractional month (1–12)
            mask_year = np.array([d.year == year for d in smooth_dates_all])
            smooth_dates_year = smooth_dates_all[mask_year]
            smooth_y  = smooth_y_all[mask_year]

            if len(smooth_dates_year) > 1:
                smooth_x_month = f.to_frac_month(smooth_dates_year)
                ax.plot(smooth_x_month, smooth_y, color=year_colors[year], linewidth=1, label=str(year))

                for qa in self.QA_LEVELS:
                    pm = pks_qa_sub == qa
                    tm = trgs_qa_sub == qa
                    if pm.any():
                        ax.scatter(f.to_frac_month(pks_x_sub[pm]), pks_y_sub[pm], color="black", s=50,
                            marker=qa_markers[qa], edgecolors="black", linewidths=0.5,
                            zorder=4, label= qa_labels[qa] if year == years[0] else None)
                    if tm.any():
                        ax.scatter(f.to_frac_month(trgs_x_sub[tm]), trgs_y_sub[tm], color="darkgray", s=50,
                            marker=qa_markers[qa], edgecolors="black", linewidths=0.5,
                            zorder=4, label=qa_labels[qa] if (year == years[0] and not pm.any()) else None)
            else:
                warnings.warn(f"Not enough data to plot for year {year}")

        # One-time axis setup after all years are plotted
        int_years = [int(y) for y in years]
        # derive months from data
        all_months = sorted(set(d.month for d in smooth_dates_all if d.year in int_years))
        month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
        ax.set_xlim(min(all_months), max(all_months))
        ax.set_xticks(all_months)
        ax.set_xticklabels([month_names[m - 1] for m in all_months])
        all_pks_y = px["pks_y"][np.array([d.year in int_years for d in px["pks_x"]])]
        all_trgs_y = px["trgs_y"][np.array([d.year in int_years for d in px["trgs_x"]])]
        if len(all_pks_y) > 0 and len(all_trgs_y) > 0:
            pks_sorted = sorted(all_pks_y)
            ymax = pks_sorted[-2] if pks_sorted[-1] > 10 and len(pks_sorted) > 1 else pks_sorted[-1]
            ax.set_ylim(sorted(all_trgs_y)[0] - 0.5, ymax + 0.5)
        textstr = f"{self.variable}\n Lake ID:{self.lakeID}\n lat, lon: {lat_val}, {lon_val}"
        ax.set_title(textstr)
        ax.grid(axis="x", linewidth=0.5)
        ax.grid(axis="y")
        ax.set_ylabel("[ug/L]")

        # Color legend for pks and trgs
        type_handles = [
            mlines.Line2D([], [], color="black", marker="o",
                        linestyle="None", markersize=8, label="Peak"),
            mlines.Line2D([], [], color="darkgray", marker="o",
                        linestyle="None", markersize=8, label="Trough"),
        ]

        leg1 = ax.legend(handles=qa_handles, loc="upper left")
        ax.add_artist(leg1)
        ax.legend(handles=type_handles, loc="upper right")

        # Colorbar for year colors
        norm = mcolors.BoundaryNorm(boundaries=range(len(years) + 1), ncolors=len(years))
        sm   = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation="vertical", pad=0.01, aspect=30)
        cbar.set_ticks([i + 0.5 for i in range(len(years))])
        cbar.set_ticklabels(years)


    def extract_pixel_timing_OHE(self, nc, i, j, vars= None, qa_var = 'pks'):
        frames = {}
        if vars is None:
            vars = ['pks','trgs','midUP','midDOWN']

        for var in vars:
            var = f.coerce_varname_to_var_x(var)
    
            var_raw = f.remove_nan(nc.variables[var][i,j,:])
            if len(var_raw) <1:
                continue
            var_dt = pd.to_datetime(var_raw, unit='s',utc=True)
            doy = var_dt.year + var_dt.day_of_year / 1000
            frames[var] = pd.Series(index=doy,dtype=bool)

        qa_var = f.parse_qa_var_from_str(qa_var)
        if qa_var is not None:
            qa_x = f.coerce_varname_to_var_x(qa_var)
            qa_x_raw = np.array(nc.variables[qa_x][i, j, :])
            qa_mask = ~np.isnan(qa_x_raw)
            if len(qa_x_raw) <1:
                pass
            qa_dt = pd.to_datetime(qa_x_raw[qa_mask], unit="s", utc=True)
            qa_doy = qa_dt.year + qa_dt.day_of_year / 1000
            qa = np.array(nc.variables[qa_var][i, j, :])[qa_mask].astype(int)
            frames[qa_var] = pd.Series(index=qa_doy, data= qa)
        return pd.DataFrame(frames)

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
        return PhenologyVisualization._extract_pixel_kde_events(
            None, None, i, j, primary_vars, secondary_vars, qa_var, arrays=arrays
        )

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


    def prep_kde_data(self, kde_df):
        """Pair each peak with its bracketing green-up advanced and green-down onset, per pixel.

        For each pixel, finds all peaks and matches each one with the most recent
        green-up advanced event before it and the earliest green-down onset event
        after it (within that same pixel). Peaks that cannot be fully bracketed are
        dropped. Results from all pixels are then pooled into a single DataFrame.

        Parameters
        ----------
        kde_df : pandas.DataFrame
            Output of assemble_kde_data. Index is year.DOY (float), columns include
            i, j (pixel indices), green_up_advanced (bool), peaks (float QA or NaN),
            green_down_onset (bool).

        Returns
        -------
        pandas.DataFrame
            Columns:
              green_up_advanced : float — year.DOY of the bracketing green-up advanced event.
              peak_qa           : int   — QA indicator of the peak (0=Good, 1=Fair, 2=Poor).
              green_down_onset  : float — year.DOY of the bracketing green-down onset event.
            Each row is one fully-bracketed peak. Duplicates across pixels are retained.
        """

        print(f"prep kde started at: {datetime.datetime.now()}")

        pixel_groups = [pixel_df for _, pixel_df in kde_df.groupby(["i", "j"])]

        with multiprocessing.Pool(processes=os.cpu_count()) as pool:
            results = pool.map(PhenologyVisualization._prep_kde_pixel_worker, pixel_groups)

        frames = [df for df in results if not df.empty]
        if not frames:
            return pd.DataFrame(columns=["primary", "qa_column", "secondary"])
        print(f"prep kde finished at: {datetime.datetime.now()}")
        return pd.concat(frames, ignore_index=True)
    

    @staticmethod
    def _prep_kde_pixel_worker(pixel_df):
        """Multiprocessing worker for a single pixel's bracket-matching step.

        Designed to be dispatched by multiprocessing.Pool in prep_kde_data.
        No initializer needed — all data is passed directly as a DataFrame.

        Parameters
        ----------
        pixel_df : pandas.DataFrame
            Single-pixel slice of the kde_df produced by assemble_kde_data.

        Returns
        -------
        pandas.DataFrame
            Columns: primary, qa_column, secondary. Empty if no fully-bracketed
            peaks exist for this pixel.
        """
        bracketed = PhenologyVisualization.find_preceding_and_following(
            pixel_df.sort_index()
        ).dropna(subset=["prev_var1_time", "next_var2_time"])
        if bracketed.empty:
            return pd.DataFrame(columns=["primary", "qa_column", "secondary"])
        bracketed = bracketed.rename(columns={
            "prev_var1_time": "primary",
            "next_var2_time": "secondary",
        })
        bracketed["qa_column"] = bracketed["qa_column"].astype(int)
        return bracketed[["primary", "qa_column", "secondary"]]

    @staticmethod
    def find_preceding_and_following(pixel_df):
        """Vectorized bracket-finder for a single pixel's long-format event DataFrame.

        For each peak row (qa_column not NaN), finds the most recent preceding
        green_up_advanced event (primary=True) and the earliest following
        green_down_onset event (secondary=True) using ffill/bfill on the year.DOY index.

        Parameters
        ----------
        pixel_df : pandas.DataFrame
            Single-pixel slice of the output of assemble_kde_data. Index is year.DOY,
            columns include primary (bool), qa_column (float), secondary (bool).

        Returns
        -------
        pandas.DataFrame
            Columns: prev_var1_time, next_var2_time, qa_column.
            Only peak rows are returned; NaN in either bracket column means
            the peak is not fully bracketed.
        """
        df = pixel_df.copy()
        df["var1_time"] = df.index.where(df["primary"].fillna(False).astype(bool))
        df["var2_time"] = df.index.where(df["secondary"].fillna(False).astype(bool))
        df["prev_var1_time"] = df["var1_time"].ffill()
        df["next_var2_time"] = df["var2_time"].bfill()
        mask = df["qa_column"].notna()
        return df.loc[mask, ["prev_var1_time", "next_var2_time", "qa_column"]]


    def sort_by_year(self, df, start_year=None, end_year=None):
        """Filter a prep_kde_data DataFrame to rows whose bloom overlaps [start_year, end_year].

        Each row is treated as a triple (green_up_advanced, peak_qa, green_down_onset).
        A row is kept when its bloom interval [adv_year, onset_year] overlaps the
        requested range. Cross-year blooms (e.g. green_up_advanced in December,
        green_down_onset in January of the next year) are included if either end
        falls within the range — the whole triple is kept rather than dropped.

        Parameters
        ----------
        df : pandas.DataFrame
            Output of prep_kde_data (columns: green_up_advanced, peak_qa, green_down_onset).
        start_year : int or None
            Earliest year to include. None means no lower bound.
        end_year : int or None
            Latest year to include. None means no upper bound.

        Returns
        -------
        pandas.DataFrame
            Filtered copy of df.
        """
        if start_year is None and end_year is None:
            return df

        adv_year = df["primary"].astype(int)
        onset_year = df["secondary"].astype(int)

        mask = pd.Series(True, index=df.index)
        if start_year is not None:
            # keep rows where the bloom ends at or after start_year
            mask &= onset_year >= start_year
        if end_year is not None:
            # keep rows where the bloom begins at or before end_year
            mask &= adv_year <= end_year

        return df[mask]


    def lake_bloom_kde(self, ax, qa_value = None, start_year = 0, end_year = 9999):
        print(f"plotting started at: {datetime.datetime.now()}")
        dir_path, file_path = self.build_kde_path()
        if os.path.isfile(file_path):
            compressed_df = pd.read_csv(file_path)
            print(file_path)
        else:
            warnings.warn("KDE events need to be calculated. Depending on the lake size this may take a while.")
            os.makedirs(dir_path, exist_ok=True)
            df = self.assemble_kde_data()
            compressed_df = self.prep_kde_data(df)
            compressed_df.to_csv(file_path, index=False)

        if qa_value is not None:
            if type(qa_value) != set:
                warnings.warn("qa_value needs to be a set")
            else:
                compressed_df = compressed_df[compressed_df['qa_column'].isin(qa_value)]

        
        if len(compressed_df) < 2:
            warnings.warn("Not enough data to plot kde")
            return
        
        years_all = np.unique(list(compressed_df["primary"].astype(int) ) + list(compressed_df["secondary"].astype(int) ))
        start, end = f.define_year_range(start_year, end_year, years_all)
        plot_df = self.sort_by_year(compressed_df, start_year=start, end_year=end)
        


        x = np.round((plot_df["primary"].values % 1) * 1000).astype(int)
        y = np.round((plot_df["secondary"].values % 1) * 1000).astype(int)

        kde = gaussian_kde(np.vstack([x, y]))
        xi = np.linspace(0, 400, 100)
        yi = np.linspace(0, 400, 100)
        Xi, Yi = np.meshgrid(xi, yi)
        Zi = kde(np.vstack([Xi.ravel(), Yi.ravel()])).reshape(Xi.shape)
        Zi_norm = Zi / Zi.max()
        cf = ax.contourf(Xi, Yi, Zi_norm, levels=np.linspace(0.05, 1.0, 20), cmap="viridis")
        plt.colorbar(cf, ax=ax)
        ax.axline((0, 0), slope=1, color="black", linewidth=1, linestyle="--")
        ax.axline((0, 365), slope=1, color="black", linewidth=1, linestyle="--")

        ax.set_xlim(0, 400)
        ax.set_ylim(0, 400)
        ax.set_xlabel("Green-up Advanced (DOY)")
        ax.set_ylabel("Green-down Onset (DOY)")
        ax.set_title(f"Green-up Advanced vs Green-down Onset\nLake ID: {self.lakeID} | {start} - {end}")
        print(f"plotting ended at: {datetime.datetime.now()}")


            




