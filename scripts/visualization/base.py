import pandas as pd
import netCDF4
import matplotlib
matplotlib.use("Agg")
import numpy as np
from pathlib import Path
import functions as f
import shapely.ops as ops
from pyproj import CRS, Transformer
import geopandas as gpd
from shapely.prepared import prep
from shapely import vectorized


class PhenologyVisualizationBase:
    """Shared data access and plotting primitives for phenology graphs.

    The graph-family subclasses at the end of this module provide focused
    entry points while this base owns the lake metadata, caches, and reusable
    plotting helpers.
    """
    shapefile_path = None
    save_format = "netcdf"  # "csv" or "netcdf" -- see set_save_format();
    # applies to spatial_aggregation() and the metric caches
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
            self.gdf = gpd.read_file(self.shapefile_path)
        self.p_path = phenology_path
        self.e_path = extract_path
        self.version = Path(self.p_path).parents[2].name.removeprefix('v')
        self.variable = Path(self.p_path).parents[0].stem
        self.lakeID = Path(self.p_path).stem
        self.info = (f"Version: {self.version} \n",
                    f"Variable: {self.variable} \n" ,
                    f"Lake ID: {self.lakeID}")
        self.methods = [method for method in dir(type(self)) if callable(getattr(type(self), method)) and not method.startswith("__")]
        self.valid_coords = f.valid_index_pairs(self.e_path)
        self.out_folder = Path(self.p_path).parents[2]
        self.data_folder = Path(self.p_path).parents[2]
        self.start_year , self.end_year = f.get_year_edges(self.e_path)
        self.aggregation_df = None
        self.aggregation_ds = None
        self._aggregation_pixel_index = None
        self.geom_shrunk = None
        self._extracted_globals = None
        self._pixel_cache = {}
        self._lake_cache = {}
        self.prep_geometry_from_shapefile()
        self.extract_nonborder_coords()
    
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


    @classmethod
    def set_save_format(cls, fmt: str):
        """Set the on-disk cache format used by spatial_aggregation() and the
        per-pixel metric caches (r2, MAD, RMSE, correlation, values_per_pixel)
        for all instances.

        Parameters
        ----------
        fmt : str
            "csv" (default): spatial_aggregation() writes a long-format table,
            one row per timestep x pixel, with lat/lon repeated on every row;
            metric caches write one row per pixel.
            "netcdf": spatial_aggregation() writes a compressed (time, pixel)
            array, lat/lon stored once per pixel, chunked for fast single-pixel
            reads; metric caches write a dense (lat, lon) grid. Much smaller on
            disk and faster to load for large lakes.
        """
        if fmt not in ("csv", "netcdf"):
            raise ValueError(f"save_format must be 'csv' or 'netcdf', got {fmt!r}")
        cls.save_format = fmt
            

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
    

    def extract_nonborder_coords(self):
        with netCDF4.Dataset(self.e_path) as nc:
            lats   = nc.variables["lat"][:]
            lons   = nc.variables["lon"][:]
        idxs = np.array(self.valid_coords)
        ii = idxs[:, 0]
        jj = idxs[:, 1]

        lons = lons[jj]  # longitudes
        lats = lats[ii]  # latitudes

        # mask of points inside geometry
        mask = vectorized.contains(self.prepped_geom.context, lons, lats)

        
        inside_indices = [
            (int(i), int(j))
            for i, j in zip(ii[mask], jj[mask])
        ]

        # collect indices + coordinates
        inside_coords = list(zip(lats[mask], lons[mask]))

        self.valid_idx_prep = inside_indices
        self.valid_coords_prep = inside_coords


    # MOVE TO PixelCalc
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
            return self.compute_and_cache_metric(metric_name="r2", col_name="r2_scores",
                                                 compute_fn=PhenologyVisualizationBase.compute_metric_score,
                                                 start=start, end=end)


    # MOVE TO PixelCalc
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
            return self.compute_and_cache_metric(metric_name="MAD", col_name="mad_scores",
                                                 compute_fn=PhenologyVisualizationBase.compute_metric_score,
                                                 start= start,end= end)


    # MOVE TO PixelCalc
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
            return self.compute_and_cache_metric(metric_name="RMSE", col_name="rmse_scores",
                                                 compute_fn=PhenologyVisualizationBase.compute_metric_score,
                                                 start=start, end=end)


    # MOVE TO PixelCalc
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
            return self.compute_and_cache_metric(metric_name="correlation", col_name="correlation_scores",
                                                 compute_fn=PhenologyVisualizationBase.compute_metric_score,
                                                 start=start, end=end)


    # MOVE TO PixelCalc
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
            return self.compute_and_cache_metric(metric_name="values_per_pixel", col_name="number_of_values",
                                                 compute_fn=PhenologyVisualizationBase.compute_metric_score,
                                                 start=start,end= end)


    def _load_extracted_globals(self, reload=False):
        """Lazily load and cache shared arrays from the extract dataset.

        On the first call, opens the extract NetCDF and stores lat, lon, the full
        time array as datenums, and the variable and QA attribute names. Subsequent
        calls return the cached dict without reopening the file.

        Returns
        -------
        dict
            Dictionary with keys: 'lat', 'lon', 't_all', 'variable', 'qa'.
        """
        if (self._extracted_globals is None) or (reload):
            self._extracted_globals = f.load_extract_meta_from_path(self.e_path)
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
            extract_vals = {
                "values": np.array(g["variable"][:, i, j]),
                "qa": np.array(g["qa"][:, i, j])
            }

            self._pixel_cache[(i,j)] = extract_vals | f.load_pixel_data_from_path(self.e_path,self.p_path)
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



    # UNSURE WHERE TO PLACE
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

    # UNSURE WHERE TO PLACE
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


    # UNSURE WHERE TO PLACE
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
    

    # UNSURE WHERE TO PLACE
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


    # UNSURE WHERE TO PLACE
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
    