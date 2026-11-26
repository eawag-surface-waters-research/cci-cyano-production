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
import re
from pathlib import Path
import csv
import statistics
import warnings
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, Patch
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr, gaussian_kde, mannwhitneyu
from csaps import csaps
import functions as f
import multiprocessing
from functools import partial
import shapely.ops as ops
from pyproj import CRS, Transformer
import geopandas
from shapely.prepared import prep
from shapely.geometry import Point
from shapely import vectorized
from numpy.lib.stride_tricks import sliding_window_view
import colorcet as cc
import time
from visualization.base import PhenologyVisualizationBase


class ComparisonVisualization(PhenologyVisualizationBase):
    """Graphs comparing phenology products or variables."""

    graph_type = "comparison"
    graph_methods = frozenset({
        "extrema_comparison",
        "plot_background_ratio_timeseries",
        "plot_background_ratio_v_self",
    })
    
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


    def extrema_comparison(self, other1,  latitude_idx, longitude_idx, ax,  peak = True,
                           aggregation= False, start = 0, end = 9999, background_pts = True, 
                           other2= None, purple_chla21= False, show_legend= False):
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

