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


class TimingVisualization(PhenologyVisualizationBase):
    """Graphs describing event timing distributions and probabilities."""

    graph_type = "timing"
    graph_methods = frozenset({
        "plot_lake_bloom_kde",
        "lake_bloom_kde_from_fp",
    })

    def plot_lake_bloom_kde(
        self,
        ax,
        probs,
        metadata,
        plt_kwargs=None,
    ):
        """
        Plot bloom probability surface from an already-loaded probability table.
        """

        interval = metadata["interval"]
        start = metadata["start_year"]
        end = metadata["end_year"]
        qa_label = metadata["qa_label"]
        x_max = metadata["x_max"]
        y_max = metadata["y_max"]

        qa_filtered_set = (
                None
                if qa_label == "all"
                else {int(value) for value in qa_label.removeprefix("qa")}
            )

        if plt_kwargs is None:
            plt_kwargs = {
                "cmap": "ocean_r",
                "levels": 20,
                "cbar": True,
            }

        plt_kwargs = dict(plt_kwargs)
        show_cbar = plt_kwargs.pop("cbar", False)

        grid = probs.pivot(index="y_low", columns="x_low", values="probability")

        x_centers = grid.columns.values + interval / 2
        y_centers = grid.index.values + interval / 2

        cf = ax.contourf(x_centers, y_centers, grid.values, **plt_kwargs, vmax=0.05)

        ax.contour(x_centers, y_centers, grid.values, levels=[0.05], colors="red", linewidths=2)

        if show_cbar:
            cbar = plt.colorbar(cf, ax=ax, label="Probability")
            cbar.ax.axhline(0.05, color="red", linewidth=2)

        if qa_filtered_set is None:
            qa_label = "All QA"
        else:
            qa_label = (
                "QA: "
                + ", ".join(
                    self.get_plot_config("qa", q)["label"]
                    for q in sorted(qa_filtered_set)
                )
            )

        ax.set_title(
            f"Bloom event probability per {interval}-day interval\n"
            f"Lake ID: {self.lakeID} | {start} - {end} | {qa_label}"
        )

        ax.set_xlabel("Green-up Advanced (DOY)")
        ax.set_ylabel("Green-down Onset (DOY)")

        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_max)

        ax.axline((0, 0), slope=1, color="black", linewidth=1, linestyle="--")
        ax.axline((0, 365), slope=1, color="black", linewidth=1, linestyle="--")

        return ax

    def lake_bloom_kde_from_fp(self,ax,file_path,plt_kwargs=None):
        """
        Read a cached BloomProb file and plot the probability surface.
        """
        
        probs, metadata = f.read_bloom_probability_with_metadata(file_path)

        return self.plot_lake_bloom_kde(
            ax=ax,
            probs=probs,
            metadata=metadata,
            plt_kwargs=plt_kwargs,
        )
