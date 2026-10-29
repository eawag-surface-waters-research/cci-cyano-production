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

    
class PhenologyVisualization(
    HeatmapVisualization,
    BoxplotVisualization,
    PixelVisualization,
    MapVisualization,
    ComparisonVisualization,
    TimingVisualization,
):
    """Backward-compatible composite of all phenology graph families.

    New code can use a focused subclass such as :class:`MapVisualization`.
    The composite remains the default because the analysis pipeline creates
    one object and passes it to several graph-producing helpers.
    """

    graph_type = "all"
    graph_methods = frozenset().union(
        HeatmapVisualization.graph_methods,
        BoxplotVisualization.graph_methods,
        PixelVisualization.graph_methods,
        MapVisualization.graph_methods,
        ComparisonVisualization.graph_methods,
        TimingVisualization.graph_methods,
    )
            
