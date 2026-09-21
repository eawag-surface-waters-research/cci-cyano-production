import os
from pathlib import Path
import logging
import netCDF4
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
import warnings

import functions as f

class PixelCalcBase:
    def __init__(self, lakeID, out_folder, variable, version, metric_name, time_splits = None):
        self.out_folder = out_folder
        self.version = version
        self.variable = variable
        self.lakeID = lakeID
        self.metric_name = metric_name

    def build_path(self):
        pass

    def read_input(self):
        pass

    def calculate(self):
        pass

    def write_output(self):
        pass

    def read_cached_metric(self):
        pass

    def run(self):
        self.build_path()
        self.read_input()
        self.calculate()
        self.write_output()