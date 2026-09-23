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
        self.e_path = os.path.join(out_folder,"extract",variable,f"{lakeID}.nc")
        self.p_path = os.path.join(out_folder,"phenology",variable,f"{lakeID}.nc")
        self.version = version
        self.variable = variable
        self.lakeID = lakeID
        self.metric_name = metric_name
        self.valid_coords = f.valid_index_pairs(self.e_path)

    def build_path(self):
        pass

    def read_input(self):
        pass

    @abstractmethod
    def calculate(self):
        pass

    @abstractmethod
    def write_output(self):
        pass

    def calculate_and_write_chunked(self,block_size=64):
        pass

    def open_cached_metric(self):
        pass

    def run(self):
        self.build_path()
        self.read_input()
        self.calculate()
        self.write_output()

    def run_chunked(self, block_size=64):
        """Run the chunked calculation and writing workflow."""
        self.build_path()
        self.read_input()
        self.calculate_and_write_chunked(block_size=block_size)

