import os
import sys
sys.path.append(os.path.abspath('../'))
sys.path.append(os.path.abspath('./scripts/'))
from pathlib import Path
import logging
import netCDF4
import re
import numpy as np
import pandas as pd
import datetime
import multiprocessing
from functools import partial
from scipy.stats import gaussian_kde
import warnings
from contextlib import contextmanager
import functions as f
from pixel_metric_calcs.base import PixelCalcBase
from pixel_metric_calcs.bloom_assemble import BloomAssemble

_GLOBALS = {}

OUTPUT_COLS = [
    "x_low",
    "x_high",
    "y_low",
    "y_high",
    "probability",
]

def _evaluate_probability_rows(kde, xi, y_rows):
    """Evaluate KDE density for a block of grid rows."""
    Xi, Yi = np.meshgrid(xi, y_rows)

    return kde(
        np.vstack([Xi.ravel(), Yi.ravel()])
    ).reshape(Xi.shape)


class BloomProb(PixelCalcBase):

    def __init__(
        self,
        *args,
        qa_value=None,
        start_year=0,
        end_year=9999,
        interval=21,
        resolution=1,
        x_max=400,
        y_max=730,
        **kwargs,
    ):

        kwargs = kwargs | {"metric_name": "bloom_probability"}
        super().__init__(*args, **kwargs)

        self.qa_value = qa_value
        self.start_year_filter = start_year
        self.end_year_filter = end_year

        self.interval = interval
        self.resolution = resolution
        self.x_max = x_max
        self.y_max = y_max


    def build_path(self):
        start_label = max(self.start_year, self.start_year_filter)
        end_label = min(self.end_year, self.end_year_filter)

        qa_label = (
            "all"
            if self.qa_value is None
            else "qa" + "".join(map(str, sorted(self.qa_value)))
        )

        out_dir = f.construct_bloom_prob_dirpath(
            self.out_folder,
            self.variable
        )
        os.makedirs(out_dir, exist_ok=True)

        filename = f.construct_bloom_prob_filename(
            lakeID = self.lakeID, start_year = start_label, end_year = end_label, qa_label = qa_label,
            interval = self.interval, resolution = self.resolution, x_max = self.x_max, y_max = self.y_max
        )

        self.save_fp = os.path.join(out_dir, filename)
        return out_dir, self.save_fp


    def read_input(self):
        """Locate bloom assembly cache file.
        Checks in 'kde_data' folder (legacy name for the data), 
        if not found in 'bloom_assembly_data' folder."""

        filename = f.construct_bloom_assembly_filename(self.lakeID)

        self.bloom_assembly_fp = f.construct_bloom_assembly_dirpath(self.out_folder, self.variable) /filename

        if not os.path.isfile(self.bloom_assembly_fp):

            self.bloom_assembly_fp = os.path.join(
                self.out_folder,
                "calculated_values",
                "kde_data",
                self.variable,
                filename,
            )

            if not os.path.isfile(self.bloom_assembly_fp):
                raise FileNotFoundError(
                    f"Could not find bloom assembly file: {filename}"
                )

        return self.bloom_assembly_fp


    def calculate(self):
        """Calculate and retain the complete result in memory."""

        fit = self._fit_bloom_kde()

        if fit is None:
            self.bloom_prob_df = pd.DataFrame()
            return self.bloom_prob_df

        self.kde, self.result_start_year, self.result_end_year, self.result_qa = fit
        self.bloom_prob_df = self._calculate_probability_grid(self.kde)

        return self.bloom_prob_df

    @staticmethod
    def _probability_dataframe(
        density,
        y_rows,
        xi,
        interval,
        resolution,
        first_window_row=0,
    ):
        """Convert a density block into probability-window rows."""
        n_cells = int(round(interval / resolution))

        if n_cells < 1:
            raise ValueError("interval/resolution must be at least one cell.")

        density = np.vstack(density) * resolution**2

        integral = np.cumsum(np.cumsum(density, axis=0),axis=1)
        integral = np.pad(integral, ((1, 0), (1, 0)))

        window_sum = (
            integral[n_cells:, n_cells:]
            - integral[:-n_cells, n_cells:]
            - integral[n_cells:, :-n_cells]
            + integral[:-n_cells, :-n_cells]
        )
        window_sum = np.clip(window_sum, 0, None)

        n_y, n_x = window_sum.shape
        y_low = np.asarray(y_rows)[
            first_window_row:first_window_row + n_y
        ]
        x_low, y_low = np.meshgrid(xi[:n_x], y_low)

        return pd.DataFrame({
            "x_low": x_low.ravel(),
            "x_high": (x_low + interval).ravel(),
            "y_low": y_low.ravel(),
            "y_high": (y_low + interval).ravel(),
            "probability": window_sum.ravel(),
        })


    def _calculate_probability_grid(self,kde):
        """Calculate the complete probability grid in memory."""
        xi = np.arange(-self.interval, self.x_max + self.interval, self.resolution)
        yi = np.arange(-self.interval, self.y_max + self.interval, self.resolution)

        density = _evaluate_probability_rows(kde, xi, yi)

        return self._probability_dataframe(
            density=density,
            y_rows=yi,
            xi=xi,
            interval=self.interval,
            resolution=self.resolution,
        )


    def _create_output(self, ds, n_events=None, block_size=64):
        """Create the BloomProbability output schema."""
        ds.createDimension("event", n_events)

        output_vars = {}

        for column in OUTPUT_COLS:
            output_vars[column] = ds.createVariable(
                column,
                "f8",
                ("event",),
                zlib=True,
                complevel=4,
                chunksizes=(
                    (min(block_size, n_events),)
                    if n_events not in (None, 0)
                    else None
                ),
            )

        ds.lake_id = str(self.lakeID)
        ds.version = str(self.version)
        ds.variable = str(self.variable)
        ds.metric_name = self.metric_name

        return output_vars


    def write_output(self):
        """Write the complete in-memory result."""
        if not hasattr(self, "bloom_prob_df"):
            raise ValueError("calculate() must be called before write_output().")

        with netCDF4.Dataset(self.save_fp, "w") as ds:
            output_vars = self._create_output(
                ds,
                n_events=len(self.bloom_prob_df),
            )

            for column, variable in output_vars.items():
                variable[:] = self.bloom_prob_df[column].to_numpy()

            ds.start_year = self.result_start_year
            ds.end_year = self.result_end_year
            ds.qa_value = (
                ""
                if self.result_qa is None
                else ",".join(map(str, sorted(self.result_qa)))
            )


    def calculate_and_write_chunked(self,block_size = 64):
        """Calculate and write probability windows in row blocks."""

        fit = self._fit_bloom_kde()

        if fit is None:
            return None

        self.kde, self.result_start_year, self.result_end_year, self.result_qa = fit

        xi = np.arange(-self.interval, self.x_max + self.interval, self.resolution)
        yi = np.arange(-self.interval, self.y_max + self.interval, self.resolution)
        n_cells = int(round(self.interval / self.resolution))

        with netCDF4.Dataset(self.save_fp, "w") as ds:
            output_vars = self._create_output(
                ds,
                n_events=None,
                block_size=block_size,
            )

            ds.start_year = self.result_start_year
            ds.end_year = self.result_end_year
            ds.qa_value = (
                ""
                if self.result_qa is None
                else ",".join(map(str, sorted(self.result_qa)))
            )

            event_start = 0

            for row_start in range(0, len(yi), block_size):
                row_stop = min(row_start + block_size, len(yi))

                # Include the following rows needed to complete the window.
                density_start = row_start
                density_stop = min(row_stop + n_cells - 1, len(yi))
                y_block = yi[density_start:density_stop]

                density = _evaluate_probability_rows(
                    self.kde,
                    xi,
                    y_block,
                )

                block_df = self._probability_dataframe(
                    density=density,
                    y_rows=y_block,
                    xi=xi,
                    interval=self.interval,
                    resolution=self.resolution,
                    first_window_row=0,
                )

                # Only write windows whose starting row belongs to this block.
                expected_rows = row_stop - row_start
                n_x = len(xi) - n_cells
                block_df = block_df.iloc[:expected_rows * n_x]

                event_stop = event_start + len(block_df)

                for column, variable in output_vars.items():
                    variable[event_start:event_stop] = (
                        block_df[column].to_numpy()
                    )

                event_start = event_stop

        return self.save_fp


    def read_full_cache(self, file_path):
        """Read a bloom-probability cache and return its data and filename metadata."""
        return f.read_bloom_probability_with_metadata(file_path)


    def _fit_bloom_kde(self):
        """
        Load/filter cached bloom events and fit a 2D gaussian_kde on
        (green-up advance DOY, green-down onset DOY).

        Returns
        -------
        tuple or None
            (kde, start, end, qa_filtered_set), or None if there isn't
            enough data to fit a KDE.
        """
        assembly_df = BloomAssemble.read_cached_metric(fp=self.bloom_assembly_fp)
        paired_df = assembly_df.set_index('year.DOY').groupby(['i','j']).apply(f.find_preceding_and_following, include_groups=False)
        paired_df = paired_df.droplevel([0,1]).dropna()
        qa_filtered_set = None
        if self.qa_value is not None:
            paired_df = paired_df[paired_df['qa_column'].isin(self.qa_value)]
            qa_filtered_set = self.qa_value
        print("Number of blooms considered:",paired_df.shape[0])

        if len(paired_df) < 2:
            warnings.warn("Not enough data to plot kde")
            return None

        paired_df["year"] = paired_df.index.astype(int)
        paired_df["doy_x"] = np.round((paired_df["prev_var1_time"] % 1) * 1000).astype(int)
        paired_df["doy_y"] = np.round((paired_df["next_var2_time"] % 1) * 1000).astype(int)

        x = paired_df["doy_x"].to_numpy()
        y = paired_df["doy_y"].to_numpy()
        y[y < x] += 365
        try:
            kde = gaussian_kde(np.vstack([x, y]))
        except np.linalg.LinAlgError:
            # too few / too degenerate (collinear or duplicate) points for a 2D KDE -
            # e.g. can happen with a restrictive qa_value filter that leaves very few events
            warnings.warn(f"Not enough distinct {self.variable} events to plot KDE for lake ID {self.lakeID}.")
            return None

        return kde, self.start_year_filter, self.end_year_filter, qa_filtered_set

if __name__ == "__main__":
    test_bloom = BloomProb(lakeID = 3500,
                               out_folder = "C:/Users/schelian/cci-cyano-production/data/v3.0",
                               variable = 'chla',
                               version = '3.0',
                               qa_value = None)
    test_bloom.build_path()
    print(test_bloom.read_input())
    test_bloom.calculate_and_write_chunked()

