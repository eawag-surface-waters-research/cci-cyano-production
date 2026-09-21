import os
from pathlib import Path
import logging
import netCDF4
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
import warnings

from base import MetricBase
import functions as f

class SpatialAgg(MetricBase):
    def build_path(self):
        out_dir = os.path.join(
                self.out_folder,
                "calculated_values",
                "spatial_aggregation_values",
                self.variable,
            )
        os.makedirs(out_dir, exist_ok=True)
        self.save_fp = os.path.join(out_dir, f"ID{self.lakeID}_background_agg.nc")

    def read_input(self):
        self.data_path = os.path.join(self.out_folder,'extract',self.variable,f"{self.lakeID}.nc")
    
    def calculate(self):
        with netCDF4.Dataset(self.data_path) as nc:
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
                # self._write_aggregation_netcdf(t_all,
                #     np.empty(0, dtype=int), np.empty(0, dtype=int),
                #     np.empty(0), np.empty(0),
                #     np.empty((ntime, 0), dtype=np.float32),
                # )
                return

            i_idx = coords[:, 0]
            j_idx = coords[:, 1]

            # indices for median_grid, which is smaller by 1 border cell each side
            ii = i_idx - 1
            jj = j_idx - 1

            lat_vals = lat[i_idx]
            lon_vals = lon[j_idx]

            n_pixels = len(coords)
            # netcdf path accumulates a (ntime, n_pixels) array and writes it in one
            # bulk call: writing timestep-by-timestep into chunks that span the full
            # time dimension would force a decompress/recompress of every chunk on
            # every iteration.
            values = np.full((ntime, n_pixels), np.nan, dtype=np.float32)

            for n in range(ntime):
                data_n = np.asarray(data_var[n], dtype=np.float32)
                qa_n = np.asarray(qa_var[n])

                # shape: (nlat-2, nlon-2, 3, 3)
                data_windows = np.lib.stride_tricks.sliding_window_view(data_n, (3, 3))
                qa_windows = np.lib.stride_tricks.sliding_window_view(qa_n, (3, 3))

                valid_mask = (data_windows != -9999) & (qa_windows == 0)

                masked = data_windows.astype(np.float32, copy=True)
                masked[~valid_mask] = np.nan

                # shape: (nlat-2, nlon-2)
                median_grid = np.nanmedian(masked, axis=(-2, -1))

                ma_values = median_grid[ii, jj]

                values[n, :] = ma_values

    def write_output(self, t_all, i_idx, j_idx, lat_vals, lon_vals, values):
        """Write spatial_aggregation() results to a compressed (time, pixel) NetCDF cache.

        lat/lon are stored once per pixel rather than once per row (the CSV's main
        source of bloat), and MA_value is chunked as (ntime, 1) so a per-pixel read
        - the only access pattern plot_background_pts uses - pulls exactly one
        contiguous chunk instead of scanning the whole file.
        """
        n_pixels = len(i_idx)
        with netCDF4.Dataset(self.save_fp, "w") as ds:
            ds.createDimension("time", len(t_all))
            ds.createDimension("pixel", n_pixels)

            ds.createVariable("time", "f8", ("time",))[:] = t_all
            ds.createVariable("pixel_i", "i4", ("pixel",))[:] = i_idx
            ds.createVariable("pixel_j", "i4", ("pixel",))[:] = j_idx
            ds.createVariable("lat", "f8", ("pixel",))[:] = lat_vals
            ds.createVariable("lon", "f8", ("pixel",))[:] = lon_vals

            chunksizes = (len(t_all), 1) if n_pixels > 0 else None
            ma_var = ds.createVariable(
                "MA_value", "f4", ("time", "pixel"),
                fill_value=np.nan, zlib=True, complevel=4,
                chunksizes=chunksizes,
            )
            ma_var[:, :] = values

    def read_cached_metric(self):
        pass