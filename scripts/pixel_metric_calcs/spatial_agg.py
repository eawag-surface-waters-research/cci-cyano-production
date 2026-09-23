import os
from pathlib import Path
import netCDF4
import numpy as np
import multiprocessing
from functools import partial
import warnings
from contextlib import contextmanager
import functions as f
from pixel_metric_calcs.base import PixelCalcBase

class SpatialAgg(PixelCalcBase):
    def __init__(self, *args, **kwargs):

        kwargs = kwargs | {'metric_name': "background_median"}
        super().__init__(*args, **kwargs)
        
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
    
    def _get_spatial_inputs(self, nc):
        """Return time values, valid interior coordinates, and pixel locations."""
        lat = np.asarray(nc.variables["lat"][:])
        lon = np.asarray(nc.variables["lon"][:])

        t_all = f.unix_to_datenum(nc.variables["time"][:])
        nlat = len(nc.dimensions["lat"])
        nlon = len(nc.dimensions["lon"])

        coords = np.asarray(self.valid_coords, dtype=int).reshape(-1,2)
        interior = (
            (coords[:, 0] >= 1)
            & (coords[:, 0] < nlat - 1)
            & (coords[:, 1] >= 1)
            & (coords[:, 1] < nlon - 1)
        )
        coords = coords[interior]

        return (
            t_all,
            coords,
            lat[coords[:, 0]],
            lon[coords[:, 1]],
        )

    @staticmethod
    def _calculate_block(data_var, qa_var, coords, start, stop):
        """Calculate spatial aggregation values for time indices [start, stop)."""
        n_pixels = len(coords)
        values = np.full(
            (stop - start, n_pixels),
            np.nan,
            dtype=np.float32,
        )

        if n_pixels == 0:
            return values

        i_idx = coords[:, 0]
        j_idx = coords[:, 1]

        offsets_i = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1], dtype=int)
        offsets_j = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1], dtype=int)

        for output_index, time_index in enumerate(range(start, stop)):
            data_n = np.asarray(data_var[time_index], dtype=np.float32)
            qa_n = np.asarray(qa_var[time_index])

            data_windows = data_n[
                i_idx[:, None] + offsets_i,
                j_idx[:, None] + offsets_j,
            ]
            qa_windows = qa_n[
                i_idx[:, None] + offsets_i,
                j_idx[:, None] + offsets_j,
            ]

            invalid = (data_windows == -9999) | (qa_windows != 0)
            data_windows[invalid] = np.nan

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                values[output_index] = np.nanmedian(
                    data_windows,
                    axis=1,
                )

        return values

    def calculate(self):
        """Calculate and retain the complete result in memory."""
        with netCDF4.Dataset(self.data_path) as nc:
            (
                self.datenum_arr,
                self.coords,
                self.lat_vals,
                self.lon_vals,
            ) = self._get_spatial_inputs(nc)

            data_var = nc.variables[getattr(nc, "variable")]
            qa_var = nc.variables[getattr(nc, "qa")]
            ntime = len(nc.dimensions["time"])

            self.values = self._calculate_block(
                data_var,
                qa_var,
                self.coords,
                0,
                ntime,
            )

    def _create_output(self, ds, t_all, coords, lat_vals, lon_vals, block_size):
        """Create dimensions and variables for a spatial aggregation cache."""
        n_pixels = len(coords)
        ntime = len(t_all)

        ds.createDimension("time", ntime)
        ds.createDimension("pixel", n_pixels)

        ds.createVariable("time", "f8", ("time",))[:] = t_all
        ds.createVariable("pixel_i", "i4", ("pixel",))[:] = coords[:, 0]
        ds.createVariable("pixel_j", "i4", ("pixel",))[:] = coords[:, 1]
        ds.createVariable("lat", "f8", ("pixel",))[:] = lat_vals
        ds.createVariable("lon", "f8", ("pixel",))[:] = lon_vals

        chunksizes = None
        if n_pixels > 0 and ntime > 0:
            chunksizes = (
                min(block_size, ntime),
                min(256, n_pixels),
            )

        return ds.createVariable(
            "MA_value",
            "f4",
            ("time", "pixel"),
            fill_value=np.nan,
            zlib=True,
            complevel=4,
            chunksizes=chunksizes,
        )

    def write_output(self):
        """Write the complete in-memory result for a small lake."""
        with netCDF4.Dataset(self.save_fp, "w") as ds:
            ma_var = self._create_output(
                ds,
                self.datenum_arr,
                self.coords,
                self.lat_vals,
                self.lon_vals,
                block_size=len(self.datenum_arr),
            )
            ma_var[:, :] = self.values

    def calculate_and_write_chunked(self, block_size=64):
        """Calculate and write aggregation values without retaining all results."""
        with netCDF4.Dataset(self.data_path) as src:
            (
                t_all,
                coords,
                lat_vals,
                lon_vals,
            ) = self._get_spatial_inputs(src)

            data_var = src.variables[getattr(src, "variable")]
            qa_var = src.variables[getattr(src, "qa")]
            ntime = len(src.dimensions["time"])

            with netCDF4.Dataset(self.save_fp, "w") as dst:
                ma_var = self._create_output(
                    dst,
                    t_all,
                    coords,
                    lat_vals,
                    lon_vals,
                    block_size=block_size,
                )

                for start in range(0, ntime, block_size):
                    stop = min(start + block_size, ntime)

                    values = self._calculate_block(
                        data_var,
                        qa_var,
                        coords,
                        start,
                        stop,
                    )
                    ma_var[start:stop, :] = values


    @contextmanager
    def open_cached_metric(self):
        ds = netCDF4.Dataset(self.save_fp, "r")
        try:
            pixel_i = np.asarray(ds.variables["pixel_i"][:])
            pixel_j = np.asarray(ds.variables["pixel_j"][:])

            self._aggregation_pixel_index = {
                (i, j): idx
                for idx, (i, j) in enumerate(zip(pixel_i, pixel_j))
            }

            yield ds

        finally:
            ds.close()