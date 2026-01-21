# -*- coding: utf-8 -*-
import os
import time
import netCDF4
import logging
import numpy as np
import xarray as xr
from rasterio import features, transform

import matplotlib.pyplot as plt

def extract(lake, parameters, files):
    try:
        start = time.time()
        logging.info(f"Starting extraction for {lake['id']}")
        output_path = os.path.join(parameters["out_folder"], "extract", parameters["variable"], "{}.nc".format(lake["id"]))
        output_temp = os.path.join(parameters["out_folder"], "extract", parameters["variable"], "{}_tmp.nc".format(lake["id"]))
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        if os.path.isfile(output_path):
            logging.info(f"File {lake['id']}.nc already exists")
            return
        if os.path.isfile(output_temp):
            os.remove(output_temp)

        with xr.open_dataset(files[0]) as ds:
            lat = ds["lat"].values
            lon = ds["lon"].values

        geom = lake.geometry
        minx, miny, maxx, maxy = geom.bounds
        lat_mask = (lat >= miny) & (lat <= maxy)
        lon_mask = (lon >= minx) & (lon <= maxx)

        lat_indices = np.where(lat_mask)[0]
        lon_indices = np.where(lon_mask)[0]

        lat_sub = lat[lat_mask]
        lon_sub = lon[lon_mask]
        lon_sub_grid, lat_sub_grid = np.meshgrid(lon_sub, lat_sub)

        if lat_sub.shape[0] == 1 and lon_sub.shape[0] == 1:
            mask = np.array([[True]])
        else:
            t = transform.from_bounds(
                west=lon_sub_grid.min(), east=lon_sub_grid.max(),
                south=lat_sub_grid.min(), north=lat_sub_grid.max(),
                width=lon_sub_grid.shape[1], height=lat_sub_grid.shape[0]
            )
            mask = np.flipud(features.rasterize(
                [(geom, 1)],
                out_shape=lon_sub_grid.shape,
                transform=t,
                fill=0,
                all_touched=True,
                dtype="uint8"
            ).astype(bool))

        summary = np.zeros_like(mask, dtype=int)

        lat_dim = lat_sub.shape[0]
        lon_dim = lon_sub.shape[0]

        with netCDF4.Dataset(output_temp, "w", format="NETCDF4") as out_nc:
            out_nc.createDimension("time", None)
            out_nc.createDimension("lat", lat_dim)
            out_nc.createDimension("lon", lon_dim)

            time_var = out_nc.createVariable("time", "f4", ("time",))
            lat_var = out_nc.createVariable("lat", "f4", ("lat",))
            lon_var = out_nc.createVariable("lon", "f4", ("lon",))
            data_var = out_nc.createVariable(parameters["variable"], "f4", ("time", "lat", "lon"), fill_value=-9999, zlib=True, complevel=4, shuffle=True)
            qa_var = out_nc.createVariable(parameters["qa"], "i1", ("time", "lat", "lon"), fill_value=0, zlib=True, complevel=4, shuffle=True)
            qa_var.long_name = "Quality Assessment Flag"
            qa_var.description = "Data quality indicator: 0 = Good, 1 = Fair, 2 = Poor, 3 = No data"
            summary_var = out_nc.createVariable("summary", "f4", ("lat", "lon"))
            qa_var.description = "Number of valid values per pixel over the full timeseries"

            lat_var[:] = lat_sub
            lon_var[:] = lon_sub

            index = 0
            start_time = time.time()
            for i, file in enumerate(files):
                with netCDF4.Dataset(file) as nc:
                    data_slice = np.array(nc.variables[parameters["variable"]][0, lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1])
                    data_slice[data_slice > 10e6] = np.nan
                    data_slice[~mask] = np.nan
                    nan_mask = np.isnan(data_slice)
                    if not np.all(nan_mask):
                        summary = summary + (~nan_mask).astype(int)
                        data_slice[nan_mask] = -9999
                        qa_slice = np.array(nc.variables[parameters["qa"]][0, lat_indices[0]:lat_indices[-1] + 1, lon_indices[0]:lon_indices[-1] + 1])
                        qa_slice[qa_slice > 10e6] = 3
                        qa_slice[nan_mask] = 3
                        data_var[index, :, :] = data_slice
                        qa_var[index, :, :] = qa_slice
                        time_var[index] = int(nc.variables["time"][0])
                        index = index + 1
                if (i + 1) % 1000 == 0:
                    elapsed = time.time() - start_time
                    avg_time = (elapsed / (i + 1) * 1000) / parameters["threads"]
                    logging.info(f"Average time per file: {avg_time:.2f} ms")
            summary_var[:] = summary

        os.rename(output_temp, output_path)
        elapsed = time.time() - start
        logging.info(f"Completed extraction for {lake['id']} in {elapsed:.2f} seconds.")
        return output_path
    except Exception as e:
        logging.error(e)
        raise