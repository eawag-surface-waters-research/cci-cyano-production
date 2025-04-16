# -*- coding: utf-8 -*-
import os
import time
import json
import glob
import netCDF4
import argparse
import logging
import numpy as np
import xarray as xr
import geopandas as gpd
from rasterio import features, transform
from concurrent.futures import ProcessPoolExecutor

from functions import set_logging, verify_arg_file, parse_args


def process_lake(args):
    start = time.time()

    logging.info(f"Starting {args['lake']['id']}")
    output_path = os.path.join(args["out_folder"], "{}.nc".format(args["lake"]["id"]))
    output_temp = os.path.join(args["out_folder"], "{}_tmp.nc".format(args["lake"]["id"]))
    if os.path.isfile(output_path):
        logging.info(f"File {args['lake']['id']}.nc already exists")
        return
    if os.path.isfile(output_temp):
        os.remove(output_temp)

    with xr.open_dataset(args["files"][0]) as ds:
        lat = ds["lat"].values
        lon = ds["lon"].values

    geom = args["lake"].geometry
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
        data_var = out_nc.createVariable(args["variable"], "f4", ("time", "lat", "lon"), fill_value=-9999, zlib=True, complevel=4, shuffle=True)
        summary_var = out_nc.createVariable("summary", "f4", ("lat", "lon"))

        lat_var[:] = lat_sub
        lon_var[:] = lon_sub

        index = 0
        start_time = time.time()
        for i, file in enumerate(args["files"]):
            with netCDF4.Dataset(file) as nc:
                t = int(nc.variables["time"][0])
                data_slice = np.array(nc.variables[args["variable"]][0, lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1])
                data_slice[data_slice > 10e6] = np.nan
                data_slice[~mask] = np.nan
                nan_mask = np.isnan(data_slice)
                if not np.all(nan_mask):
                    summary = summary + (~nan_mask).astype(int)
                    data_slice[nan_mask] = -9999
                    data_var[index, :, :] = data_slice
                    time_var[index] = t
                    index = index + 1
            if (i + 1) % 1000 == 0:
                elapsed = time.time() - start_time
                avg_time = (elapsed / (i + 1) * 1000) / args["threads"]
                logging.info(f"Average time per file: {avg_time:.2f} ms")

        summary_var[:] = summary

    os.rename(output_temp, output_path)
    elapsed = time.time() - start
    logging.info(f"Completed {args['lake']['id']} in {elapsed:.2f} seconds.")


def main(args, log=False):
    set_logging(log)
    args = parse_args(args)

    logging.info("Reading lake shapefile from: {}".format(args["shapefile"]))
    gdf = gpd.read_file(args["shapefile"])

    if args["start_index"] and args["end_index"]:
        gdf = gdf.iloc[args["start_index"]:args["end_index"]]
    elif args["lakes"]:
        gdf = gdf[gdf['id'].isin(args["lakes"])]

    files = glob.glob(f"{args['images']}/**/*.nc", recursive=True)
    files.sort()

    lakes = [{"lake": row, "variable": args['variable'], "files": files, "out_folder": args['out_folder'],
              "threads": args['threads']} for i, (_, row)
             in enumerate(gdf.iterrows())]

    if args["threads"] == 1:
        for lake in lakes:
            process_lake(lake)
    else:
        with ProcessPoolExecutor(max_workers=args["threads"]) as executor:
            executor.map(process_lake, lakes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse data from satellite images')
    parser.add_argument('--file', '-f', type=verify_arg_file, help='Name of the argument file in /args')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    args = parser.parse_args()
    with open(args.file) as f:
        file_args = json.load(f)
    main(file_args, args.logs)