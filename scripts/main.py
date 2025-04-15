# -*- coding: utf-8 -*-
import os
import glob
import json
import argparse
import logging
import xarray as xr
import geopandas as gpd
import numpy as np
from rasterio import features, transform
from functions import set_logging, verify_arg_file, parse_args, chunked


def main(args, log=False):
    set_logging(log)
    args = parse_args(args)

    logging.info("Reading lake shapefile from: {}".format(args["shapefile"]))
    gdf = gpd.read_file(args["shapefile"])
    gdf = gdf[gdf['id'] == 18]

    if args["start_index"] and args["end_index"]:
        gdf = gdf.iloc[args["start_index"]:args["end_index"]]

    logging.info("Creating pixel mask")
    files = glob.glob(f"{args['images']}/**/*.nc", recursive=True)
    files = files[-1000:]
    files.sort()

    ds = xr.open_dataset(files[0])
    lat = ds["lat"].values
    lon = ds["lon"].values

    mask_dict = {}

    for idx, lake in gdf.iterrows():
        try:
            geom = lake.geometry
            minx, miny, maxx, maxy = geom.bounds
            lat_mask = (lat >= miny) & (lat <= maxy)
            lon_mask = (lon >= minx) & (lon <= maxx)

            lat_indices = np.where(lat_mask)[0]  # Indices for latitudes
            lon_indices = np.where(lon_mask)[0]

            print(f"Lat {lat_indices[0]},{lat_indices[-1]}")
            print(f"Lon {lon_indices[0]},{lon_indices[-1]}")


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
                mask = features.rasterize(
                    [(geom, 1)],
                    out_shape=lon_sub_grid.shape,
                    transform=t,
                    fill=0,
                    all_touched=True,
                    dtype="uint8"
                ).astype(bool)

            mask_dict[lake["id"]] = {"lat": lat[lat_mask], "lon": lon[lon_mask], "mask": np.flipud(mask)}
        except:
            raise
    ds.close()

    logging.info("Processing data")
    ds = xr.open_mfdataset(files, combine='nested', concat_dim='time')
    logging.info(f"Loaded chunk")
    for id in mask_dict.keys():
        logging.info("Extracting {}".format(id))
        data = ds["chla_mean"].sel(lat=mask_dict[id]["lat"], lon=mask_dict[id]["lon"]).where(mask_dict[id]["mask"])
        valid_time_steps = data.notnull().any(dim=['lat', 'lon'])
        if valid_time_steps.any():
            out = os.path.join(args["out"], "{}/2020.nc".format(id))
            os.makedirs(os.path.dirname(out), exist_ok=True)
            data.sel(time=valid_time_steps).to_netcdf(out, mode='w', unlimited_dims=["time"])
        else:
            logging.info("No valid")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse data from satellite images')
    parser.add_argument('--file', '-f', type=verify_arg_file, help='Name of the argument file in /args')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    args = parser.parse_args()
    with open(args.file) as f:
        file_args = json.load(f)
    main(file_args, args.logs)