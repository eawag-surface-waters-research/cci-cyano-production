# -*- coding: utf-8 -*-
import os
import sys
import yaml
import glob
import json
import time
import argparse
import logging
import pandas as pd
import xarray as xr
import geopandas as gpd
import numpy as np
from rasterio import features, transform
import matplotlib.pyplot as plt
from functions import set_logging, verify_arg_file, parse_args, process_lake


def main(args, log=False):
    set_logging(log)
    args = parse_args(args)

    logging.info("Reading lake shapefile from: {}".format(args["shapefile"]))
    gdf = gpd.read_file(args["shapefile"])

    logging.info("Creating pixel mask")
    files = glob.glob(f"{args['images']}/**/*.nc", recursive=True)
    files.sort()

    files = ["/media/runnalja/My Passport/Documents/Phenology/data/CCI_lakes_nc/v2.0.2/2020/06/ESACCI-LAKES-L3S-LK_PRODUCTS-MERGED-20200615-fv2.0.2.nc"]

    ds = xr.open_dataset(files[0])
    lat = ds["lat"].values
    lon = ds["lon"].values

    for idx, lake in gdf.iterrows():
        try:
            geom = lake.geometry
            minx, miny, maxx, maxy = geom.bounds
            lat_mask = (lat >= miny) & (lat <= maxy)
            lon_mask = (lon >= minx) & (lon <= maxx)

            data_sub = ds["chla_mean"].sel(lat=lat[lat_mask], lon=lon[lon_mask])
            lon_sub, lat_sub = np.meshgrid(data_sub["lon"], data_sub["lat"])

            t = transform.from_bounds(
                west=lon_sub.min(), east=lon_sub.max(),
                south=lat_sub.min(), north=lat_sub.max(),
                width=lon_sub.shape[1], height=lat_sub.shape[0]
            )

            mask = features.rasterize(
                [(geom, 1)],
                out_shape=lon_sub.shape,
                transform=t,
                fill=0,
                all_touched=True,
                dtype="uint8"
            ).astype(bool)

            # 3. Mask the data
            lake_data = data_sub.where(np.flipud(mask))
        except:
            print(lake["id"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Simstrat on an operational basis')
    parser.add_argument('--file', '-f', type=verify_arg_file, help='Name of the argument file in /args')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    args = parser.parse_args()
    with open(args.file) as f:
        file_args = json.load(f)
    main(file_args, args.logs)
