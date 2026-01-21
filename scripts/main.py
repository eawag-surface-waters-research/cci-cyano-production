# -*- coding: utf-8 -*-
import json
import glob
import argparse
import itertools
import logging
import geopandas as gpd
from concurrent.futures import ProcessPoolExecutor

from functions import set_logging, verify_arg_file, parse_args
from extract import extract
from phenology import phenology

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

    lakes = [row for i, (_, row) in enumerate(gdf.iterrows())]

    if args["extract"]:
        logging.info("Extracting {} from {}".format(args["variable"], args["images"]))
        if args["threads"] == 1:
            for lake in lakes:
                extract(lake, args, files)
        else:
            with ProcessPoolExecutor(max_workers=args["threads"]) as executor:
                executor.map(extract, lakes, itertools.repeat(args), itertools.repeat(files))
    else:
        logging.info("Skipping extraction step.")

    if args["phenology"]:
        logging.info("Calculating pixelwise phenology metrics")
        if args["threads"] == 1:
            for lake in lakes:
                phenology(lake, args)
        else:
            with ProcessPoolExecutor(max_workers=args["threads"]) as executor:
                executor.map(phenology, lakes, itertools.repeat(args))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse data from satellite images')
    parser.add_argument('--file', '-f', type=verify_arg_file, help='Name of the argument file in /args')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    args = parser.parse_args()
    with open(args.file) as f:
        file_args = json.load(f)
    main(file_args, args.logs)