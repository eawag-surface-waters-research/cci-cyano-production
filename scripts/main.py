# -*- coding: utf-8 -*-
import json
import glob
import argparse
import itertools
import logging
import geopandas as gpd
from concurrent.futures import ProcessPoolExecutor
import os
import matplotlib as plt

from functions import set_logging, verify_arg_file, parse_args
from extract import extract
from phenology import phenology
from analysis import PhenologyEDA

def main(args, log=False, threads=1, parallel="lake", batch_size=100):
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
        if parallel == "pixels" or threads == 1:
            for lake in lakes:
                extract(lake, args, files, threads=threads)
        else:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                executor.map(extract, lakes, itertools.repeat(args), itertools.repeat(files))
    else:
        logging.info("Skipping extraction step.")

    if args["phenology"]:
        logging.info("Calculating pixelwise phenology metrics")
        if parallel == "pixels" or threads == 1:
            for lake in lakes:
                phenology(lake, args, threads=threads, batch_size=batch_size)
        else:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                executor.map(phenology, lakes, itertools.repeat(args),
                             itertools.repeat(1), itertools.repeat(batch_size))
                
    if args["analysis"]:
        logging.info("Starting Analysis")
        PhenologyEDA.set_shapefile_path(args["shapefile"])
        for lake in lakes:
            e_path = os.path.join(args["out_folder"], "extract", args["variable"], f"{lake['id']}.nc")
            p_path = os.path.join(args["out_folder"], "phenology", args["variable"], f"{lake['id']}.nc")
            if not os.path.isfile(e_path) or not os.path.isfile(p_path):
                logging.warning(f"Skipping lake {lake['id']}: extract or phenology file missing")
                continue
            logging.info(f"Analysing lake {lake['id']}")
            eda = PhenologyEDA(e_path, p_path)
            eda.r2_scores()
            eda.MAD_scores()
            eda.RMSE_scores()
            eda.correlation_scores()
            eda.values_per_pixel()
        logging.info(f"Analysis lake {lake['id']} complete")

        if args["plot"]:
            if not args["analysis"]:
                raise ValueError("Analysis must be true for plots to run.")
            logging.info("Starting Plots")
           




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse data from satellite images')
    parser.add_argument('--file', '-f', type=verify_arg_file, help='Name of the argument file in /args')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    parser.add_argument('--threads', '-t', help="Number of threads", default=1)
    parser.add_argument('--parallel', '-p', help="Run parallelisation on lakes or pixels", choices=['lakes', 'pixels'], default='lake')
    parser.add_argument('--batch-size', '-b', help="Number of pixels per batch for phenology I/O", type=int, default=100)
    args = parser.parse_args()
    with open(args.file) as f:
        file_args = json.load(f)
    main(file_args, log=args.logs, threads=int(args.threads), parallel=args.parallel, batch_size=args.batch_size)