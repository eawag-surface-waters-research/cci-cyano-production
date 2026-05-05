# -*- coding: utf-8 -*-
import json
import glob
import argparse
import itertools
import logging
import geopandas as gpd
from concurrent.futures import ProcessPoolExecutor
import os
import matplotlib.pyplot as plt

from functions import set_logging, verify_arg_file, parse_args, save_maps, save_pixel_plots, create_summary
from extract import extract
from phenology import phenology
from visualization import PhenologyVisualization

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
        PhenologyVisualization.set_shapefile_path(args["shapefile"])
        lake_analysis_folder = os.path.join(os.path.dirname(os.path.dirname(args["out_folder"])), "lake_analysis")
        for lake in lakes:
            e_path = os.path.join(args["out_folder"], "extract", args["variable"], f"{lake['id']}.nc")
            p_path = os.path.join(args["out_folder"], "phenology", args["variable"], f"{lake['id']}.nc")
            if not os.path.isfile(e_path) or not os.path.isfile(p_path):
                logging.warning(f"Skipping lake {lake['id']}: extract or phenology file missing")
                continue
            logging.info(f"Analysing lake {lake['id']}")
            eda = PhenologyVisualization(e_path, p_path)
            lake_name = eda.ID_to_name(lake['id']).replace(" ", "")
            lake_str = f"ID{lake['id']}_{lake_name}"
            eda.out_folder = os.path.join(lake_analysis_folder, lake_str)
            eda.r2_scores()
            eda.MAD_scores()
            eda.RMSE_scores()
            eda.correlation_scores()
            eda.values_per_pixel()
            logging.info(f"Analysis lake {lake['id']} complete")

            if args["maps"]:
                logging.info(f"Starting maps for lake {lake['id']}")
                save_maps(eda, lake_analysis_folder, lake_str, start = args["start"], end=args["end"])
                logging.info(f"Maps for {lake['id']} complete")

            if args["pixel_plots"]:
                lake_id_str = str(lake['id'])
                with open(args["pixel_plots"]) as f:
                    pixel_dict = json.load(f)
                if lake_id_str in pixel_dict:
                    logging.info(f"Starting pixel plots for lake {lake['id']}")
                    save_pixel_plots(eda, pixel_dict[lake_id_str], lake_analysis_folder, lake_str, end=args["end"], start=args["start"], split_start=args["split_start"], split_end=args["split_end"], aggregation=args["aggregation"])
                    create_summary(eda, pixel_dict[lake_id_str], lake_analysis_folder, lake_str)
                    logging.info(f"Pixel plots for lake {lake['id']} complete")
                else:
                    logging.info(f"WARNING: lake {lake['id']} is not in the pixel dictionary")
                    continue
                    
            if args["comparison"]:
                logging.info("Starting Comparison Plots")
                if not args["pixel_plots"]:
                    logging.warning("Skipping comparison: pixel_plots arg is required")
                else:
                    PhenologyVisualization.set_shapefile_path(args["shapefile"])
                    with open(args["pixel_plots"]) as f:
                        pixel_dict = json.load(f)

                    _class_paths = {
                        "chla21":        ("v2.1", "chla_mean"),
                        "chla31":        ("v3.1", "chla"),
                        "phycocyanin31": ("v3.1", "phycocyanin"),
                    }

                    for lake in lakes:
                        lake_id = lake['id']
                        lake_id_str = str(lake_id)
                        if lake_id_str not in pixel_dict:
                            continue

                        # Build the three PhenologyVisualization instances
                        instances = {}
                        for class_name in args["comparison_classes"]:
                            base_ver, folder = _class_paths[class_name]
                            if base_ver == os.path.basename(args["out_folder"]):
                                base_dir = args["out_folder"]
                            else:
                                base_dir = os.path.join(os.path.dirname(args["out_folder"]), base_ver)
                            e_path = os.path.join(base_dir, "extract", folder, f"{lake_id}.nc")
                            p_path = os.path.join(base_dir, "phenology", folder, f"{lake_id}.nc")
                            if not os.path.isfile(e_path) or not os.path.isfile(p_path):
                                logging.warning(f"Comparison: missing files for {class_name}, lake {lake_id}")
                                instances[class_name] = None
                            else:
                                instances[class_name] = PhenologyVisualization(e_path, p_path)

                        chla21 = instances.get("chla21")
                        chla31 = instances.get("chla31")
                        phyco  = instances.get("phycocyanin31")
                        if not all([chla21, chla31, phyco]):
                            logging.warning(f"Skipping comparison for lake {lake_id}: one or more instances missing")
                            continue

                        name     = chla21.ID_to_name(lake_id)
                        name_str = name.replace(" ", "")
                        lake_str = f"ID{lake_id}_{name_str}"
                        agg      = "_agg" if args["aggregation"] else ""
                        pt       = args["comparison_plot_types"]
                        bg       = args["background_pts"]
                        purple   = args["purple_chla21"]

                        for index in pixel_dict[lake_id_str]:
                            i, j      = index
                            save_path = os.path.join(lake_analysis_folder, lake_str,
                                                    "plots", "pixel_plots", "comparisons", f"{i}_{j}")
                            os.makedirs(save_path, exist_ok=True)

                            if "chla21 vs chla31" in pt:
                                fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                                chla21.extrema_comparison(chla31, i, j, ax, aggregation=args["aggregation"])
                                fig.savefig(os.path.join(save_path, f"extrema_{chla21.variable}_v{chla21.version.replace('.','')}_VS_{chla31.variable}_v{chla31.version.replace('.','')}_fullts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "chla21 vs chla31_split" in pt:
                                fig, axs = plt.subplots(1, 2, figsize=(20, 8))
                                chla21.extrema_comparison(chla31, i, j, axs[0], aggregation=args["aggregation"], end=args["split_end"])
                                chla21.extrema_comparison(chla31, i, j, axs[1], aggregation=args["aggregation"], start=args["split_start"])
                                fig.savefig(os.path.join(save_path, f"extrema_{chla21.variable}_v{chla21.version.replace('.','')}_VS_{chla31.variable}_v{chla31.version.replace('.','')}_split_ts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "chla21 vs phyco" in pt:
                                fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                                chla21.extrema_comparison(phyco, i, j, ax, aggregation=args["aggregation"])
                                fig.savefig(os.path.join(save_path, f"extrema_{chla21.variable}_v{chla21.version.replace('.','')}_VS_{phyco.variable}_v{phyco.version.replace('.','')}_fullts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "chla21 vs phyco_split" in pt:
                                fig, axs = plt.subplots(1, 2, figsize=(20, 8))
                                chla21.extrema_comparison(phyco, i, j, axs[0], aggregation=args["aggregation"], end=args["split_end"])
                                chla21.extrema_comparison(phyco, i, j, axs[1], aggregation=args["aggregation"], start=args["split_start"])
                                fig.savefig(os.path.join(save_path, f"extrema_{chla21.variable}_v{chla21.version.replace('.','')}_VS_{phyco.variable}_v{phyco.version.replace('.','')}_split_ts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "chla31 vs phyco" in pt:
                                fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                                chla31.extrema_comparison(phyco, i, j, ax, aggregation=args["aggregation"])
                                fig.savefig(os.path.join(save_path, f"extrema_{chla31.variable}_v{chla31.version.replace('.','')}_VS_{phyco.variable}_v{phyco.version.replace('.','')}_fullts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "chla31 vs phyco_split" in pt:
                                fig, axs = plt.subplots(1, 2, figsize=(20, 8))
                                chla31.extrema_comparison(phyco, i, j, axs[0], aggregation=args["aggregation"], end=args["split_end"])
                                chla31.extrema_comparison(phyco, i, j, axs[1], aggregation=args["aggregation"], start=args["split_start"])
                                fig.savefig(os.path.join(save_path, f"extrema_{chla31.variable}_v{chla31.version.replace('.','')}_VS_{phyco.variable}_v{phyco.version.replace('.','')}_split_ts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "triple" in pt:
                                fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                                chla21.extrema_comparison(chla31, i, j, ax, aggregation=args["aggregation"],
                                                        background_pts=bg, other2=phyco, purple_chla21=purple)
                                fig.savefig(os.path.join(save_path, f"extrema_{chla21.variable}_v{chla21.version.replace('.','')}_VS_{chla31.variable}_v{chla31.version.replace('.','')}_VS_{phyco.variable}_v{phyco.version.replace('.','')}_fullts{agg}.png"), dpi=600)
                                plt.close(fig)

                            if "triple_split" in pt:
                                fig, axs = plt.subplots(1, 2, figsize=(20, 8))
                                chla21.extrema_comparison(chla31, i, j, axs[0], aggregation=args["aggregation"],
                                                        end=args["split_end"], background_pts=bg, other2=phyco, purple_chla21=purple)
                                chla21.extrema_comparison(chla31, i, j, axs[1], aggregation=args["aggregation"],
                                                        start=args["split_start"], background_pts=bg, other2=phyco, purple_chla21=purple)
                                suffix = "_purple_chla21" if purple else ""
                                fig.savefig(os.path.join(save_path, f"extrema_{chla21.variable}_v{chla21.version.replace('.','')}_VS_{chla31.variable}_v{chla31.version.replace('.','')}_VS_{phyco.variable}_v{phyco.version.replace('.','')}_split_ts{agg}{suffix}.png"), dpi=600)
                                plt.close(fig)

                            logging.info(f"Comparison plots for lake {lake_id} pixel ({i},{j}) complete")
                


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