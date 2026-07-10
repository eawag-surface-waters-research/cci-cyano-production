# -*- coding: utf-8 -*-
import json
import glob
import argparse
import itertools
import logging
import geopandas as gpd
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
from datetime import datetime, timezone

from functions import set_logging, verify_arg_file, parse_args, save_maps, save_pixel_plots, create_summary, save_comparison_plots, save_timing_plots, write_provenance, sanitize_filename
from extract import extract
from phenology import phenology
from visualization import PhenologyVisualization

def main(args, log=False, threads=1, parallel="lake", batch_size=100, args_file=None):
    set_logging(log)
    args = parse_args(args)
    run_id = datetime.now(timezone.utc).isoformat().replace(":", "-")

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
        if args["provenance"]:
            write_provenance(args["out_folder"], "extract", args, args_file=args_file,
                             extra={"lakes": [int(lake["id"]) for lake in lakes], "threads": threads, "parallel": parallel}, run_id=run_id)
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
        if args["provenance"]:
            write_provenance(args["out_folder"], "phenology", args, args_file=args_file,
                             extra={"lakes": [int(lake["id"]) for lake in lakes], "threads": threads,
                                    "parallel": parallel, "batch_size": batch_size}, run_id=run_id)

    if args["analysis"]:
        logging.info("Starting Analysis")
        PhenologyVisualization.set_shapefile_path(args["shapefile"])
        lake_analysis_folder = os.path.join(os.path.dirname(os.path.dirname(args["out_folder"])), "lake_analysis")
        if args["provenance"]:
            write_provenance(args["out_folder"], "analysis", args, args_file=args_file,
                             extra={"lakes": [int(lake["id"]) for lake in lakes]}, run_id=run_id)
        for lake in lakes:
            e_path = os.path.join(args["out_folder"], "extract", args["variable"], f"{lake['id']}.nc")
            p_path = os.path.join(args["out_folder"], "phenology", args["variable"], f"{lake['id']}.nc")
            if not os.path.isfile(e_path) or not os.path.isfile(p_path):
                logging.warning(f"Skipping lake {lake['id']}: extract or phenology file missing")
                continue
            logging.info(f"Analysing lake {lake['id']}")
            eda = PhenologyVisualization(e_path, p_path)
            lake_name = sanitize_filename(eda.ID_to_name(lake['id']).replace(" ", ""))
            lake_str = f"ID{lake['id']}_{lake_name}"
            eda.out_folder = Path(os.path.join(lake_analysis_folder, lake_str))
            eda.r2_scores(args["time_splits"])
            eda.MAD_scores(args["time_splits"])
            eda.RMSE_scores(args["time_splits"])
            eda.correlation_scores(args["time_splits"])
            eda.values_per_pixel(args["time_splits"])
            logging.info(f"Analysis lake {lake['id']} complete")

            if args["maps"]:
                logging.info(f"Starting maps for lake {lake['id']}")
                save_maps(eda, lake_analysis_folder, lake_str, time_splits = args["time_splits"])
                logging.info(f"Maps for {lake['id']} complete")

            if args["pixel_plots"]:
                lake_id_str = str(lake['id'])
                with open(args["pixel_plots"]) as f:
                    pixel_dict = json.load(f)
                if lake_id_str in pixel_dict:
                    logging.info(f"Starting pixel plots for lake {lake['id']}")
                    save_pixel_plots(eda, pixel_dict[lake_id_str], lake_analysis_folder, lake_str, time_splits= args["time_splits"], aggregation=args["aggregation"])
                    create_summary(eda, pixel_dict[lake_id_str], lake_analysis_folder, lake_str, time_splits= args["time_splits"])
                    logging.info(f"Pixel plots for lake {lake['id']} complete")
                else:
                    logging.info(f"WARNING: lake {lake['id']} is not in the pixel dictionary")

            if args["timing_plots"]:
                logging.info(f"Starting timing plots for lake {lake['id']}")
                save_timing_plots(eda, lake_analysis_folder, lake_str, time_splits=args["time_splits"], kde_qa=args["kde_qa"])
                logging.info(f"Timing plots for lake {lake['id']} complete")

            if args["comparison"]:
                logging.info("Starting Comparison Plots")
                if not args["pixel_plots"]:
                    logging.warning("Skipping comparison: pixel_plots arg is required")
                elif lake_id_str not in pixel_dict:
                    logging.warning(f"Skipping comparison for lake {lake['id']}: not in pixel dictionary")
                else:
                    _class_paths = {
                        "chla21":        ("v2.1", "chla_mean"),
                        "chla3":        ("v3.0", "chla"),
                        "phycocyanin3": ("v3.0", "phycocyanin"),
                    }
                    instances = {}
                    for class_name in args["comparison_classes"]:
                        base_ver, folder = _class_paths[class_name]
                        base_dir = args["out_folder"] if base_ver == os.path.basename(args["out_folder"]) \
                                   else os.path.join(os.path.dirname(args["out_folder"]), base_ver)
                        e = os.path.join(base_dir, "extract",   folder, f"{lake['id']}.nc")
                        p = os.path.join(base_dir, "phenology", folder, f"{lake['id']}.nc")
                        instances[class_name] = PhenologyVisualization(e, p) if (os.path.isfile(e) and os.path.isfile(p)) else None
                        if instances[class_name] is None:
                            logging.warning(f"Comparison: missing files for {class_name}, lake {lake['id']}")

                    if not all(instances.get(c) for c in args["comparison_classes"]):
                        logging.warning(f"Skipping comparison for lake {lake['id']}: one or more instances missing")
                    else:
                        save_comparison_plots(
                            instances, pixel_dict[lake_id_str], lake_analysis_folder, lake_str,
                            time_splits=args["time_splits"],
                            comparison_plot_types=args["comparison_plot_types"],
                            aggregation=args["aggregation"],
                            background_pts=args["background_pts"],
                            purple_chla21=args["purple_chla21"],
                            ratio_qa_source=args["ratio_qa_source"]
                        )
                        logging.info(f"Comparison plots for lake {lake['id']} complete")
                


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
    main(file_args, log=args.logs, threads=int(args.threads), parallel=args.parallel,
         batch_size=args.batch_size, args_file=os.path.splitext(os.path.basename(args.file))[0])