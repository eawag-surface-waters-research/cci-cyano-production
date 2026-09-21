import argparse
import numpy as np
import logging
import os
from pathlib import Path
from functions import unix_to_datetime, unix_to_datenum, datenum_to_datetime, remove_nan, sanitize_filename, ID_to_name
from visualization import PhenologyVisualization

def postprocess(lake, p, threads = 1, batch_size=100):
    PhenologyVisualization.set_shapefile_path(p["shapefile"])
    PhenologyVisualization.set_save_format(p["save_format"])
    lake_analysis_folder = os.path.join(os.path.dirname(os.path.dirname(p["out_folder"])), "lake_analysis")
    e_path = os.path.join(p["out_folder"], "extract", p["variable"], f"{lake['id']}.nc")
    p_path = os.path.join(p["out_folder"], "phenology", p["variable"], f"{lake['id']}.nc")
    if not os.path.isfile(e_path) or not os.path.isfile(p_path):
        logging.warning(f"Skipping lake {lake['id']}: extract or phenology file missing")
        return
    logging.info(f"Analysing lake {lake['id']}")
    eda = PhenologyVisualization(e_path, p_path)
    lake_name = sanitize_filename(ID_to_name(eda.gdf,lake['id']).replace(" ", ""))
    lake_str = f"ID{lake['id']}_{lake_name}"
    eda.out_folder = Path(os.path.join(lake_analysis_folder, lake_str))
    eda.r2_scores(p["time_splits"])
    eda.MAD_scores(p["time_splits"])
    eda.RMSE_scores(p["time_splits"])
    eda.correlation_scores(p["time_splits"])
    eda.values_per_pixel(p["time_splits"])
    eda.calculate_bloom_probabilities_from_kde()
    if p['aggregation']:
        eda.spatial_aggregation()
    logging.info(f"Analysis lake {lake['id']} complete")
    return eda.out_folder