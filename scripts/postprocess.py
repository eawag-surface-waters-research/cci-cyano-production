import argparse
import numpy as np
import logging
import os

from pathlib import Path
from functions import unix_to_datetime, unix_to_datenum, datenum_to_datetime, remove_nan, sanitize_filename, ID_to_name
from visualization import PhenologyVisualization
from pixel_metric_calcs.fit_metrics import FitMetric
from pixel_metric_calcs.bloom_assemble import BloomAssemble
from pixel_metric_calcs.bloom_prob import BloomProb
from pixel_metric_calcs.spatial_agg import SpatialAgg

def postprocess(lake, p, threads = 1, batch_size=100):
    e_path = os.path.join(p["out_folder"], "extract", p["variable"], f"{lake['id']}.nc")
    p_path = os.path.join(p["out_folder"], "phenology", p["variable"], f"{lake['id']}.nc")
    if not os.path.isfile(e_path) or not os.path.isfile(p_path):
        logging.warning(f"Skipping lake {lake['id']}: extract or phenology file missing")
        return
    version_str = Path(p['out_folder']).name.removeprefix('v')
    logging.info(f"Postprocessing lake {lake['id']}")
    calc_kwargs = {
        'lakeID' : lake['id'],
        'out_folder' : p['out_folder'],
        'variable' : p['variable'],
        'version' : version_str
    }
    if p['aggregation']:
        SpatialAgg(calc_kwargs).run_chunked()
    for start_year, end_year in p["time_split"]:
        calc_kwargs['start_year'] = start_year
        calc_kwargs['end_year'] = end_year
        FitMetric(calc_kwargs).run_chunked()
        BloomAssemble(calc_kwargs).run_chunked()
        BloomProb(calc_kwargs).run_chunked()

    logging.info(f"Postprocessing lake {lake['id']} complete")