import argparse
import numpy as np
import logging
import os

from pathlib import Path
from pixel_metric_calcs.fit_metrics import FitMetric
from pixel_metric_calcs.bloom_assemble import BloomAssemble
from pixel_metric_calcs.bloom_prob import BloomProb
from pixel_metric_calcs.spatial_agg import SpatialAgg

def postprocess(lake, p, threads = 1, batch_size=100):

    # Path handling
    out_folder = Path(p["out_folder"])
    e_path = out_folder / "extract" / p["variable"] / f"{lake['id']}.nc"
    p_path = out_folder / "phenology" / p["variable"] / f"{lake['id']}.nc"
    if not os.path.isfile(e_path) or not os.path.isfile(p_path):
        logging.warning(f"Skipping lake {lake['id']}: extract or phenology file missing")
        return
    version_str = Path(p['out_folder']).name.removeprefix('v')
    logging.info(f"Postprocessing lake {lake['id']}")

    # Collect kwargs to pass to each subclass.
    calc_kwargs = {
        'lakeID' : lake['id'],
        'out_folder' : p['out_folder'],
        'variable' : p['variable'],
        'version' : version_str
    }

    # Postprocessing components
    if p['aggregation']:
        SpatialAgg(**calc_kwargs).run_chunked()
        logging.info(f"Completed aggregation for lake {lake['id']}")
 
    BloomAssemble(**calc_kwargs).run_chunked()
    logging.info(f"Completed bloom assembly for lake {lake['id']}")
    
    for start_year, end_year in p["time_splits"]:
        split_kwargs = {**calc_kwargs, "start_year": start_year, "end_year": end_year}

        FitMetric(**split_kwargs).run_chunked()
        logging.info(f"Completed fit metrics for lake {lake['id']} for years {start_year}-{end_year}")

        for qa_spec in p['kde_qa']:
            BloomProb(**split_kwargs,qa_value = qa_spec).run_chunked()
            logging.info(f"Completed bloom probability for lake {lake['id']} for years {start_year}-{end_year} and QA {qa_spec}")

    logging.info(f"Postprocessing lake {lake['id']} complete")