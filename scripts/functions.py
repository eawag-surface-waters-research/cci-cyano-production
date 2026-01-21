import os
import sys
import logging
from datetime import datetime

def parse_args(args):
    default_args = {
        "variable": "chla_mean", #
        "qa": "lwlr_quality_flag",
        "shapefile": "metadata.shp",
        "start_index": False,
        "end_index": False,
        "lakes": [],
        "out_folder": "",
        "threads": 1,
        "extract": True,
        "phenology": True,
        "qa_filter": True, # Only accept qa_flag = 0
        "spline_min_phase_length": 14,
        "spline_min_relative_amplitude": 0,
        "spline_min_phase_data": 0,
        "spline_data_gap_size": 31,
        "spline_date_gap_size_buffer": 0,
        "subs_peak_win_size": 365,
        "subs_peak_ampl_frac_list": [0.05, 0.35]
    }
    return default_args | args

def chunked(iterable, n):
    for i in range(0, len(iterable), n):
        yield iterable[i:i + n]

def set_logging(save):
    if save:
        filename = datetime.now().strftime("%Y%m%d_%H%M%S") + ".log"
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        logging.basicConfig(
            filename=os.path.join(repo, "logs", filename),
            filemode='a',
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO
        )
    else:
        logging.basicConfig(
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO
        )

def verify_arg_file(value):
    arg_folder = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "../args"))
    for file in os.listdir(arg_folder):
        if os.path.splitext(file)[0] == value or file == value:
            return os.path.join(arg_folder, file)
    raise ValueError("Argument file {} not found in the args folder.".format(value))