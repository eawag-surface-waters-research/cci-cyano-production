import os
import sys
import logging
from datetime import datetime

def process_lake(lake, args):
    print(lake)

def parse_args(args):
    default_args = {
        "variable": "chla_mean",
        "shapefile": "metadata.shp",
        "start_index": False,
        "end_index": False,
        "min_date": "",
        "max_date": ""
    }
    return default_args | args

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