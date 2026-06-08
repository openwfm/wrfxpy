import sys

from fmda.fuel_moisture_da import execute_da_step, retrieve_mesowest_observations
from fmda.fuel_moisture_model import FuelMoistureModel
from fmda.var_wisdom import get_wisdom
from ingest.grib_file import GribFile
from ingest.HRRRA import HRRRA
from ingest.HRRR import HRRR
from utils import Dict, utc_to_esmf, esmf_to_utc
from vis.postprocessor import scalar_field_to_raster, vector_field_to_raster, scatter_to_raster
from fwi.fire_weather_indices import calculate_svp, calculate_eta
from ssh_shuttle import send_product_to_server

import pandas as pd
import netCDF4
import numpy as np
import json
import sys
import logging
import os
import os.path as osp
from datetime import datetime, timedelta, timezone
from typing import Sequence

# setup environment
sys_cfg = Dict(json.load(open("etc/conf.json")))
cfg = Dict(json.load(open("etc/fmda_cycler.json")))
meso_token = json.load(open("etc/tokens.json"))["mesowest"]

def get_hrrr_var_list(target_date):
    hrrr_vars = {
        datetime(2018, 7, 12): { # HRRRv3
            "rain": 634, "snow": 635, "t2": 616, "rh": 620, "psfc": 607, 
            "soil_t_0": 560, "soil_t_1": 562, "soil_t_2": 564, 
            "soil_t_3": 566, "soil_t_4": 568, "soil_t_5": 570, 
            "soil_t_6": 572, "soil_t_7": 574, "soil_t_8": 576, 
            "soil_moist_0": 561, "soil_moist_1": 563, "soil_moist_2": 565, 
            "soil_moist_3": 567, "soil_moist_4": 569, "soil_moist_5": 571, 
            "soil_moist_6": 573, "soil_moist_7": 575, "soil_moist_8": 577,
            "snowh": 615, "spfh": 618, "massden": None, "u10": 621, "v10": 622, 
            "ws": 623, "ust": 642, "znt": 641, "swdown": 661
        },
        datetime(2020, 12, 2): { # HRRRv4
            "rain": 635, "snow": 636, "t2": 616, "rh": 620, "psfc": 607, 
            "soil_t_0": 560, "soil_t_1": 562, "soil_t_2": 564, 
            "soil_t_3": 566, "soil_t_4": 568, "soil_t_5": 570, 
            "soil_t_6": 572, "soil_t_7": 574, "soil_t_8": 576, 
            "soil_moist_0": 561, "soil_moist_1": 563, "soil_moist_2": 565, 
            "soil_moist_3": 567, "soil_moist_4": 569, "soil_moist_5": 571, 
            "soil_moist_6": 573, "soil_moist_7": 575, "soil_moist_8": 577,
            "snowh": 615, "spfh": 618, "massden": 621, "u10": 622, "v10": 623, 
            "ws": 624, "ust": 643, "znt": 642, "swdown": 664
        }
    }
    # Sort the datetime keys
    sorted_dates = sorted(hrrr_vars.keys())

    for i, start_date in enumerate(sorted_dates):
        # If it's the last entry or the target date is before the next start
        if i == len(sorted_dates) - 1 or target_date < sorted_dates[i + 1]:
            if target_date >= start_date:
                return hrrr_vars[start_date]

    return {}  # If no matching config is found

    
## Import fmda_interval from hrrr_cycler
from hrrr_cycler import fmda_cycle_interval

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print(f"Invalid arguments. {len(sys.argv)} was given but 2 expected")
        print(('Usage: %s  <config_path>' % sys.argv[0]))
        print("Example: python src/hrrr_interval.py etc/hrrr_interval.yaml")
        sys.exit(-1)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    conf = Dict(json.load(open(sys.argv[1])))    
    from_utc = esmf_to_utc(conf.from_utc).replace(tzinfo=None)
    to_utc   = esmf_to_utc(conf.to_utc).replace(tzinfo=None)
    fmda_conf_path = conf.fmda_conf_path
    logging.info("Running `%s`,\n    from: %s to: %s", "fmda_cycle_interval", from_utc, to_utc)
    logging.info("FMDA Config File: %s", fmda_conf_path or "etc/fmda_cycler.json")

    fmda_cycle_interval(from_utc, to_utc, conf_path=fmda_conf_path)

