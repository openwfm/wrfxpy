# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from fmda.fuel_moisture_da import execute_da_step, retrieve_mesowest_observations
from fmda.fuel_moisture_model import FuelMoistureModel
from fmda.var_wisdom import get_wisdom
from ingest.grib_file import GribFile
from ingest.HRRRA import HRRRA
from ingest.HRRR import HRRR
from utils import Dict, ensure_dir, utc_to_esmf, delete, force_copy, read_yml
from vis.postprocessor import scalar_field_to_raster, vector_field_to_raster, scatter_to_raster
from fwi.fire_weather_indices import calculate_svp, calculate_eta
from ssh_shuttle import send_product_to_server
from fmda.moisture_rnn_operational import OperationalRNNPredictor

import pandas as pd
import netCDF4
from netCDF4 import num2date
import numpy as np
import json
import sys 
import logging
import os
import os.path as osp
from pathlib import Path
import joblib
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

def write_postprocess(
        mf, postproc_path, cycle_dir, esmf_cycle, name, 
        raster_png, coords, cb_png=None, levels=None, alpha=None
    ):
    """
    Write postprocessing files.

    :param post: the UTC cycle time
    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :return: the postprocessing path
    """
    raster_name = f"{cycle_dir}-{name}-raster.png"
    cb_name = f"{cycle_dir}-{name}-raster-cb.png"
    with open(osp.join(postproc_path, raster_name), "wb") as f:
        f.write(raster_png)
    mf["1"][esmf_cycle][name] = { "raster" : raster_name, "coords" : coords }
    if cb_png is not None:
        with open(osp.join(postproc_path, cb_name), "wb") as f:
            f.write(cb_png) 
        mf["1"][esmf_cycle][name].update({ "colorbar": cb_name })
    if levels is not None:
        mf["1"][esmf_cycle][name].update({ "levels" : levels })
    if alpha is not None:
        mf["1"][esmf_cycle][name].update({ "alpha" : alpha })


def postprocess_cycle(cycle, region_cfg, wksp_path, fcst_hour, bounds=None):
    """
    Build rasters from the computed fuel moisture.

    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :param bounds: bounding box of the post-processing
    :return: the postprocessing path
    """
    model_path = compute_model_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    cycle_id = compute_fmda_id(cycle, region_cfg.code)
    postproc_path = compute_postproc_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    if fcst_hour > 0:
        cycle_id += f"-f{fcst_hour:02d}"
        prev_cycle = cycle
        prev_fcst_hour = fcst_hour - 1
        prev_postproc_path = compute_postproc_path(prev_cycle, region_cfg.code, wksp_path, prev_fcst_hour)
    else:
        # get the most recent forecast product or at the end previous real-time
        prev_cycle = cycle - timedelta(hours=1)
        prev_fcst_hour = region_cfg.forecast_length
        prev_postproc_path = compute_postproc_path(prev_cycle, region_cfg.code, wksp_path, prev_fcst_hour)
        while not osp.exists(prev_postproc_path) and prev_fcst_hour > 0:
            prev_fcst_hour -= 1
            prev_postproc_path = compute_postproc_path(prev_cycle, region_cfg.code, wksp_path, prev_fcst_hour)
    
    if prev_fcst_hour > 0:
        prev_cycle_id = compute_fmda_id(prev_cycle, region_cfg.code) + f"-f{prev_fcst_hour:02d}"
    else:
        prev_cycle_id = compute_fmda_id(prev_cycle, region_cfg.code)
    
    manifest_name = cycle_id + ".json"
    complete_manifest_name = f"fmda-{region_cfg.code}.json"
    if not is_cycle_computed(cycle, region_cfg, wksp_path, fcst_hour) and not osp.exists(prev_postproc_path):
        logging.warning(
            f"CYCLER postprocessing misses information for cycle {cycle} and forecast hour {fcst_hour}"
        )
        return None

    scalar_vars = [ "HGT", "T2", "RH", "PRECIP", "SNOWH", "EQUILd FM", "EQUILw FM", "WINDSPD", "SMOKE" ]
    vector_vars = [ "WINDVEC" ]

    esmf_cycle = utc_to_esmf(cycle + timedelta(hours=fcst_hour))
    mf = { "1" : {esmf_cycle : {}}}
    ensure_dir(osp.join(postproc_path, manifest_name))
    
    if not is_cycle_computed(cycle, region_cfg, wksp_path, fcst_hour):
        logging.info(
            f"CYCLER copying postprocessing from cycle {prev_cycle}-f{prev_fcst_hour} " 
            f"to cycle {cycle}-f{fcst_hour}"
        )
        prev_manifest_name = prev_cycle_id + ".json"
        prev_manifest_path = osp.join(prev_postproc_path, prev_manifest_name)
        if not osp.exists(prev_manifest_path):
            logging.error("CYCLER previous post-processing could not be found")
            return None
        prev_mf = json.load(open(prev_manifest_path, "r")) 
        prev_esmf_cycle = utc_to_esmf(prev_cycle + timedelta(hours=fcst_hour))
        for name in prev_mf["1"][prev_esmf_cycle].keys():
            prev_raster_name = prev_mf["1"][prev_esmf_cycle][name]["raster"]
            prev_cb_name = prev_mf["1"][prev_esmf_cycle][name]["colorbar"]
            raster_name = f"{cycle_id}-{name}-raster.png"
            cb_name = f"{cycle_id}-{name}-raster-cb.png"
            coords = prev_mf["1"][prev_esmf_cycle][name]["coords"]
            alpha = prev_mf["1"][prev_esmf_cycle][name].get("alpha",None)
            force_copy(
                osp.join(prev_postproc_path, prev_raster_name),osp.join(postproc_path, raster_name)
            )
            force_copy(
                osp.join(prev_postproc_path, prev_cb_name),osp.join(postproc_path, cb_name)
            )
            if alpha:
                mf["1"][esmf_cycle][name] = {
                    "raster" : raster_name, "coords" : coords, 
                    "colorbar" : cb_name, "alpha" : alpha
                }
            else:
                mf["1"][esmf_cycle][name] = {
                    "raster" : raster_name, "coords" : coords, 
                    "colorbar" : cb_name 
                }
    else:
        if bounds is None:
            bounds = (
                region_cfg.bbox[1], region_cfg.bbox[3], 
                region_cfg.bbox[0], region_cfg.bbox[2]
            )
        # read in the longitudes and latitudes
        geo_path = osp.join(wksp_path, "{}-geo.nc".format(region_cfg.code))
        logging.info(f"CYCLER reading longitudes and latitudes from NetCDF file {geo_path}")
        d = netCDF4.Dataset(geo_path)
        lats = d.variables["XLAT"][:,:]
        lons = d.variables["XLONG"][:,:]
        d.close()
        # read and process model variables
        with netCDF4.Dataset(model_path) as d:
            for name in scalar_vars:
                if name in d.variables.keys():
                    wisdom = get_wisdom(name).copy()
                    raster_png, coords, cb_png, levels = scalar_field_to_raster(
                        d.variables[name][:,:], lats, lons, wisdom
                    )
                    write_postprocess(
                        mf, postproc_path, cycle_id, esmf_cycle, name, 
                        raster_png, coords, cb_png, levels, 0.5
                    )
            for name in vector_vars:
                wisdom = get_wisdom(name).copy()
                if region_cfg.code == 'CONUS':
                    wisdom.update({"ref": 20})
                else:
                    wisdom.update({"ref": 10})
                c1_name, c2_name = wisdom["components"]
                if c1_name in d.variables.keys() and c2_name in d.variables.keys():
                    c1 = d.variables[c1_name][:,:]
                    c2 = d.variables[c2_name][:,:]
                    raster_png, coords = vector_field_to_raster(
                        c1, c2, lats, lons, wisdom
                    )
                    write_postprocess(
                        mf, postproc_path, cycle_id, esmf_cycle, name, 
                        raster_png, coords, alpha = 0.5
                    ) 
            fuel_classes = [(0, "1-hr DFM"), (1, "10-hr DFM"), (2, "100-hr DFM"), (3, "1000-hr DFM")]
            for i,name in fuel_classes:
                wisdom = get_wisdom("dfm").copy()
                fm_wisdom = wisdom
                fm_wisdom["name"] = f"Estimated {name}"
                raster_png, coords, cb_png, levels = scalar_field_to_raster(
                    d.variables["FMC_GC"][:,:,i], lats, lons, fm_wisdom
                )
                write_postprocess(
                    mf, postproc_path, cycle_id, esmf_cycle, name, 
                    raster_png, coords, cb_png, levels, 0.5
                )
               
            # Add other variables (Fire Weather Indices)
            # Get Temperature [K]
            T = d.variables["T2"][:]
            # Get Relative humidity [1]
            rh = d.variables["RH"][:] / 100  
            # Get Wind speed [m/s]
            ws = d.variables["WINDSPD"][:]
            # Calculate Saturated vapor presure [Pa]
            pws = calculate_svp(T) 
            # Calculate Vapor pressure deficit [hPa]         
            vpd = (1 - rh) * pws / 100
            # Fuel moisture (1-hour FM)
            fm1 = d.variables["FMC_GC"][:, :, 0] * 100
            # Calculate Moisture damping coefficient
            eta = calculate_eta(fm1)
            # Convert Wind speed to mph
            ws_mph = ws * 2.23694
            
            ### Fosberg Index (unitless) ###
            # Compute index
            ffwi = (eta * np.sqrt(1 + ws_mph**2)) / 0.3002
            # Visualization
            fosberg_wisdom = get_wisdom("FFWI").copy()
            raster_png, coords, cb_png, levels = scalar_field_to_raster(
                ffwi, lats, lons, fosberg_wisdom
            )
            write_postprocess(
                mf, postproc_path, cycle_id, esmf_cycle, "FFWI", 
                raster_png, coords, cb_png, levels, 0.5
            )
            
            ### HDW Index (kPa m s-1) ###
            # Compute index
            hdw = vpd * ws
            # Visualization
            hdw_wisdom = get_wisdom("HDWI").copy()
            raster_png, coords, cb_png, levels = scalar_field_to_raster(
                hdw, lats, lons, hdw_wisdom
            )
            write_postprocess(
                mf, postproc_path, cycle_id, esmf_cycle, "HDWI", 
                raster_png, coords, cb_png, levels, 0.5
            )
        if osp.exists("src/ingest/SynopticDB"):
            from ingest.SynopticDB.SynopticDB import SynopticDB
            import pandas as pd
            db = SynopticDB("ingest/SynopticDB")
            db.params["startDatetime"] = cycle - timedelta(hours=1)
            db.params["endDatetime"] = cycle
            db.params["minLongitude"] = bounds[0]
            db.params["maxLongitude"] = bounds[1]
            db.params["minLatitude"] = bounds[2] 
            db.params["maxLatitude"] = bounds[3]
            db.params["vars"] = ["fuel_moisture"]
            db.params["makeFile"] = False
            df,st = db.query_db()
            meso_wisdom = get_wisdom("dfm").copy()
            meso_wisdom["name"] = "MesoWest 10-hr DFM"
            meso_wisdom["bbox"] = bounds
            meso_wisdom["text"] = False
            if not (isinstance(df, pd.DataFrame) and len(df)) or fcst_hour > 0:
                logging.info("postprocess_cycle - missing Synoptic data, skipping")
                raster_png, coords, cb_png, levels = scatter_to_raster(
                    np.array([]), np.array([]), np.array([]), meso_wisdom
                )
            else:
                st = st.set_index("STID")
                data = df.groupby("STID").mean(numeric_only=True).join(st[["LONGITUDE","LATITUDE"]])
                fm = np.array(data["FUEL_MOISTURE_VALUE"]).astype(float) / 100.
                raster_png, coords, cb_png, levels = scatter_to_raster(
                    fm, np.array(data["LATITUDE"]).astype(float), 
                    np.array(data["LONGITUDE"]).astype(float), meso_wisdom
                ) 
            name = "MESO 10-hr DFM"
            write_postprocess(
                mf, postproc_path, cycle_id, esmf_cycle, name, 
                raster_png, coords, cb_png, levels, 1.0
            )

    prev_complete_manifest_path = osp.join(prev_postproc_path, complete_manifest_name)
    complete_manifest_path = osp.join(postproc_path, complete_manifest_name)
    manifest_path = osp.join(postproc_path, manifest_name)
    logging.info(f"writing manifest file {manifest_path}")
    json.dump(mf, open(manifest_path, "w"), indent=1, separators=(",", ":"))
    logging.debug(json.dumps(mf))
    if osp.exists(prev_complete_manifest_path):
        complete_mf = json.load(open(prev_complete_manifest_path, "r"))
        complete_mf["1"].update(mf["1"])
        json.dump(complete_mf, open(complete_manifest_path, "w"), indent=1, separators=(",", ":"))
    else:
        json.dump(mf, open(complete_manifest_path, "w"), indent=1, separators=(",", ":"))

    return postproc_path

def postprocess_cycle_rnn(cycle, region_cfg, wksp_path, fcst_hour, bounds=None):
    """
    Build rasters from the computed fuel moisture RNN prediction. 
    NOTE: As of June 29 2026, this is running after postprocessing onbly in forecast mode. This doubles up reading of weather vars like temp and ws
    Skipping over other covariates, starting with just FM10 raster, then other fuel classes and Indices

    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :param bounds: bounding box of the post-processing
    :return: the postprocessing path
    """
    model_path = compute_model_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    rnn_path = compute_rnn_path(cycle, region_cfg.code, wksp_path, fcst_hour=0) # Saving all relative to 0 hr cycle for now, in future if we do cyclical RNN prediction this will change
    model_path = compute_model_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    cycle_id = compute_fmda_id(cycle, region_cfg.code)
    if fcst_hour >0:
        cycle_id += f"-f{fcst_hour:02d}"
    postproc_path = compute_postproc_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    manifest_name = cycle_id + ".json"
    esmf_cycle = utc_to_esmf(cycle + timedelta(hours=fcst_hour))
    mf = { "1" : {esmf_cycle : {}}}
    ensure_dir(osp.join(postproc_path, manifest_name))
    # TODO: add check for already existing postproc here
    if False:
        pass
    else:
        if bounds is None:
            bounds = (
                region_cfg.bbox[1], region_cfg.bbox[3],
                region_cfg.bbox[0], region_cfg.bbox[2]
            )
        # read in the longitudes and latitudes
        geo_path = osp.join(wksp_path, "{}-geo.nc".format(region_cfg.code))
        logging.info(f"CYCLER reading longitudes and latitudes from NetCDF file {geo_path}")
        gd = netCDF4.Dataset(geo_path)
        lats = gd.variables["XLAT"][:,:]
        lons = gd.variables["XLONG"][:,:]       
        # read weather vars for FFWI
        with netCDF4.Dataset(model_path) as dat: 
            T = dat.variables["T2"][:]
            rh = dat.variables["RH"][:]/100
            ws = dat.variables["WINDSPD"][:]        
        # read and process RNN Predictions
        with netCDF4.Dataset(rnn_path) as rnn:
            fuel_classes = [(0, "1-hr DFM"), (1, "10-hr DFM"), (2, "100-hr DFM"), (3, "1000-hr DFM")]
            fuel_classes = [(0, "1-hr DFM"), (1, "10-hr DFM")]
            shortname = ["FM1", "FM10"]
            for i,name in fuel_classes:
                wisdom = get_wisdom("dfm").copy()
                fm_wisdom = wisdom
                fm_wisdom["name"] = f"Estimated {name} (RNN)"
                time_var = rnn.variables["time"]
                times = num2date(
                    time_var[:],
                    units=time_var.units,
                    calendar=time_var.calendar,
                )
                target_time = cycle + timedelta(hours=fcst_hour)
                tindex = np.where(times == target_time)[0].item()
                raster_png, coords, cb_png, levels = scalar_field_to_raster(
                    rnn.variables[shortname[i]][:,:,tindex]*1/100, lats, lons, fm_wisdom
                )
                write_postprocess(
                    mf, postproc_path, cycle_id, esmf_cycle, f"{name} (RNN)",
                    raster_png, coords, cb_png, levels, 0.5
                )
            # Compute FFWI
            fm1 = rnn.variables['FM1'][:,:,tindex]*1/100
            eta = calculate_eta(fm1)
            ws_mph = ws * 2.23694
            ffwi = (eta * np.sqrt(1 + ws_mph**2)) / 0.3002
            # Visualization
            fosberg_wisdom = get_wisdom("FFWI").copy()
            raster_png, coords, cb_png, levels = scalar_field_to_raster(
                ffwi, lats, lons, fosberg_wisdom
            )
            write_postprocess(
                mf, postproc_path, cycle_id, esmf_cycle, "FFWI (RNN)",
                raster_png, coords, cb_png, levels, 0.5
            )
            




def compute_fmda_id(cycle, region_code):
    """
    Construct a fmda id unique for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :return: a unique fmda id
    """
    time_stamp = cycle.strftime("%Y%m%d")
    fmda_id = f"fmda-{region_code}-{time_stamp}-{cycle.hour:02d}"
    return fmda_id

def compute_cycle_path(cycle, region_code, wksp_path):
    """
    Construct a relative path to the cycle path for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :return: a relative path (w.r.t. workspace and region) of the cycle path
    """
    fmda_id = compute_fmda_id(cycle, region_code)
    year_month_folder = cycle.strftime("%Y%m")
    return osp.join(wksp_path, region_code, year_month_folder, fmda_id)

def compute_model_path(cycle, region_code, wksp_path, fcst_hour=0, ext="nc"):
    """
    Construct a relative path to the fuel moisture model file
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :param fcst_hour: forecast hour
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    fmda_id = compute_fmda_id(cycle, region_code)
    cycle_path = compute_cycle_path(cycle, region_code, wksp_path)
    if fcst_hour != 0:
        filename = f"{fmda_id}-f{fcst_hour:02d}.{ext}"
    else:
        filename = f"{fmda_id}.{ext}" 
    return osp.join(cycle_path, filename)


def compute_rnn_path(cycle, region_code, wksp_path, fcst_hour=0, ext="nc"):
    """
    Construct a relative path to the fuel moisture RNN predictions file
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :param fcst_hour: forecast hour
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    fmda_id = compute_fmda_id(cycle, region_code)
    cycle_path = compute_cycle_path(cycle, region_code, wksp_path)
    ymd = cycle.strftime('%Y%m%d')
    hr = cycle.strftime('%H')
    to_ymd = to_utc.strftime('%Y%m%d')
    to_hr = to_utc.strftime('%H')
    filename = f"rnn_preds_{ymd}-{hr}_{to_ymd}-{to_hr}.{ext}"
    return osp.join(cycle_path, filename)

def compute_postproc_path(cycle, region_code, wksp_path, fcst_hour=0, ext="nc"):
    """
    Construct a relative path to the post-processing folder
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :param fcst_hour: forecast hour
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    year_month_folder = cycle.strftime("%Y%m")
    cycle_id = compute_fmda_id(cycle, region_code)
    if fcst_hour != 0:
        fmda_id = f"{cycle_id}-f{fcst_hour:02d}"
    else:
        fmda_id = f"{cycle_id}" 
    return osp.join(wksp_path, year_month_folder, fmda_id)

def find_region_indices(glat,glon,minlat,maxlat,minlon,maxlon):
    """
    Find the indices i1:i2 (lat dimension) and j1:j2 (lon dimension)
    that contain the desired region (minlat-maxlat,minlon-maxlon).

    :param glat: the grid latitudes
    :param glon: the grid longitudes
    :param minlat: the minimum latitude
    :param maxlat: the maximum latitude
    :param minlon: the minimum longitude
    :param maxlon: the maximum longitude
    :return: dim 0 min/max indices and dim1 min/max indices
    """
    i1, i2, j1, j2 = 0, glat.shape[0], 0, glat.shape[1]
    done = False
    while not done:
        done = True
        tmp = np.where(np.amax(glat[:, j1:j2],axis=1) < minlat)[0]
        if len(tmp):
            tmp = tmp[-1]
        else:
            tmp = i1
        if i1 != tmp:
            i1 = tmp
            done = False
        tmp = np.where(np.amin(glat[:, j1:j2],axis=1) > maxlat)[0]
        if len(tmp):
            tmp = tmp[0]
        else:
            tmp = i2
        if i2 != tmp:
            i2 = tmp
            done = False
        tmp = np.where(np.amax(glon[i1:i2,:],axis=0) < minlon)[0]
        if len(tmp):
            tmp = tmp[-1]
        else:
            tmp = j1
        if j1 != tmp:
            j1 = tmp
            done = False
        tmp = np.where(np.amin(glon[i1:i2,:],axis=0) > maxlon)[0]
        if len(tmp):
            tmp = tmp[0]
        else:
            tmp = j2
        if j2 != tmp:
            j2 = tmp
            done = False
    return i1,i2,j1,j2


def compute_hrrr_bounds(bbox):
    """
    Compute bounds from HRRR data even when HRRR data is not available from terrain static data
    
    :param bbox: the bounding box of the data
    :return: a tuple containing bound coordinates (min_lon,max_lon,min_lat,max_lat)
    """
    ds = netCDF4.Dataset("static/hrrr.terrainh.nc")
    lats,lons = ds["XLAT_M"][0], ds["XLONG_M"][0]
    i1, i2, j1, j2 = find_region_indices(lats, lons, bbox[0], bbox[2], bbox[1], bbox[3])
    lats,lons = lats[i1:i2,j1:j2], lons[i1:i2,j1:j2]
    return (lons.min(), lons.max(), lats.min(), lats.max())


def load_hrrr_data(grib_file, bbox):
    """
    Load relevant GRIB fields and return them
    
    :param grib_file: path to HRRR grib file
    :param bbox: the bounding box of the data
    :return: a dictionary with variables
    """
    gf = GribFile(grib_file)
    lats, lons = gf[1].latlons()
    date_str = str(gf[1].grb.date)
    date = datetime.strptime(date_str, '%Y%m%d')
    hrrr_vars = get_hrrr_var_list(date)
    if len(hrrr_vars) == 0:
        logging.warning("HRRR version is not provided, empty list of variables")
    # bbox format: minlat, minlon, maxlat, maxlon
    i1, i2, j1, j2 = find_region_indices(lats, lons, bbox[0], bbox[2], bbox[1], bbox[3])
    lats = lats[i1:i2,j1:j2] 
    lons = lons[i1:i2,j1:j2]
    hgt = np.ma.array(netCDF4.Dataset("static/hrrr.terrainh.nc")["HGT_M"][0])[i1:i2,j1:j2]
    data = {"lats": lats, "lons": lons, "hgt": hgt}
    
    for v,idx in hrrr_vars.items():
        if idx is not None:
            logging.info(gf[idx])
            var = np.ma.array(gf[idx].values())[i1:i2,j1:j2]
            logging.info("{} min {} max {}".format(v, np.min(var), np.max(var)))
            data.update({v: var})

    return data


def compute_equilibria(T, H):
    """
    Compute the drying and wetting equilibrium given temperature and relative humidity.
    
    :param T: the temperature at 2 meters in K
    :param H: the relative humidity in percent
    :return: a tuple containing the drying and wetting equilibrium
    """
    d = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    w = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    d *= 0.01
    w *= 0.01
    return d, w


def fmda_advance_region(cycle, cfg, grib_files, wksp_path, lookback_length, fcst_hour, meso_token, acquire=False):
    """
    Advance the fuel moisture estimates in the region specified by the configuration.
    The function assumes that the fuel moisture model has not been advanced to this
    cycle yet and will overwrite any previous computations.
    
    Control flow:
    
    1) read in HRRR variables
    2) check if there is a stored FM model for previous cycle
    2a) yes -> load it, advance one time-step, perform DA
    2b) no -> compute equilibrium, use background covariance to do DA
    3) store model
    
    :param cycle: the datetime indicating the processed cycle in UTC
    :param cfg: the configuration dictionary specifying the region
    :param grib_files: path to HRRR grib files to retrieve variables for this cycle (or previous)
    :param wksp_path: the workspace path for the cycler
    :param lookback_length: number of cycles to search before we find a computed cycle
    :param fcast_hour: if in forecast mode, the forecasting hour (0 otherwise)
    :param meso_token: the mesowest API access token or a list of them
    :param acquire: should the SynopticDB be updated? Normally only if CONUS
    :return: the model advanced and assimilated at the current cycle
    """
    min_num_obs = 10
    max_fm10_value = 0.5
    run_postprocessing = cfg.get("run_postprocessing", True)
    logging.info(f"hrrr_cycler.fmda_advance_region: cycle {cycle} and forecasts hour {fcst_hour}")
    model = None
    if fcst_hour == 0:
        prev_cycle = cycle - timedelta(hours=1)
        prev_fcst_hour = fcst_hour
        prev_model_path = compute_model_path(prev_cycle, cfg.code, wksp_path, prev_fcst_hour)
    else:
        prev_cycle = cycle
        prev_fcst_hour = fcst_hour - 1
        prev_model_path = compute_model_path(prev_cycle, cfg.code, wksp_path, prev_fcst_hour)
    if not osp.exists(prev_model_path):
        logging.info(
            f"CYCLER cannot find model from previous cycle {prev_cycle} "
            f"and forecasts hour {prev_fcst_hour}"
        )
        logging.info(f"CYCLER lookback length is {lookback_length}")
        if lookback_length > 0:
            model = fmda_advance_region(
                prev_cycle, cfg, grib_files, wksp_path, 
                lookback_length - 1, prev_fcst_hour, meso_token
            )
        elif fcst_hour > 0:
            model = fmda_advance_region(
                prev_cycle, cfg, grib_files, wksp_path, 
                lookback_length, prev_fcst_hour, meso_token
            )
    else:
        logging.info(f"CYCLER found previous model for cycle {prev_cycle} and forecasts hour {fcst_hour}.")
        model = FuelMoistureModel.from_netcdf(prev_model_path)
    
    grib_file = grib_files[lookback_length]
    # retrieve the variables and make sure they are available (we should not be here if they are not)
    if not osp.exists(grib_file):
        logging.warning("CYCLER could not find useable cycle.")
        logging.error(e)
        if run_postprocessing:
            logging.warning("CYCLER copying previous post-processing.")
            try:
                bounds = compute_hrrr_bounds(cfg.bbox)
                pp_path = postprocess_cycle(cycle, cfg, wksp_path, fcst_hour, bounds)   
                if pp_path != None:
                    if "shuttle_remote_host" in sys_cfg:
                        sim_code = "fmda-" + cfg.code
                        try:
                            send_product_to_server(
                                sys_cfg, pp_path, sim_code, sim_code,
                                sim_code + ".json", cfg.region_id + " FM"
                            )
                        except Exception as e:
                            logging.warning(
                                f"CYCLER failed sending to server. {sys_cfg['shuttle_remote_host']=}"
                            )
                            logging.warning("CYCLER exception {}".format(e))                    
            except Exception as e:
                logging.warning("CYCLER exception {}".format(e))
                logging.error("CYCLER skipping region {} for cycle {}".format(cfg.region_id,str(cycle)))
        sys.exit(1) 
    
    logging.info(f"CYCLER loading HRRR data from {grib_file}.")
    data = load_hrrr_data(grib_file, cfg.bbox)
    lats = data["lats"]
    lons = data["lons"]
    hgt = data["hgt"]
    Ed, Ew = compute_equilibria(data["t2"], data["rh"])
    
    rain = data["rain"][:,:] + 0
    # remove rain that is too small to make any difference 
    rain[rain < 0.01] = 0
    # remove bogus rain that is too large 
    rain[rain > 1e10] = 0
    # remove masked rain values
    rain[rain.mask] = 0

    dom_shape = data["t2"].shape
    # store the lons/lats for this domain
    geo_path = osp.join(wksp_path, "{}-geo.nc".format(cfg.code))
    if not osp.isfile(geo_path):
        logging.info(f"CYCLER initializing new file {geo_path}.")
        geo_path = ensure_dir(geo_path)
        d = netCDF4.Dataset(geo_path, "w", format="NETCDF4")
        d.createDimension("south_north", dom_shape[0])
        d.createDimension("west_east", dom_shape[1])
        xlat = d.createVariable("XLAT", "f4", ("south_north", "west_east"))
        xlat[:,:] = lats
        xlong = d.createVariable("XLONG", "f4", ("south_north", "west_east"))
        xlong[:,:] = lons
        d.close()
    else:
        logging.info(f"CYCLER file already exists:  {geo_path}.")

    # check if we must start from equilibrium
    if model is None:
        logging.info(f"CYCLER initializing from equilibrium for cycle {cycle}.")
        # setup model parameters    
        Tk = np.array([1.0, 10.0, 100.0, 1000.0]) * 3600
        m0 = np.expand_dims(0.5 * (Ed + Ew), axis=2)
        # background covariance
        P0 = np.diag([0.01, 0.01, 0.01, 0.01, 0.001, 0.001])
        model = FuelMoistureModel(m0[:,:,[0, 0, 0, 0]], Tk, P0)
    else:
        logging.info(f"CYCLER advancing model one hour to cycle {cycle}.")
        # always 1 hr step in HRRR
        dt = 3600
        # the process noise matrix
        Q = np.diag([1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-5])
        model.advance_model(Ed, Ew, rain, dt, Q)

    logging.info(f"CYCLER retrieving fm-10 observations for cycle {cycle}.")
    # no assimilation
    if fcst_hour > 0:
        logging.info("CYCLER forecasting mode, skipping data assimilation")
    else:
        # perform assimilation with mesowest observations
        tm_start = cycle - timedelta(minutes=30)
        tm_end = cycle + timedelta(minutes=30)
        if cfg.code == "CONUS" or acquire == True:
            fm10 = retrieve_mesowest_observations(meso_token, tm_start, tm_end, lats, lons, hgt, True)
        else:
            fm10 = retrieve_mesowest_observations(meso_token, tm_start, tm_end, lats, lons, hgt, False)
        
        logging.info(f"CYCLER filtering valid data from {len(fm10)} times")
        # filter fm10 values for statistics
        valid_times = [z for z in fm10.keys() if abs((z - cycle).total_seconds()) < 3600]
        fm10_filter = {}
        obs_valid_now = []
        for z in valid_times:
            vobs = [f for f in fm10[z] if f.obs_val > 0. and f.obs_val < max_fm10_value]
            obs_valid_now.extend(vobs)
            fm10_filter.update({z: vobs})
        fm10v = [obs.get_value() for obs in obs_valid_now]
        fm10 = fm10_filter
        
        if len(obs_valid_now) > min_num_obs:
            logging.info(
                f"CYCLER retrieved {len(fm10v)} valid observations at {len(valid_times)} unique times, "
                f"min/mean/max [{np.amin(fm10v):.3f}/{np.mean(fm10v):.3f}/{np.amax(fm10v):.3f}]."
            )
            # run the data assimilation step
            covs = [np.ones(dom_shape), hgt, lats, lons]
            covs_names = ["const", "hgt", "lat", "lon"]
            if np.any(rain > 0.01):
                covs.append(rain)
                covs_names.append("rain")
            if np.any(data["snow"] > 0.01):
                covs.append(data["snow"])
                covs_names.append("snow")
            other_covs = [
                "t2", "rh", "psfc", "soil_t_0", "soil_t_1", "soil_t_2", 
                "soil_t_3", "soil_t_4", "soil_t_5", "soil_t_6", "soil_t_7", 
                "soil_t_8", "soil_moist_0", "soil_moist_1", "soil_moist_2", 
                "soil_moist_3", "soil_moist_4", "soil_moist_5", "soil_moist_6", 
                "soil_moist_7", "soil_moist_8", "snowh", "spfh", "u10", "v10", 
                "ws", "ust", "znt", "swdown"
            ]
            for cov in other_covs:
                covs_names.append(cov)
                covs.append(data[cov])
                
            execute_da_step(model, cycle, covs, covs_names, fm10, use_lstsq=True)
    
    # make geogrid files for WPS; datasets and lines to add to GEOGRID.TBL
    geo_path = compute_model_path(cycle, cfg.code, wksp_path, fcst_hour, ext="geo")
    index = {
        "projection": "lambert",
        "dx" : 3000.0,
        "dy" : -3000.0,
        "truelat1" : 38.5,
        "truelat2" : 38.5,
        "stdlon" : 262.5,
        "radius" : 6370000.0
    }
    model.to_geogrid(geo_path, index, lats, lons)

    # make wps format files for WPS
    time_tag = cycle.strftime("%Y-%m-%d_%H") + f"f{fcst_hour:02d}"
    model.to_wps_format(osp.dirname(geo_path), index, lats, lons, time_tag)
    
    # store the new model  
    model_path = compute_model_path(cycle, cfg.code, wksp_path, fcst_hour)
    logging.info("CYCLER writing model variables to:  %s." % model_path)
    data.update({
        "EQUILd FM": Ed, "EQUILw FM": Ew, "PRECIP": rain, "HGT": hgt
    })
    rename_vars = {
        "t2": "T2", "rh": "RH", "snowh": "SNOWH", "ws": "WINDSPD",
        "u10": "U10", "v10": "V10", "SMOKE": "massden"
    }
    for orig_var_name,new_var_name in rename_vars.items():
        if orig_var_name in data:
            data[new_var_name] = data.pop(orig_var_name)
    model.to_netcdf(
        ensure_dir(model_path), data
    )

    if run_postprocessing:
        # create visualization and send results
        bounds = (lons.min(), lons.max(), lats.min(), lats.max())
        pp_path = postprocess_cycle(cycle, cfg, wksp_path, fcst_hour, bounds)   
        if pp_path != None:
            if "shuttle_remote_host" in sys_cfg:
                sim_code = "fmda-" + cfg.code
                try:
                    send_product_to_server(
                        sys_cfg, pp_path, sim_code, sim_code, 
                        sim_code + ".json", cfg.region_id + " FM"
                    )
                except Exception as e:
                    logging.warning(
                        f"CYCLER failed sending to server. {sys_cfg['shuttle_remote_host']=}"
                    )
                    logging.warning("CYCLER exception {}".format(e))
    
    return model
    
    
def is_cycle_computed(cycle, cfg, wksp_path, fcst_hour=0):
    """
    Check if the fuel model file exists (has been computed) for the
    cycle <cycle> and region configuration <cfg>.
    
    :param cycle: the cycle datetime in UTC
    :param cfg: the region configuration wrapped in a Dict for convenience
    :param wksp_path: the workspace path for the cycler
    :return: True if the model file has been found, False otherwise
    """
    path = compute_model_path(cycle, cfg.code, wksp_path, fcst_hour=fcst_hour)
    return osp.isfile(path)


def fmda_cycle_interval(start_cycle, end_cycle, conf_path=None):
    """
    Run historical cycle of fuel moisture estimates using FMDA for all the regions 
    specified by the configuration.
    
    :param start_cycle: initial time to create FMDA estimates
    :param end_cycle: final time to create FMDA estimates
    """
    if conf_path is None:
        conf = cfg
    else:
        conf = Dict(json.load(open(conf_path)))

    lookback_length = 0
    forecast_length = 0
    run_postprocessing = conf.get("run_postprocessing", True)
    hrrra = HRRRA(sys_cfg)
    cycle = start_cycle
    tstep = 0
    while cycle <= end_cycle:
        gribs = hrrra.retrieve_gribs(cycle, cycle)
        grib_files_anl = gribs["grib_files"]
        for region_id,region_cfg in conf.regions.items():
            logging.info(f"CYCLER processing region {region_id} for {cycle}")
            wrapped_cfg = Dict(region_cfg)
            wrapped_cfg.update({"region_id": region_id})
            wrapped_cfg.update({"forecast_length": forecast_length}) 
            wrapped_cfg.update({"run_postprocessing": run_postprocessing})
            # Process real-time
            if not is_cycle_computed(cycle, wrapped_cfg, conf.workspace_path):
                logging.info(f"CYCLER real-time processing for region {region_id} at cycle {cycle}")
                try:
                    fmda_advance_region(
                        cycle, wrapped_cfg, grib_files_anl, 
                        conf.workspace_path, lookback_length, 0, 
                        meso_token, True
                    )
                except Exception as e:
                    logging.warning(
                        f"CYCLER failed processing for region {region_id} at cycle {cycle}"
                    )
                    logging.warning("CYCLER exception {}".format(e))
            else:
                logging.info(
                    f"CYCLER already completed processing for region {region_id} "
                    f"at cycle {cycle}, skipping ..."
                )
        cycle += timedelta(hours=1)
        tstep += 1

def parse_bbox(args: Sequence[str]) -> tuple[float, float, float, float]:
    """
    Convert input arguments for bbox into numeric list
    """
    if len(args) != 4:
        raise ValueError(f"Expected 4 bbox values, got {len(args)}")

    try:
        return tuple(float(x) for x in args)
    except ValueError as e:
        raise ValueError(f"Invalid bbox values: {args}") from e


def run_checks():
    """
    Series of checks to run before executing cycler. Should stop process before spending a long time just to get nothing
    """
    import shutil

    if shutil.which("aws") is None:
        logging.error(
            "AWS CLI not found. Please install AWS CLI and ensure 'aws' is on your PATH."
        )
        return False    
    if not osp.exists("static/hrrr.terrainh.nc"):
        logging.error("Static HRRR terrain data doesn't exist at: static/hrrr.terrainh.nc")
    return


def get_rnn_dir(models_root, region_code):
    """
    Return the model directory for a GACC region code. Directory contains config, weights, and scaler

    Defaults to Rocky Mountain weights when the region code is unknown.
    """
    region_dirs = {
        "GBCC": "gb23-24",
        "RMCC": "rocky23-24",
        "SWCC": "sw23-24",
        "NRCC": "nr23-24",
    }
    default_code = "RMCC"
    code = str(region_code).upper()

    if code not in region_dirs:
        logging.warning(
            "No model weights directory configured for region code %s. "
            "Defaulting to Rocky Mountain weights (%s).",
            region_code,
            region_dirs[default_code],
        )
        code = default_code
    else:
        logging.info(
            "Using model weights directory for region code %s: %s",
            code,
            region_dirs[code],
        )
    return osp.join(models_root, region_dirs[code])

# Namelist converter
source_to_target = {
    'Ed': 'EQUILd FM',
    'Ew': 'EQUILw FM',
    'solar': 'swdown',
    'wind': 'WINDSPD',
    'elev': 'HGT',
    'rain': 'PRECIP',
    'lat': 'lats',
    'lon': 'lons',
    'solar': 'swdown'
} 

def build_analysis_paths(code, ts, wksp_dir = "wksp"):
    """
    """

    ym = ts.strftime("%Y%m")
    ymd = ts.strftime("%Y%m%d")
    hh = ts.strftime("%H")
    hpath = f"fmda-{code}-{ymd}-{hh}"
    filename = f"fmda-{code}-{ymd}-{hh}.nc"
    return osp.join(wksp_dir, code, ym, hpath, filename)

def build_fcst_paths(code, from_utc, fcst_hours, wksp_dir="wksp"):
    """
    """
    path0 = Path(build_analysis_paths(code, from_utc, wksp_dir))
    stem = path0.stem  # fmda-FIRE-20260518-16

    paths = [
        path0.with_name(f"{stem}-f{h:02d}.nc")
        for h in range(1, fcst_hours + 1)
    ]

    return [path0] + paths

def predict_auto_batch(model,
                       X,
                       batch_sizes=(16384, 8192, 4096, 2048, 1024, 512, 256, 128, 32),
                       verbose=1):
    """
    Predict using the largest batch size that fits in memory.

    NOTE: at this step for non-stateful model, batch size in predict is just a performance issue. The bigger the faster
    """
    last_exception = None

    for bs in batch_sizes:
        try:
            if verbose:
                print(f"Trying predict batch_size={bs}")
            preds = model.predict(X, batch_size=bs, verbose=verbose)
            if verbose:
                print(f"Success with batch_size={bs}")
            return preds
        except (MemoryError, tf.errors.ResourceExhaustedError) as e:
            last_exception = e
            if verbose:
                print(f"Failed with batch_size={bs}")

    raise RuntimeError(
        "All batch sizes failed during prediction."
    ) from last_exception

def to_netcdf(path, arr, valid_times, fuel_vars=("FM1", "FM10")):
    d = netCDF4.Dataset(path, "w", format="NETCDF4")
    
    arr = np.asarray(arr)
    ny, nx, nt, nfuel = arr.shape
    if len(fuel_vars) != nfuel: raise ValueError( f"len(fuel_vars)={len(fuel_vars)} must match last dim nfuel={nfuel}" )


    d.createDimension("south_north", ny)
    d.createDimension("west_east", nx)
    d.createDimension("time", nt)

    time = d.createVariable("time", "f8", ("time",))
    time.units = "hours since 1970-01-01 00:00:00 UTC"
    time.calendar = "standard"
    time[:] = netCDF4.date2num(
        valid_times.astype("datetime64[ms]").astype(object),
        units=time.units,
        calendar=time.calendar,
    )
    for i, var_name in enumerate(fuel_vars): 
        v = d.createVariable( var_name, "f4", ("south_north", "west_east", "time"), zlib=True, complevel=4, ) 
        v.long_name = f"{var_name} fuel moisture prediction" 
        v.coordinates = "time" 
        v[:] = arr[..., i]

    d.Conventions = "CF-1.8"

    d.close()


def warp_weights(weights0, bi_warp, bf_warp):
    """
    Given LSTM layer weights and time-warp parameters, return a new list
    of time-warped LSTM weights without modifying the input weights.
    """
    # Copy all arrays to avoid mutating the originals
    w_warped = [w.copy() for w in weights0]
    # Bias vector (Keras LSTM layout: [i, f, c, o])
    b = w_warped[2]
    # Infer number of LSTM units from bias length
    if b.ndim != 1 or b.shape[0] % 4 != 0:
        raise ValueError("Unexpected LSTM bias shape.")
    lstm_units = b.shape[0] // 4
    # Input gate biases (i)
    b[0:lstm_units] += bi_warp
    # Forget gate biases (f)
    b[lstm_units:2 * lstm_units] += bf_warp

    return w_warped

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    if len(sys.argv) == 1:
            mode = None
    if len(sys.argv) > 1:
        mode = sys.argv[1]
        if mode in ["a", "A", "f", "F"]:
            mode = mode.lower()
        else:
            mode = None
    if len(sys.argv) == 3:
        code = sys.argv[2]
        for k,region in cfg.regions.items():
            if region["code"] == code:
                cfg.regions = {
                    k: region
                }
                break
    elif len(sys.argv) == 6:
        mode = sys.argv[1]
        code = "FIRE"
        cfg.regions = {
            "Fire domain" : {
                "code" : code,
                "bbox" : parse_bbox(sys.argv[2:6])
            }
        }
        try:
            os.remove(osp.join(cfg.workspace_path,code+"-geo.nc"))
        except Exception as e:
            logging.warning(e)
        try:
            delete(osp.join(cfg.workspace_path,code))
        except Exception as e:
            logging.warning(e)
    
    if mode is None or len(cfg.regions) < 1:
        print("Usage: to use domains configured in etc/fmda_cycler.json")
        print(f"{sys.argv[0]} mode code")
        print("The supported modes are: analysis (a) and forecast (f)")
        print("To use a custom domain named FIRE by giving a bounding box:")
        print("./hrrr_cycler.sh mode lat1 lon1 lat2 lon2")
        print("Example: ./hrrr_cycler.sh a 42 -124.6 49 -116.4")
        exit(1) 

    # get parameters from configuration
    lookback_length = cfg.get("lookback_length", 24)
    forecast_length = cfg.get("forecast_length", 48)
    period_hours = cfg.get("period_hours", 6)
    run_postprocessing = cfg.get("run_postprocessing", True)
    run_rnn = cfg.get("run_rnn", False)
    # get more readable mode
    mode_name = "analysis" if mode == "a" else "forecast"
    # current time
    now = datetime.now(timezone.utc)
    cycle = (now - timedelta(minutes=59)).replace(minute=0, second=0, microsecond=0, tzinfo=None)
    # print statements
    logging.info(
        f"CYCLER activated at {now}, will attempt cycle at {cycle} with mode {mode_name}"
    )
    logging.info(f"regions: {json.dumps(cfg.regions)}")

    ### REAL-TIME SECTION
    # Get HRRR real-time data
    try:
        from_utc = cycle - timedelta(hours=lookback_length)
        to_utc = cycle
        hrrra = HRRRA(sys_cfg)
        anl_gribs = hrrra.retrieve_gribs(from_utc, to_utc) 
        grib_files_anl = anl_gribs["grib_files"]
    except:
        logging.warning(f"CYCLER could not find useable cycle {cycle}.")
        logging.warning("CYCLER copying previous post-processing.")
        for region_id,region_cfg in cfg.regions.items():
            wrapped_cfg = Dict(region_cfg)
            wrapped_cfg.update({"region_id": region_id})
            wrapped_cfg.update({"forecast_length": forecast_length}) 
            wrapped_cfg.update({"run_postprocessing": run_postprocessing})
            if run_postprocessing:
                try:
                    bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                    pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, bounds=bounds)
                    if pp_path != None:
                        if "shuttle_remote_host" in sys_cfg:
                            sim_code = "fmda-" + wrapped_cfg.code
                            try:
                                send_product_to_server(
                                    sys_cfg, pp_path, sim_code, sim_code,
                                    sim_code + ".json", cfg.region_id + " FM"
                                )
                            except Exception as e:
                                logging.warning(
                                    f"CYCLER failed sending to server. {sys_cfg['shuttle_remote_host']=}"
                                )
                                logging.warning("CYCLER exception {}".format(e))                        
                except Exception as e:
                    logging.warning("CYCLER exception {}".format(e))
                    logging.error(f"CYCLER skipping region {region_id} for cycle {cycle}")
        sys.exit(1)
    logging.info(f"have necessary HRRR data for cycle {cycle}.")
    
    # Start processing every region
    for region_id,region_cfg in cfg.regions.items():
        logging.info(f"CYCLER processing region {region_id} for {cycle}")
        wrapped_cfg = Dict(region_cfg)
        wrapped_cfg.update({"region_id": region_id})
        wrapped_cfg.update({"forecast_length": forecast_length}) 
        wrapped_cfg.update({"run_postprocessing": run_postprocessing})
        # Process real-time
        if not is_cycle_computed(cycle, wrapped_cfg, cfg.workspace_path):
            logging.info(f"CYCLER real-time processing for region {region_id} at cycle {cycle}")
            try:
                fmda_advance_region(
                    cycle, wrapped_cfg, grib_files_anl, 
                    cfg.workspace_path, lookback_length, 0, meso_token
                )
            except Exception as e:
                logging.warning(
                    f"CYCLER failed real-time processing for region {region_id} at cycle {cycle}"
                )
                logging.warning("CYCLER exception {}".format(e))
                logging.warning("CYCLER copying previous post-processing or re-trying.")
                if run_postprocessing:
                    try:
                        bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                        pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, 0, bounds=bounds)   
                        if pp_path != None:
                            if "shuttle_remote_host" in sys_cfg:
                                sim_code = "fmda-" + wrapped_cfg.code
                                try:
                                    send_product_to_server(
                                        sys_cfg, pp_path, sim_code, sim_code,
                                        sim_code + ".json", cfg.region_id + " FM"
                                    )
                                except Exception as e:
                                    logging.warning(
                                        f"CYCLER failed sending to server. {sys_cfg['shuttle_remote_host']=}"
                                    )
                                    logging.warning("CYCLER exception {}".format(e))
                    except Exception as e:
                        logging.error(
                            f"CYCLER skipping region {region_id} for cycle {cycle} and mode {mode_name}"
                        )
        else:
            logging.info(
                f"CYCLER already completed real-time processing for region {region_id} "
                f"at cycle {cycle}, skipping ..."
            )
        
    ### FORECASTING SECTION
    if mode == "f":
        logging.info(f"CYCLER forecasting region {region_id} at cycle {cycle}") 
        hrrr = HRRR(sys_cfg)
        shift_hours = cycle.hour % period_hours
        cycle_start = (cycle - timedelta(hours=shift_hours)).replace(tzinfo=timezone.utc)
        from_utc = cycle.replace(tzinfo=timezone.utc)
        to_utc = (cycle_start + timedelta(hours=forecast_length)).replace(tzinfo=timezone.utc)
        fcst_hour = 1
        tmp_utc = from_utc + timedelta(hours=fcst_hour)
        while tmp_utc <= to_utc:
            # Get HRRR forecast data
            try:
                fct_gribs = hrrr.retrieve_gribs(tmp_utc, tmp_utc, cycle_start=cycle_start)
                grib_files_fct = fct_gribs["grib_files"]
            except:
                logging.warning("CYCLER could not find useable cycle.")
                logging.warning("CYCLER copying previous post-processing.")
                for region_id,region_cfg in cfg.regions.items():
                    logging.info(f"CYCLER processing region {region_id} for {cycle}")
                    wrapped_cfg = Dict(region_cfg)
                    wrapped_cfg.update({"region_id": region_id})
                    wrapped_cfg.update({"forecast_length": forecast_length}) 
                    wrapped_cfg.update({"run_postprocessing": run_postprocessing})
                    if run_postprocessing:
                        try:
                            bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                            pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, fcst_hour, bounds)
                            if pp_path != None:
                                if "shuttle_remote_host" in sys_cfg:
                                    sim_code = "fmda-" + wrapped_cfg.code
                                    try:
                                        send_product_to_server(
                                            sys_cfg, pp_path, sim_code, sim_code,
                                            sim_code + ".json", cfg.region_id + " FM"
                                        )
                                    except Exception as e:
                                        logging.warning(
                                            f"CYCLER failed sending to server. {sys_cfg['shuttle_remote_host']=}"
                                        )
                                        logging.warning("CYCLER exception {}".format(e))                                
                        except Exception as e:
                            logging.warning(f"CYCLER exception {e}")
                            logging.error(
                                f"CYCLER skipping region {region_id} for cycle {cycle} "
                                f"and forecast hour {fcst_hour}"
                            )
                continue
            logging.info(
                f"have necessary HRRR data for forecasting cycle {cycle} "
                f"at forecast hour {fcst_hour}."
            )
            # Start processing every region
            for region_id,region_cfg in cfg.regions.items():
                logging.info(f"CYCLER processing region {region_id} for {cycle}")
                wrapped_cfg = Dict(region_cfg)
                wrapped_cfg.update({"region_id": region_id})
                wrapped_cfg.update({"forecast_length": forecast_length}) 
                wrapped_cfg.update({"run_postprocessing": run_postprocessing})
                # Processing forecast
                if not is_cycle_computed(cycle, wrapped_cfg, cfg.workspace_path, fcst_hour=fcst_hour):
                    logging.info(
                        f"CYCLER forecasting region {region_id} at cycle {cycle} and "
                        f"forecast hour {fcst_hour}"
                    ) 
                    try:
                        fmda_advance_region(
                            cycle, wrapped_cfg, grib_files_fct,
                            cfg.workspace_path, 0, fcst_hour, meso_token
                        )
                    except Exception as e:
                        logging.warning(
                            f"CYCLER failed forecasting for region {region_id} at cycle {cycle} and "
                            f"forecast hour {fcst_hour}"
                        )
                        logging.warning(f"CYCLER exception {e}")
                        logging.warning("CYCLER copying previous post-processing or re-trying.")
                        if run_postprocessing:
                            try:
                                bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                                pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, fcst_hour, bounds)
                                if pp_path != None:
                                    if "shuttle_remote_host" in sys_cfg:
                                        sim_code = "fmda-" + wrapped_cfg.code
                                        try:
                                            send_product_to_server(
                                                sys_cfg, pp_path, sim_code, sim_code,
                                                sim_code + ".json", cfg.region_id + " FM"
                                            )
                                        except Exception as e:
                                            logging.warning(
                                                f"CYCLER failed sending to server. {sys_cfg['shuttle_remote_host']=}"
                                            )
                                            logging.warning("CYCLER exception {}".format(e)) 
                            except Exception as e:
                                logging.error(
                                    f"CYCLER skipping region {region_id} for cycle {cycle} and "
                                    f"forecast hour {fcst_hour}"
                                )
                else:
                    logging.info(
                        f"CYCLER already completed forecasting processing for region {region_id} "
                        f"at cycle {cycle} and forecast hour {fcst_hour}, skipping ..."
                    )

            fcst_hour += 1
            tmp_utc = from_utc + timedelta(hours=fcst_hour)

        # RNN Forecast
        if run_rnn:
            logging.info(f"Running RNN Forecast")
            for region_id,region_cfg in cfg.regions.items():
                wrapped_cfg = Dict(region_cfg)
                wrapped_cfg.update({"region_id": region_id})
                wrapped_cfg.update({"forecast_length": forecast_length})
                #rdir = get_rnn_dir(models_root=cfg.rnn_model_dir, region_code=wrapped_cfg.region_id)
                rdir = cfg.rnn_model_dir
                logging.info(f"Running RNN Forecast with trained model: {rdir}")
                region_params = Dict(read_yml(osp.join(rdir, "params.yaml")))
                region_params.update({'timesteps': None}) # Pred model uses flexible time dimension
                region_scaler = joblib.load(osp.join(rdir, "scaler.joblib"))
                region_weights_path = osp.join(rdir, 'rnn.weights.h5')                
                rnn = OperationalRNNPredictor.from_weights(region_params, region_weights_path)
                # Get geographic info
                geopath = osp.join(cfg.workspace_path, f"{wrapped_cfg.code}-geo.nc")
                geodat = netCDF4.Dataset(geopath)
                lats = geodat.variables["XLAT"][:,:]
                lons = geodat.variables["XLONG"][:,:]
                # Get timeseries of weather variables
                times = pd.date_range(from_utc, from_utc+timedelta(hours=forecast_length), freq="1h")
                if mode == "a":
                    ncpaths = [build_analysis_paths(wrapped_cfg.code, ts, cfg.workspace_path) for ts in times]
                elif mode == "f":
                    ncpaths = build_fcst_paths(wrapped_cfg.code, from_utc, forecast_length, cfg.workspace_path)
                ncpaths_exist = [osp.exists(path) for path in ncpaths]
                ncpaths = np.array(ncpaths)[ncpaths_exist]
                valid_times = times[ncpaths_exist]
                logging.info("Found %d nc hour files out of %d total.", len(ncpaths), len(ncpaths_exist))
                arrs = []
                features = [source_to_target.get(f, f) for f in region_params.features_list]
                ncvars = [feat for feat in features if feat not in {"lats", "lons", "hod", "doy"}]
                for ncp in ncpaths:
                    logging.info(f"Processing file {ncp}")
                    with netCDF4.Dataset(ncp) as tmp_dat:
                        arrs.append(np.stack([tmp_dat.variables[v][:] for v in ncvars], axis=-1)) 
                arr = np.stack(arrs, axis=2) # Make shape (ny, nx, ntime, nfeats)
                # Add derived time features and lat/lon
                ## NOTE: valid_times is real times, physical times used to extract time features
                ## ncpaths might be forecast times f01, f02 relative to start
                doy   = valid_times.dayofyear.to_numpy()
                doy3d = np.broadcast_to(doy[None, None, :, None], (*arr.shape[:2], arr.shape[2], 1))
                hod   = valid_times.hour.to_numpy()
                hod3d = np.broadcast_to(hod[None, None, :, None], (*arr.shape[:2], arr.shape[2], 1))
                lat3d = np.broadcast_to(lats[:, :, None, None], (*arr.shape[:3], 1))
                lon3d = np.broadcast_to(lons[:, :, None, None], (*arr.shape[:3], 1))
                arr = np.concatenate([arr, doy3d, hod3d, lat3d, lon3d], axis=-1)
                vlist = ncvars + ["doy", "hod"] + ["lats", "lons"]
                # Structure Inputs
                idx = [vlist.index(f) for f in features if f in vlist]
                X_gridded = arr[:, :, :, idx]  # shape (ny, nx, ntime, nfeat)
                # Change Eq units (x100 for %)
                iEd, iEw = region_params.features_list.index("Ed"), region_params.features_list.index("Ew")
                X_gridded[:,:,:,iEd] *= 100
                X_gridded[:,:,:,iEw] *= 100
                # Reshape to 2d table to apply scaler, flatten (xy) dimensions. 
                # Target shape (nx*ny*ntime, nfeats)
                X_flat = X_gridded.reshape(-1, X_gridded.shape[-1])  # shape (ny*ny*ntime, nfeat)
                X_scaled = region_scaler.transform(X_flat)    # shape (ny*nx*ntime, nfeat)
                # Reshape to 3d (ny*nx, ntime, nfeat)
                ny, nx, ntimes, nfeatures = X_gridded.shape
                nbatch = ny*nx
                assert X_scaled.shape[0] == nbatch * ntimes
                assert X_scaled.shape[1] == nfeatures
                X = X_scaled.reshape(nbatch, ntimes, nfeatures)
                logging.info(f"Predictor Data Shape: {X_gridded.shape}")
                logging.info(f"    (ny, nx): {(ny, nx)}\n    ntimes: {ntimes}\n    nfeatures:{nfeatures}")
                logging.info(f"RNN Input Shape: {X.shape}")
                # Predict FM10, utility function tries largest batch size param for perfomrance
                logging.info(f"Predicting FM10")
                preds10 = predict_auto_batch(rnn, X)
                # Predict FM1 with twarped
                ## NOTE Hard coding for seed 29 with bi and bf warps, make flexible TODO
                fm1_info = Path(osp.join(cfg.transfer_dir, "fm1_median_rep_report.txt")).read_text().splitlines()
                bs = {'bi': 5.0, 'bf': -1.25}
                #weights10 = rnn.get_layer("lstm").get_weights()
                lstm_layer = next(layer for layer in rnn.layers if layer.name.startswith("lstm"))
                weights10 = lstm_layer.get_weights()
                logging.info(f"Predicting FM1")
                logging.info(f"Time-warping with bi={bs['bi']}, bf={bs['bf']}")
                weights1 = warp_weights(weights10, bi_warp = bs["bi"], bf_warp = bs["bf"])
                lstm_layer.set_weights(weights1)
                #rnn.get_layer("lstm").set_weights(weights1) 
                preds1 = predict_auto_batch(rnn, X)
                preds_gridded = np.concatenate([
                            preds1.reshape(ny, nx, ntimes, 1),
                            preds10.reshape(ny, nx, ntimes, 1), 
                            ], axis=-1
                        )


                # Write output in the right directory
                outdir = Path(compute_model_path(cycle, wrapped_cfg.code, cfg.workspace_path, 0)).parent
                filename = f"rnn_preds_{from_utc.strftime('%Y%m%d')}-{from_utc.strftime('%H')}_{to_utc.strftime('%Y%m%d')}-{to_utc.strftime('%H')}.nc"
                logging.info(f"Writing RNN Predictions for forecast period to: {osp.join(outdir, filename)}")
                to_netcdf(osp.join(outdir, filename), preds_gridded, valid_times.tz_localize(None))

                # Run postprocessing with RNN Predictions
                # Reuse postprocess_cycle looping over fcst_hours
                bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                for fhr in range(0, fcst_hour-1):
                    pp_path = postprocess_cycle_rnn(cycle, wrapped_cfg, cfg.workspace_path, fhr, bounds)


    # done
    logging.info(f"CYCLER cycle {cycle} complete with mode {mode_name}.")
