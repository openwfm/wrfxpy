# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict
import subprocess
import json
from datetime import datetime, timedelta
import pickle
import os
import os.path as osp
from pathlib import Path
import pandas as pd
import numpy as np
import sys
from osgeo import gdal, osr
from scipy.interpolate import griddata, RegularGridInterpolator
from synoptic.services import stations_timeseries

# Script used to generate starndardized dictionaries with fields necessary to run FMDA 
# User inputs location in the form of a RAWS STID code, start time, and end time
# Currently, combines 1) RAWS observations from either stash or MesoWest 2) HRRR data interpolated to station lat lon
# Future, add satellite data fields 


# Variables used throughout
sys_cfg = Dict(json.load(open('etc/conf.json')))
rawspath = "/home/hirschij/data" # path for raws data stash
hrrrpath = "/home/hirschij/data/hrrr/geotiff_files" # path for atmospheric data stash


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def within_past_year(dt):
    # Given a datetime object, return a logical indicating whether or not the input date is within the past year

    # Get the current date
    current_date = datetime.now()
    
    # Calculate the date one year ago
    one_year_ago = current_date - timedelta(days=365)
    
    # Check if the input datetime is within the past year
    return dt > one_year_ago

def read_raws_stash(stid, t, rawspath=rawspath):
    # Given RAWS station id and time, return fmda data from stash pickle file 
    
    # TODO: attempt to read from stash first, and if not present attempt MesoWest

    ## Format database stash file
    year = str(t.year)
    day = f"{t.timetuple().tm_yday:03d}" # stash organized with day of year out of 365
    hour = t.strftime("%H")
    f = osp.join(rawspath, "MesoDB", year, day, f"{year}{day}{hour}.pkl")

    ## Read from Stash if exists
    if osp.exists(f):
        return pd.read_pickle(f)
    else: 
        print(f"File does not exist: {f}")
        return None

def retrieve_raws_stash(stid, tstart, tend):
    ## Read Station location Info
    st = pd.read_pickle(osp.join(rawspath, "MesoDB/stations.pkl"))
    if stid in st.index.tolist():
        st = st.loc[stid]
    else:
        raise ValueError(f"Station {stid} not found in MesoDB/stations.pkl")
    
    ## Get hourly timestamps given input dates
    dates = pd.date_range(start=tstart,end=tend, freq="1H")

    # Initialize FM array as NaN
    fm = np.full(len(dates), np.nan, dtype=np.float64) # initialize fm array as NaN type float 64
    times = dates.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy() # initialize actual RAWS time as dates, replaced with actual time if data present
    for i in range(0, len(dates)):
        di = dates[i]
        dat = read_raws_stash(stid, di)
        dat = dat[dat.STID == stid] # limit to target STID
        #print(f"TEST: {dat}")
        if (dat.shape[0] == 1):
                #fm[i] = dat["fm10"][0]
                fm[i] = np.float64(dat.fm10.iloc[0])
                #times[i] = dat["datetime"][0].strftime('%Y-%m-%dT%H:%M:%SZ')
                times[i] = dat.datetime.iloc[0].strftime('%Y-%m-%dT%H:%M:%SZ')

    ## Read Station location Info
    st = pd.read_pickle(osp.join(rawspath, "MesoDB/stations.pkl"))
    st = st.loc[stid]

    ## Format Return Dictionaries
    loc = {
        "STID": stid,
        'lat' : st['LATITUDE'],
        'lon' : st['LONGITUDE'],
        'elev': st["ELEVATION"]
    }
    
    raws = {
        'time_raws': times,
        'fm': fm,
        'hours':len(fm)
    }

    return loc, raws

def retrieve_raws_api(stid, tstart, tend):
    # Given RAWS station ID, start and end times, return formatted FMDA dictionaries for location data and RAWS sensor data
    # Inputs:
    # stid: (str) RAWS station ID
    # tstart: (datetime) start time object
    # tend: (datetime) end time object
    # Returns: dataframe of data, output of synopticpy stations_timeseries

    params = dict(
        stid=[stid],
        vars=["air_temp", "relative_humidity", "precip_accum", "fuel_moisture", "wind_speed", "solar_radiation"],
        start=tstart,
        end=tend
    )

    df = stations_timeseries(**params)
    
    loc, raws = format_raws_df(df)

    return loc, raws

def format_raws_df(df):
    # Given input dataframe (the output of retrieve_raws_api), return formatted dictionary
    # Inputs:
    # df: (dataframe)
    # Returns: fmda dictionary

    ## Format Return Dictionaries
    loc = {
        "STID": df.attrs["STID"],
        'lat' : df.attrs['latitude'],
        'lon' : df.attrs['longitude'],
        'elev': df.attrs["ELEVATION"]
    }
    
    ## Extract times from dataframe index
    times = df.index.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy() # convert index to utc time
    ## Convert dataframe to dictionary
    raws = df.to_dict(orient = "list")
    
    # Convert lists to NumPy arrays
    raws = {key: np.array(value) for key, value in raws.items()}

    raws["time_raws"]=times
    raws["hours"]=len(times)
    
    ## Convert C to K 
    if df.attrs["UNITS"]["air_temp"] == "Celsius":
        print("Converting RAWS temp from C to K")
        raws["air_temp"] = raws["air_temp"]+273.15

    return loc, raws

def build_raws_dict(stid, tstart, tend, rawspath=rawspath, hrrrpath=hrrrpath):
    # Given RAWS station ID, start and end times, return formatted FMDA dictionaries for location data and RAWS sensor data
    # Inputs:
    # stid: (str) RAWS station ID
    # tstart: (datetime) start time object
    # tend: (datetime) end time object
    # rawspath: path to raws data stash MesoDB
    # hrrrpath: path to hrrr stash of fmda bands saved from grib files 
    #
    # Returns: two dictionaries, one for location specific data and one with observations of FMDA variables over given time period

    print("Retrieving RAWS data")
    if not (tstart <= tend):
        ValueError("Invalid time strings, start not before end.")
    
    if within_past_year(tstart):
        print(f"Time start {tstart} is within 1 year, retrieving all FMDA data with SynpoticPy")

        loc, raws = retrieve_raws_api(stid, tstart, tend)
    else:
        print(f"Time start {tstart} is NOT within 1 year, retrieving fuel moisture from RAWS stash")
        
        loc, raws = retrieve_raws_stash(stid, tstart, tend)

    
    ## Get hourly timestamps given input dates
    dates = pd.date_range(start=tstart,end=tend, freq="1H")
    raws["time"]=dates.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy()
    
    ## Standardize Inputs
    if "precip_accum" in raws.keys():
        raws["precip_accum"] = format_precip(raws["precip_accum"])
    
    name_mapping = {"air_temp":"temp", "fuel_moisture":"fm", "relative_humidity":"rh", "precip_accum":"rain","solar_radiation":"solar", "wind_speed":"wind_speed"}
    
    old_keys = [*raws.keys()]
    new_keys = []
    for key in old_keys:
        new_keys.append(name_mapping.get(key, key))
    old_keys = [*raws.keys()]
    new_keys = []
    for key in old_keys:
        new_keys.append(name_mapping.get(key, key))
    raws = dict(zip(new_keys, list(raws.values())))

    
    ## Calculate Equilibria if available
    if {"rain", "temp"} & set(raws.keys()):
        print("Calculating Equilibrium Moisture from RAWS data")
        raws['Ed'] = 0.924*raws['rh']**0.679 + 0.000499*np.exp(0.1*raws['rh']) + 0.18*(21.1 + 273.15 - raws['temp'])*(1 - np.exp(-0.115*raws['rh']))
        raws['Ew'] = 0.618*raws['rh']**0.753 + 0.000454*np.exp(0.1*raws['rh']) + 0.18*(21.1 + 273.15 - raws['temp'])*(1 - np.exp(-0.115*raws['rh']))
    
    return loc, raws

def format_precip(precipa):
    rain=np.array(precipa, dtype = 'float64')
    rain = np.diff(rain) # first difference to convert accumulated to hourly
    rain = np.insert(rain, 0, [np.NaN]) # add NaN entry to account for diff
    # Highest ever recorded hourly rainfall in inches is about 16: https://www.weather.gov/owp/hdsc_world_record
    rain[rain > 100] = np.NaN # filter out erroneously high
    rain[rain < 0] = np.NaN # filter out negative, results from diff function after precipa goes to zero
    return rain

def ts_at(interp_x, interp_y, values, method = "linear"):
    # Python implementation on regular grid of Jan methodology from https://github.com/openwfm/wrf-fire-matlab/blob/master/vis/ts_at.m
    interp_pts = np.array([interp_y, interp_x])

    # Get nearest neighbor
    center_x = round(interp_x)
    center_y = round(interp_y)

    # Build 3x3 grid around center, NOTE: xy flip in GDAL
    grid = np.meshgrid(np.array([center_y-1, center_y, center_y+1]),
            np.array([center_x-1, center_x, center_x+1]))
    grid = np.array([grid[0].flatten(), grid[1].flatten()]).T
    # Subset values
    value9 = values[
        grid[:,0],
        grid[:,1]]
    value9=value9.reshape(3,3)

    # print(f"Using method: {method}")

    interp = RegularGridInterpolator([np.array([center_y-1, center_y, center_y+1]),np.array([center_x-1, center_x, center_x+1])], value9)

    return interp(interp_pts, method=method)

def build_hrrr_dict(tstart, tend, lon, lat, hrrrpath=hrrrpath, method = 'linear', fmt = "%Y%m%d%H%M"):
    # tstart: (datetime)    start time
    # tend: (datetime)     end time
    # lon, lat      scalars, coordinates to interpolate to
    # hrrrpath      path to stored GRIB data, bands stored as geotiff files
    # method        only l1=nearest neighbor
    # fmt           date string format of tstart_str, tend_str

    # method (str): interpolation method, passed to RegularGridInterpolator
    # Internal Functions, will only reasonably be used within here

    def get_vals(tpath, pixel_x, pixel_y, method=method):
        # tpath (str): absolute path to tiff file
        # pixel_x, pixel_y grid coordinates after geotransform (would be index if at node)

        ds = gdal.Open(tpath)
        if(ds.RasterCount>1):
            print('Not Implemented for multiple Raster bands')
            sys.exit(-1)
        band = ds.GetRasterBand(1)
        data = band.ReadAsArray()

        # Interpolate
        vals = ts_at(pixel_x, pixel_y, data, method = method)

        return vals


    lon = np.float64(lon)
    lat = np.float64(lat)
    # Get dates
    dates = pd.date_range(start=tstart,end=tend, freq="1H")


    # Assemble time series, Start with gust Band, then reuse those projection calculations later
    bandnum=585
    d = dates[0]
    day_file = d.strftime("%Y%m%d") # HRRR data stash is in this format
    hour = d.strftime("%H")
    tpath = osp.join(hrrrpath, day_file, f"hrrr.t{hour}z.wrfprsf00.{bandnum}.tif")
    print("Opening: "+tpath)
    ds = gdal.Open(tpath)
    gt = ds.GetGeoTransform()
    gp = ds.GetProjection()
    

    if(ds.RasterCount>1):
        print('Not Implemented for multiple Raster bands')
        sys.exit(-1)


    # Get Projection info
    point_srs = osr.SpatialReference()
    point_srs.ImportFromEPSG(4326) # hardcode for lon/lat
    # GDAL>=3: make sure it's x/y
    # see https://trac.osgeo.org/gdal/wiki/rfc73_proj6_wkt2_srsbarn
    point_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    file_srs = osr.SpatialReference()
    file_srs.ImportFromWkt(gp)
    ct = osr.CoordinateTransformation(point_srs, file_srs)
    
    # Project given lon/lat to raster
    print(f"Target lon/lat: ({lon}, {lat})")
    mapx, mapy, z = ct.TransformPoint(np.float64(lon), np.float64(lat))
    gt_inv = gdal.InvGeoTransform(gt)
    pixel_x, pixel_y = gdal.ApplyGeoTransform(gt_inv, mapx, mapy)
    print(f"Projected pixel: ({pixel_x}, {pixel_y})")

    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()

    # Get wind (separate to not reload/transform data)
    wind_speed = np.zeros(len(dates)) # Array filled in loop
    time_hrrr = [""]*len(dates) # Initialize empty string list
    wind_speed[0] = get_vals(tpath, pixel_x, pixel_y)
    time_hrrr[0]=datetime.strftime(datetime.strptime(day_file, "%Y%m%d").replace(hour = int(hour)), "%Y-%m-%d %H:%M:%S")

    for i in range(1, len(dates)):
        d = dates[i]
        day_file = d.strftime("%Y%m%d") # HRRR data stash is in this format
        hour = d.strftime("%H")
        tpath = osp.join(hrrrpath, day_file, f"hrrr.t{hour}z.wrfprsf00.{bandnum}.tif")
        wind_speed[i] = get_vals(tpath, pixel_x, pixel_y)
        time_hrrr[i]=datetime.strftime(datetime.strptime(day_file, "%Y%m%d").replace(hour = int(hour)), "%Y-%m-%d %H:%M:%S")

    # Get other fields
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def build_timeseries(bandnum):
        loop_arr = np.zeros(len(dates)) # Array filled in loop
        for i in range(0, len(dates)):
            d = dates[i]
            day_file = d.strftime("%Y%m%d") # HRRR data stash is in this format
            hour = d.strftime("%H")
            tpath = osp.join(hrrrpath, day_file, f"hrrr.t{hour}z.wrfprsf00.{bandnum}.tif")
            loop_arr[i] = get_vals(tpath, pixel_x, pixel_y)
        return loop_arr

    temp = build_timeseries(616) # temp
    # Convert to Kelvin if in C, HRRR documentation says it should be in K but many series are not
    # Checking whether any value in temp array is less than 150, which is not a feasible K temp (equals -123 C)
    if np.any(temp < 150):
        print("Converting HRRR data temp from C to K")
        temp += 273.15

    rh = build_timeseries(620) # rh

    rain = build_timeseries(629) # accumulated precip

    solarDS = build_timeseries(661) # downward short wave solar

    solarDL = build_timeseries(662) # downward long wave solar

    solarDL = build_timeseries(662) # downward long wave solar

    solarUS = build_timeseries(663) # upward short wave solar

    solarUL = build_timeseries(664) # upward long wave solar

    # calculate Equilibriums from rh and temp
    print("Calculating Equilibrium Moisture from HRRR data")
    Ed = 0.924*rh**0.679 + 0.000499*np.exp(0.1*rh) + 0.18*(21.1 + 273.15 - temp)*(1 - np.exp(-0.115*rh))
    Ew = 0.618*rh**0.753 + 0.000454*np.exp(0.1*rh) + 0.18*(21.1 + 273.15 - temp)*(1 - np.exp(-0.115*rh))


    hrrr1 = {
        'time': dates.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy(),
        'time_hrrr': time_hrrr,
        'rain' : rain,
        'rh' : rh,
        'temp' : temp,
        'Ed' : Ed,
        'Ew' : Ew,
        'wind_speed' : wind_speed,
        'solarDS' : solarDS,
        'solarDL' : solarDL,
        'solarUS' : solarUS,
        'solarUL' : solarUL,
    }

    return hrrr1



# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <STID> <target_file>' % sys.argv[0]))
        print("Example: python src/ingest/build_fmda_dict.py 202106070800 202106300900 NWRU1 ~/testfile.pickle")
        sys.exit(-1)
    
    start = sys.argv[1]
    end = sys.argv[2]
    stid = sys.argv[3]
    outfile = sys.argv[4]
    
    # Build empty dictionary to fill with functions below
    t0 = datetime.strptime(start, "%Y%m%d%H%M")
    t1 = datetime.strptime(end, "%Y%m%d%H%M")
    title = f"{stid}_{t0.year}{t0.strftime('%m')}"
    out_dict = {title: dict.fromkeys(["loc", "RAWS", "HRRR"])}
    
    print(f"Building FMDA Dictionary at RAWS Site {stid} from {start} to {end}")
    print("~"*50)
    print("Hard-coded paths in the code. (Make more flexible in the future)")
    print(f"RAWS Stash location: {rawspath}")
    print(f"HRRR FMDA Bands Stash location: {hrrrpath}")
    
    print("~"*50)
    # get fm observations, site location data from RAWS, and atmospheric data from groundlevel sensors if available
    dict1, dict2 = build_raws_dict(stid, t0, t1)
    
    out_dict[title]["loc"] = dict1
    out_dict[title]["RAWS"] = dict2
    
    print("~"*50)
    out_dict[title]["HRRR"] = build_hrrr_dict(t0, t1, out_dict[title]['loc']['lon'], out_dict[title]['loc']['lat'])

    #print(d)
    print('Writing json to: ' + outfile)
    with open(outfile, 'wb') as handle:
        pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    







