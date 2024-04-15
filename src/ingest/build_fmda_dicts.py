# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict, hash2
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
import ast
from osgeo import gdal, osr
from scipy.interpolate import griddata, RegularGridInterpolator
from synoptic.services import stations_timeseries, stations_metadata

# Script used to generate starndardized dictionaries with fields necessary to run FMDA 
# User inputs location in the form of a lat/lon bounding box, start time, and end time
# Currently, combines 1) RAWS observations from either stash or SynopticPy 2) HRRR data interpolated to station lat lon
# Future Work: add satellite data fields 

# NOTE: bbox format in SynopticPy is [lonmin,latmin,lonmax,latmax]
#       bbox format in wrfxpy rtma_cycler_all is [latmin, lonmin, latmax, lonmax]

# Variables used throughout
sys_cfg = Dict(json.load(open('etc/conf.json')))
tokens = Dict(json.load(open('etc/tokens.json')))
meso_token = tokens.mesowest
rawspath = "/home/hirschij/data" # path for raws data stash
hrrrpath = "/data001/projects/hirschij/data/hrrr/geotiff_files/" # path for atmospheric data stash

# NOTE: choosing to exclude solar bands 'DLWRF', 'USWRF', 'ULWRF'
# Downward shortwave is expected theoretically to be the most useful solar field 
# RAWS have downward shortwave sensors, so these could be compared to model fields
band_df_hrrr = pd.DataFrame({
    'Band': [616, 620, 624, 628, 629, 661, 561, 612, 643],
    'hrrr_name': ['TMP', 'RH', "WIND", 'PRATE', 'APCP',
                  'DSWRF', 'SOILW', 'CNWAT', 'GFLUX'],
    'dict_name': ["temp", "rh", "wind", "rain", "precip_accum",
                 "solar", "soilm", "canopyw", "groundflux"],
    'descr': ['2m Temperature [K]', 
              '2m Relative Humidity [%]', 
              '10m Wind Speed [m/s]'
              'surface Precip. Rate [kg/m^2/s]',
              'surface Total Precipitation [kg/m^2]',
              'surface Downward Short-Wave Radiation Flux [W/m^2]',
              'surface Total Precipitation [kg/m^2]',
              '0.0m below ground Volumetric Soil Moisture Content [Fraction]',
              'Plant Canopy Surface Water [kg/m^2]',
              'surface Ground Heat Flux [W/m^2]']
})

# These are the variables I could find that correspond to the HRRR vars above, look into "variables()" from SynopticPy for more
# RAWS var 'Vertical Heat_Flux' might correspond to ground flux?
# Soil moisture has sparse availability and might never show up in final dictioanries
raws_vars = ["air_temp", "relative_humidity", "precip_accum", "fuel_moisture", "wind_speed", "solar_radiation", "soil_moisture"]

utc_format = "%Y-%m-%dT%H:%M:%SZ" # date format for UTC

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
    
    ## Get hourly timestamps given input times
    times = pd.date_range(start=tstart,end=tend, freq="1H")

    # Initialize FM array as NaN
    fm = np.full(len(times), np.nan, dtype=np.float64) # initialize fm array as NaN type float 64
    times = times.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy() # initialize actual RAWS time as dates, replaced with actual time if data present
    for i in range(0, len(times)):
        di = times[i]
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

def format_raws_df(df, tstart, tend):
    # Given input dataframe (the output of retrieve_raws_api), return formatted dictionary
    # Inputs:
    # df: (dataframe)
    # tstart: (datetime)
    # tend: (datetime)
    # Returns: tuple of dictionaries, location data and raws data (loc, raws)

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
    if "air_temp" in df.columns:
        if df.attrs["UNITS"]["air_temp"] == "Celsius":
            print("Converting RAWS temp from C to K")
            raws["air_temp"] = raws["air_temp"]+273.15
    
    ## Calculate Hourly Precipitation from accumulated
    if "precip_accum" in df.columns:
        print("Calculating hourly precipitation")
        raws["rain"] = format_precip(raws["precip_accum"])

    ## Format Names
    name_mapping = {"air_temp":"temp", "fuel_moisture":"fm", "relative_humidity":"rh", "precip_accum":"rain","solar_radiation":"solar", "wind_speed":"wind", "precip_accum":"precip_accum", "soil_moisture":"soil_moisture"}
    old_keys = [*raws.keys()]
    new_keys = []
    for key in old_keys:
        new_keys.append(name_mapping.get(key, key))
    old_keys = [*raws.keys()]
    new_keys = []
    for key in old_keys:
        new_keys.append(name_mapping.get(key, key))
    raws = dict(zip(new_keys, list(raws.values())))

    ## Add array of times requested, often different from returned time by a couple mins
    times = pd.date_range(start=tstart,end=tend, freq="1H")
    raws["time"]=times.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy()
    
    ## Calculate Equilibria if available
    if 'rh' in raws and 'temp' in raws:
        print("Calculating Equilibrium Moisture from RAWS data")
        raws['Ed'] = 0.924*raws['rh']**0.679 + 0.000499*np.exp(0.1*raws['rh']) + 0.18*(21.1 + 273.15 - raws['temp'])*(1 - np.exp(-0.115*raws['rh']))
        raws['Ew'] = 0.618*raws['rh']**0.753 + 0.000454*np.exp(0.1*raws['rh']) + 0.18*(21.1 + 273.15 - raws['temp'])*(1 - np.exp(-0.115*raws['rh']))
    else:
        print(f"Fields needed to calculate Equilibrium missing: {set(key for key in ['rh', 'temp'] if key not in raws)}")
        print("Equilibrium Moisture not calculated")

    return loc, raws

def format_precip(precipa):
    rain=np.array(precipa, dtype = 'float64')
    rain = np.diff(rain) # first difference to convert accumulated to hourly
    rain = np.insert(rain, 0, [np.NaN]) # add NaN entry to account for diff
    # Highest ever recorded hourly rainfall in inches is about 16: https://www.weather.gov/owp/hdsc_world_record
    rain[rain > 100] = np.NaN # filter out erroneously high
    rain[rain < 0] = np.NaN # filter out negative, results from diff function after precipa goes to zero
    return rain

# Function to return nested dictionary, 
# Top-level keys is station ID with start YYYYmm
# Next-level keys is location data and RAWS sensor data
def build_raws_dict(sts, tstart, tend):
    # Inputs:
    # sts: (df) dataframe of station data, output of stations_metadata
    out_dict = {} # set up return dictionary

    for st in sts:

        print("~"*50)
        print(f"Collecting RAWS data for {st}")
        params["stid"] = [st]
        try:
            dat = stations_timeseries(verbose="HIDE", **params)
    
            if "fuel_moisture" in dat.columns:
                print("Collected FMC data")
                loc, raws = format_raws_df(dat,tstart, tend)
                title = f"{st}_{t0.year}{t0.strftime('%m')}"
                out_dict[title] = {"loc":loc, "RAWS": raws}
            else:
                print("No FMC found for this station and time")
        except AssertionError as e:
            # Error handling behavior
            print("AssertionError caught:", e)
            
    return out_dict

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

def get_projection_info(ds, epsg = 4326):
    # Given a geotiff file (a HRRR band), 
    # return info necessary to transform lat/lon coords to the file structure
    # Inputs: 
    # ds: (osgeo.gdal.Dataset)
    # epsg: (int) default 4326 for lon/lat
    # Return: (tuple) with fields (ct, g_inv)
        # ct: (osgeo.osr.CoordinateTransformation)
        # gt_inv: (tuple) output of gdal.InvGeoTransform, also could be found with gdalinfo on command line
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
    gt_inv = gdal.InvGeoTransform(gt)

    return ct, gt_inv

def build_hrrr_path(d, band, fhour):
    # Inputs:
    # d: (datetime)
    # band: (int) HRRR band number
    # fhour: (int) forecast hour, eg "00" or "01"
    # Returns: (str) filepath to geotiff file
    day_file = d.strftime("%Y%m%d") # HRRR data stash is in this format
    hour = d.strftime("%H")
    tpath = osp.join(hrrrpath, day_file, f"hrrr.t{hour}z.wrfprsf{fhour}.{band}.tif")
    return tpath

def ftimes(n):
    # Helper function to generate list of strings used to identify forecast times
    # Input: integer n of number of forecast hours
    # Output: list of strings of form ["00", "01", "02", ...]
    
    # Generate numbers from 0 to n
    n = int(n)
    numbers = np.arange(n + 1)
    # Convert numbers to strings with leading zeros
    strings = np.char.zfill(numbers.astype(str), 2)
    # Convert NumPy array to a Python list
    l = strings.tolist()
    return l

def build_hrrr_dict(tstart, tend, dat, fhours, method = 'linear'):
    # tstart: (datetime)    start time
    # tend: (datetime)     end time
    # dat: (dict) dictionary, output of build_raws_dict
    # fhours: (int) number of hours into forecast
    # method: (str) interpolation method, passed to ts_at
    
    # Get times array
    times = pd.date_range(start=tstart,end=tend, freq="1H")
    # Forecast hours string formatting
    ft = ftimes(fhours)

    # Get Projection data from first band from band_df_hrrr,
    # reuse projection info for other bands
    # NOTE: this results in 1 extra read of geotiff files, but doing it for clarity
    d = times[0]
    band = band_df_hrrr.Band[0]
    tpath = build_hrrr_path(d, band, "00")
    print(f"Opening: {tpath}")
    if not osp.exists(tpath): 
        raise FileNotFoundError(f"The file '{tpath}' does not exist.")
    ds = gdal.Open(tpath)
    ct, gt_inv = get_projection_info(ds)
    ds = None # close connection
    print(f"Projection info collected: {gt_inv}")
    
    # Set up dictionary entries and projected pixel values
    # Format HRRR subdictionary with key for each band in band_df_hrrr, np.nan as placeholder
    # Add pixel_x and pixel_y to loc subdirectory to use with interpolation
    for k in dat.keys():
        dat[k]["HRRR"] = {"time": times.strftime(utc_format).to_numpy()} # Initialize HRRR subdir with times
        for ts in ft:
            dat[k]["HRRR"][f"f{ts}"]={}
            for name in band_df_hrrr.dict_name:
                dat[k]["HRRR"][f"f{ts}"][name] = np.full(len(times), np.nan, dtype=np.float64) # Initialize time series with np.nan
        lon = dat[k]["loc"]["lon"]
        lat = dat[k]["loc"]["lat"]
        print(f"Interpolating to RAWS {dat[k]['loc']['STID']}, target lat/lon {lat, lon}")
        mapx, mapy, z = ct.TransformPoint(np.float64(lon), np.float64(lat))
        pixel_x, pixel_y = gdal.ApplyGeoTransform(gt_inv, mapx, mapy)
        print(f"Projected pixel: ({pixel_x}, {pixel_y})")
        dat[k]["loc"]["pixel_x"] = pixel_x
        dat[k]["loc"]["pixel_y"] = pixel_y
        print("")
    
    # Loop over bands and forecast hours, build time series for each station in dictionary
    for index, row in band_df_hrrr.iterrows():
        print("~"*50)
        band = row["Band"]
        dict_name = row["dict_name"]
        print(f"Building Time Series for band: {band}, {row['descr']}")
        for i in range(0, len(times)):
            d = times[i]
            print(f"Collecting Data for date: {d}")
            for ts in ft:
                tpath = build_hrrr_path(d, band, ts)
                print(f"Opening: {tpath}")
                ds = gdal.Open(tpath)
                data = ds.GetRasterBand(1).ReadAsArray()
                # Loop over dictionary keys to interpolte for each station loc
                for k in dat:
                    pixel_x = dat[k]["loc"]["pixel_x"]
                    pixel_y = dat[k]["loc"]["pixel_y"]
                    interp_val = ts_at(pixel_x, pixel_y, data, method = method)
                    # Fill appropriate array value with interpolated value
                    dat[k]["HRRR"][f"f{ts}"][dict_name][i] = interp_val
                ds = None # close connection
            
    # Convert temp C to K
    # NOTE: this assumes simple celcius temps, which is what we have seen in HRRR 
        # But weather data may come in other formats
        # we should check the original grib files gdalinfo for temp being units C within the code
    # Then, Calculate Equilibria Moisture from HRRR Data and add description field
    print("~"*50)
    print("Calculating moisture Equilibria from rh and temp")
    for k in dat:
        for ts in ft:
            rh = dat[k]["HRRR"][f"f{ts}"]["rh"]
            temp = dat[k]["HRRR"][f"f{ts}"]["temp"]
            if np.any(temp < 150):
                temp += 273.15
                # Below if statement to only print once
                if k == [*out_dict.keys()][0]:
                    print("Converting HRRR data temp from C to K") 
            dat[k]["HRRR"][f"f{ts}"]["Ed"] = 0.924*rh**0.679 + 0.000499*np.exp(0.1*rh) + 0.18*(21.1 + 273.15 - temp)*(1 - np.exp(-0.115*rh))
            dat[k]["HRRR"][f"f{ts}"]["Ew"] = 0.618*rh**0.753 + 0.000454*np.exp(0.1*rh) + 0.18*(21.1 + 273.15 - temp)*(1 - np.exp(-0.115*rh))    
            dat[k]["HRRR"][f"f{ts}"]["descr"] = f"Source: HRRR data from 3d pressure model, linear grid interpolated to RAWS location"
    return dat

def parse_bbox(box_str):
    try:
        # Use ast.literal_eval to safely parse the string representation
        # This will only evaluate literals and avoids security risks associated with eval
        box = ast.literal_eval(box_str)
        # Check if the parsed box is a list and has four elements
        if isinstance(box, list) and len(box) == 4:
            return box
        else:
            raise ValueError("Invalid bounding box format")
    except (SyntaxError, ValueError) as e:
        print("Error parsing bounding box:", e)
        sys.exit(-1)
        return None

def call_curl(url):
    try:
        # Run the curl command and capture its output
        result = subprocess.run(['curl', url], capture_output=True, text=True, check=True)
        # Decode the JSON output
        json_data = json.loads(result.stdout)
        return json_data
    except subprocess.CalledProcessError as e:
        print("Error executing curl command:", e)
        return None
    except json.JSONDecodeError as e:
        print("Error decoding JSON:", e)
        return None

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 6:
        print(f"Invalid arguments. {len(sys.argv)} was given but 5 expected")
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <bbox> <forecast_hours> <target_file>' % sys.argv[0]))
        print("Example: python src/ingest/build_fmda_dicts.py 202401010000 202401010200 '[37,-105,39,-103]' 1 ~/testfile.pickle")
        print("bbox format should match rtma_cycler: [latmin, lonmin, latmax, lonmax]")
        sys.exit(-1)
    
    start = sys.argv[1]
    end = sys.argv[2]
    bbox = parse_bbox(sys.argv[3])
    forecast_hours = sys.argv[4]
    outfile = sys.argv[5]
    
    # Build empty dictionary to fill with functions below
    t0 = datetime.strptime(start, "%Y%m%d%H%M")
    t1 = datetime.strptime(end, "%Y%m%d%H%M")

    print(f"Building FMDA Dictionary for RAWS Sites within {bbox}, from {t0} to {t1}")
    print(f"Number of Forecast Hours: {forecast_hours}")
    print("~"*50)
    print("Hard-coded paths in the code. (Make more flexible in the future)")
    print(f"RAWS Stash location: {rawspath}")
    print(f"HRRR FMDA Bands Stash location: {hrrrpath}")
    print("Target HRRR Data Info:")
    print(band_df_hrrr)
    print("~"*50)
    
    # Convert bbox from wrfxpy format to synoptic
    bbox2 = [bbox[1], bbox[0], bbox[3], bbox[2]]
    
     # Get station availability within bbox                                                                                  # NOTE: SynopticPy has bug where stations_metadata where it doesn't properly use time parameters for station availability in a bbox. See Issue #55 SynopticPy on github
    url = f"https://api.synopticdata.com/v2/stations/metadata?bbox={bbox2[0]},{bbox2[1]},{bbox2[2]},{bbox2[3]}&vars=fuel_moisture&obrange={start},{end}&token={meso_token}"
    print(f"Attempting Synoptic retrieval from URL: https://api.synopticdata.com/v2/stations/metadata?bbox={bbox2[0]},{bbox2[1]},{bbox2[2]},{bbox2[3]}&vars=fuel_moisture&obrange={start},{end}&token=HIDDEN")
    command = f"curl -X GET '{url}'"
    sts_json = call_curl(url)
    sts = pd.DataFrame(sts_json["STATION"], index=[i["STID"] for i in sts_json["STATION"]])
    print(sts)
    sts = sts.transpose()

    print(f"Number of Stations within bbox: {sts.shape[0]}")
    # Parameter dictionary used within SynopticPy
    params = dict(
        stid=["PLACEHOLDER"], # change this in the loop
        vars=raws_vars,
        start=t0,
        end= t1+timedelta(hours=1) # add an hour since it doesn't include end date exactly
    )
    
    print("Retrieving RAWS Data")
    out_dict = build_raws_dict(sts, t0, t1)

    print("~"*50)
    print("Retrieving HRRR Data and Interpolating")
    out_dict = build_hrrr_dict(t0, t1, out_dict, forecast_hours)
    
    print(f"Writing output to {outfile}")
    # print(f"Output hash: {hash2(out_dict)}") 
    ## NOTE: the interpolation procedure returns very slightly different pixel values given the same lon/lat input
    ## This makes exact hashing not possible, but values are still the same to many decimal places
    ## TODO: check for reproducibility
    
    with open(outfile, 'wb') as handle:
        pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

