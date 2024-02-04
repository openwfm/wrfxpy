# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict
import subprocess
import json
from datetime import datetime
import pickle
import os
import os.path as osp
from pathlib import Path
import pandas as pd
import numpy as np
import sys
from osgeo import gdal, osr
from scipy.interpolate import griddata, RegularGridInterpolator

sys_cfg = Dict(json.load(open('etc/conf.json')))

# DataFrame to track needed Bands for FMDA and metadata
## Note there is a surface temp, but no surface RH. 
## Sticking with 2m temp and rh for consistency since those are combined to EQ
band_df_hrrr = pd.DataFrame({
    'Band': [585, 616, 620, 628, 629, 661, 662, 663, 664],
    'hrrr_name': ['GUST', 'TMP', 'RH', 'PRATE', 'APCP',
                  'DSWRF', 'DLWRF', 'USWRF', 'ULWRF'],
    'descr': ['surface Wind Speed (Gust) [m/s]',
              '2 m Temperature [K]', 
              '2 m Relative Humidity [%]', 
              'surface Precip. Rate [kg/m^2/s]',
              'surface Total Precipitation [kg/m^2]',
              'surface Downward Short-Wave Radiation Flux [W/m^2]',
              'surface Downward Long-Wave Rad. Flux [W/m^2]',
              'surface Upward Short-Wave Radiation Flux [W/m^2]',
              'Upward Long-Wave Rad. Flux [W/m^2]']
})


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def extract_raws_stash(t, stash_path = "/storage/math/NSF1/farguella/MesoDB.tar.gz", dest_path = "wksp"):
    # Given time string in UTC, extract RAWS stash pickle file to temporary path.
    # Inputs:
    # t: (pandas Timestamp) time in UTC, output of pandas date_range
    # stash_path: (str) absolute path to RAWS stash
    # dest_path: (str) path to extract files, default the wrfxpy/wksp folder
    # Returns: file path of extracted pickle file
    
    ## Format database stash file
    year = str(t.year)
    day = f"{t.timetuple().tm_yday:03d}" # stash organized with day of year out of 365
    hour = t.strftime("%H")
    stash_file = osp.join("MesoDB", year, day, f"{year}{day}{hour}.pkl") # build target file relative to stash path
    
    ## Extract stations file if not present
    if not osp.exists(osp.join(dest_path, "MesoDB/stations.pkl")):
        command = f"tar -xzvf {stash_path} -C {dest_path} MesoDB/stations.pkl" # Format tar subprocess
        subprocess.call(command, shell=True)
        print(f"Executing command: {command}")

    ## If file does not exist, extract
    r = osp.join(dest_path, stash_file)
    if not osp.exists(r):
        if stash_path.endswith("tar.gz"):
            command = f"tar -xzvf {stash_path} -C {dest_path} {stash_file}" # Format tar subprocess
            print(f"Executing command: {command}")
            subprocess.call(command, shell=True)
            print(f"Checking file {stash_file} exists at {dest_path}")
            osp.exists(r)
        else: raise ValueError("Invalid file path. The file must end with '.tar.gz'.")

    return r



def build_raws(tstart_str, tend_str, stid, datapath, fmt = "%Y%m%d%H%M", temppath = "wksp"):
    # Given times, stid, and path, return a dictionary of RAWS data fields (fm10, time, lat/lon, elevation) 

    print('~'*25)
    print('Collecting RAWS station: '+stid)

    ## Format with datetime
    tstart = datetime.strptime(tstart_str, fmt)
    tend = datetime.strptime(tend_str, fmt)
    dates = pd.date_range(start=tstart,end=tend, freq="1H")
    
    if not (tstart <= tend): 
        ValueError("Invalid time strings, start not before end.")

    # Initialize FM array as NaN
    fm = np.full(len(dates), np.nan, dtype=np.float64) # initialize fm array as NaN type float 64
    times0 = dates.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy() # passed to output dict
    times = dates.strftime('%Y-%m-%dT%H:%M:%SZ').to_numpy() # initialize actual RAWS time as dates, replaced with actual time if available 
    for i in range(0, len(dates)):
        d = dates[i]
        p = extract_raws_stash(d) # extract RAWS stash pickle file from DB
        dat = pd.read_pickle(p) # Read hourly pickle file
        dat = dat[dat.STID == stid] # limit to target STID
        if (dat.shape[0] == 1):
                fm[i] = dat["fm10"][0]
                times[i] = dat["datetime"][0].strftime('%Y-%m-%dT%H:%M:%SZ')
    
    ## Read STID Info
    st = pd.read_pickle(osp.join(temppath, "MesoDB/stations.pkl"))
    st = st.loc[stid]
    
    dict1 = {
        'time': times0,
        'time_raws': times,
        'STID' : stid,
        'fm' : fm,
        'title' : 'RAWS Station '+stid,
        'hours':len(fm),
        'lat' : st['LATITUDE'],
        'lon' : st['LONGITUDE'],
        'elev': st["ELEVATION"]
    }
    return dict1

def build_atm(tstart_str, tend_str, lon, lat, atmpath, atm_source="HRRR"):

    if atm_source == "HRRR": atm1 = build_hrrr(tstart_str, tend_str, lon, lat, atmpath)

    return atm1


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

def build_hrrr(tstart_str, tend_str, lon, lat, hrrrpath, method = 'linear', fmt = "%Y%m%d%H%M"):
    # tstart_str    start time as string 
    # tend_str      end time as string
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
    

    # Get dates
    tstart = datetime.strptime(tstart_str, fmt)
    tend = datetime.strptime(tend_str, fmt)
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
        print("Converting temp from C to K")
        temp += 273.15

    rh = build_timeseries(620) # rh
    
    rain = build_timeseries(629) # accumulated precip

    solarDS = build_timeseries(661) # downward short wave solar

    solarDL = build_timeseries(662) # downward long wave solar

    solarUS = build_timeseries(663) # upward short wave solar

    solarUL = build_timeseries(664) # upward long wave solar

    # calculate Equilibriums from rh and temp
    Ed = 0.924*rh**0.679 + 0.000499*np.exp(0.1*rh) + 0.18*(21.1 + 273.15 - temp)*(1 - np.exp(-0.115*rh))
    Ew = 0.618*rh**0.753 + 0.000454*np.exp(0.1*rh) + 0.18*(21.1 + 273.15 - temp)*(1 - np.exp(-0.115*rh))
    

    hrrr1 = {
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

def build_dictionary(tstart_str, tend_str, stid, rawspath, atmpath, dict_name = "test_dict", h2=200, fmt = "%Y%m%d%H%M"):
    dict1 = build_raws(tstart_str, tend_str, stid, rawspath, fmt)
    
    atm1 = build_atm(tstart_str, tend_str, dict1['lon'], dict1['lat'], atmpath)
    
    dict1['time_hrrr'] = atm1['time_hrrr'] ## MAKE flexible to do other data sources besides HRRR
    dict1['id'] = dict_name
    dict1['rain'] = atm1['rain']
    dict1['rh'] = atm1['rh']
    dict1['temp'] = atm1['temp']
    dict1['Ed'] = atm1['Ed']
    dict1['Ew'] = atm1['Ew']
    dict1['wind'] = atm1['wind_speed']
    dict1['solarDS'] = atm1['solarDS']
    dict1['solarDL'] = atm1['solarDL']
    dict1['solarUS'] = atm1['solarUS']
    dict1['solarUL'] = atm1['solarUL']
    dict1['h2']=h2 # default number of hours for training period, we have been testing h2=200

    return dict1

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 7:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <STID> <rawspath> <atmpath> <target_file>' % sys.argv[0]))
        print("Example: python src/ingest/build_hrrr_dict.py 202106070800 202106300900 NWRU1 /storage/math/NSF1/farguella/ /home/hirschij/data/hrrr/geotiff_files/ ~/testfile.pickle")
        sys.exit(-1)
    start = sys.argv[1]
    end = sys.argv[2]
    stid = sys.argv[3]
    rawspath = sys.argv[4] # path for raws data
    atmpath = sys.argv[5] # path for atmospheric data
    outfile = sys.argv[6]

    print('Building RAWS Training Dictionary')
    print('Start: ' + str(start))
    print('End: ' + str(end))

    dict1 = build_dictionary(start, end, stid, rawspath, atmpath)
    
    # Dictionary must be nested with case names in this way
    # TODO: make this script flexible to get more STIDs without horrible inputs
    out_dict = {}
    out_dict[stid]=dict1
    
    print(out_dict.keys())

    print('Writing json to: ' + outfile)
    with open(outfile, 'wb') as handle:
        pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

