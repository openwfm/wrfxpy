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

def build_raws(tstart_str, tend_str, stid, datapath, fmt = "%Y%m%d%H%M"):
    # Given times, stid, and path, return a dictionary of RAWS data fields (fm10, time, lat/lon, elevation) 

    print('~'*25)
    print('Collecting RAWS station: '+stid)

    ## Format with datetime
    tstart = datetime.strptime(tstart_str, fmt)
    tend = datetime.strptime(tend_str, fmt)
    dates = pd.date_range(start=tstart,end=tend, freq="1H")
    
    # Given datapath of RAWS stach files, create Dataframe with file path, time 1 and time 2
    df=pd.DataFrame({
        'File': os.listdir(datapath)
        })
    # Remove stations.pkl file
    mask = ~df['File'].str.contains('stations')
    df = df[mask]
    # Extract times from file name
    df['t0']=df['File'].apply(lambda f: Path(f).stem.split('_')[0])
    df['t1']=df['File'].apply(lambda f: Path(f).stem.split('_')[1])
    
    # Convert to Datetime
    df['d0']=df['t0'].apply(lambda t: datetime.strptime(t, fmt))
    df['d1']=df['t1'].apply(lambda t: datetime.strptime(t, fmt))
   
    # Arrange dataframe by date
    df = df.sort_values(by=['d0'])
    df = df.reset_index()

    # Find Files that cover given time period
    # latest date before/equal to tstart to
    # earliest date after/equal to tend
    ind1=np.where(df['d0'].eq(df['d0'][df['d0'] <= tstart].max()))[0][0]
    ind2=np.where(df['d0'].eq(df['d0'][df['d0'] >= tend].min()))[0][0]
    
    df2=df.iloc[ind1:ind2]
    df2=df2.reset_index()
    
    ## Open Pickle Files
    p = osp.join(datapath, df2['File'][0])
    f = open(p, 'br')
    data = pickle.load(f)

    ## Subset to given stid
    raws1 = data[data["STID"].eq(stid)]

    for i in range(1, len(df2)):
        p = osp.join(datapath, df2['File'][i])
        f = open(p, 'br')
        data = pickle.load(f)
        # Add to dataframe for given stid
        raws1 = pd.concat([raws1, data[data["STID"].eq(stid)]])
    
    raws1 = raws1.drop_duplicates()
    # Get Station data
    f = open(osp.join(datapath, 'stations.pkl'), 'br')
    st = pickle.load(f)
    st = st[st.index == stid]

    print(st)
    
    dict1 = {
        'time_raws': raws1['datetime'].to_numpy(),
        'STID' : raws1['STID'].unique()[0],
        'fm' : raws1['fm10'].to_numpy(),
        'title' : 'RAWS Station '+stid,
        'hours':len(raws1['datetime']),
        'lat' : st['LATITUDE'].to_numpy()[0],
        'lon' : st['LONGITUDE'].to_numpy()[0],
        'other': {'elev': st['ELEVATION']}
    }
    return dict1

def build_atm(tstart_str, tend_str, lon, lat, atmpath, atm_source="HRRR"):

    if atm_source == "HRRR": atm1 = build_hrrr(tstart_str, tend_str, lon, lat, atmpath)

    return atm1

def interp_l1(x, y):
    x = round(x)
    y = round(y)
    return x, y


def build_hrrr(tstart_str, tend_str, lon, lat, hrrrpath, method = 'l1', fmt = "%Y%m%d%H%M"):
    # tstart_str    start time as string 
    # tend_str      end time as string
    # lon, lat      scalars, coordinates to interpolate to
    # hrrrpath      path to stored GRIB data, bands stored as geotiff files
    # method        only l1=nearest neighbor
    # fmt           date string format of tstart_str, tend_str

    # method (str): interpolation method (NOTE: currently only L1 nearest neighbors aka rounding) 
    # Internal Functions, will only reasonably be used within here

    def get_vals(tpath, pixel_x, pixel_y, method='l1'):
        # tpath (str): absolute path to tiff file
        # pixel_x, pixel_y grid coordinates after geotransform (would be index if at node)

        ds = gdal.Open(tpath)
        if(ds.RasterCount>1):
            print('Not Implemented for multiple Raster bands')
            sys.exit(-1)
        band = ds.GetRasterBand(1)
        data = band.ReadAsArray()
        
        # Interpolate
        if method == "l1":
            x, y = interp_l1(pixel_x, pixel_y)
        
        return data[y, x]
    

    # Get dates
    tstart = datetime.strptime(tstart_str, fmt)
    tend = datetime.strptime(tend_str, fmt)
    dates = pd.date_range(start=tstart,end=tend, freq="1H")
    # dates = dates[0:100] ### NOTE: THIS LINE IS ONLY FOR TESTING

    # Assemble time series
    bandnum=585 # Start with gust Band, then reuse those projection calculations later
    tpath = osp.join(hrrrpath, "20210304",f"hrrr.t00z.wrfprsf00.{bandnum}.tif")
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
    for i in range(0, len(dates)):
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
        print("Example: python src/ingest/build_hrrr_dict.py 202106070800 202106300900 NWRU1 /storage/math/NSF1/jmandel/RAWSData/haguespeak/meso /home/hirschij/data/hrrr/geotiff_files/ ~/testfile.pickle")
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
    out_dict['case1']=dict1
    
    print(out_dict.keys())

    print('Writing json to: ' + outfile)
    with open(outfile, 'wb') as handle:
        pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

