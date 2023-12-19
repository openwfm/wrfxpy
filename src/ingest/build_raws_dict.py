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
from MesoPy import Meso

sys_cfg = Dict(json.load(open('etc/conf.json')))
tokens = Dict(json.load(open('etc/tokens.json')))

meso_token = tokens.mesowest
m=Meso(meso_token)

vars='air_temp,relative_humidity,precip_accum,fuel_moisture,wind_speed,solar_radiation' # Variables needed to run FMDA

# Data Retrieval Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def check_data_availability(time_start, state):
    # Query with MesoPy for 1 hour of data to collect STIDs for RAWS that have complete data
    # Currently requires a state string to limit the search, but this could be removed
    # Inputs: 
    # time_start: (str) string of time to start query
    # state: (str) state string, eg "CO"
    # Returns:
    # DataFrame with STIDs and indicators of data availability (whether variable is present at all)
    
    time_s2    = f"{time_start[0:-3]}{int(time_start[-3])+1}{time_start[-2:len(time_start)]}"  # small time increment used to get station ids
    
    print(f"Checking data availability for {state} from {time_start} to {time_s2}")
    # Get one hour of data
    meso_obss = m.timeseries(start=time_start,end=time_s2, state=state,
                                 showemptystations = '0', vars=vars)
    # Set up DF to view data availability
    station_df = pd.DataFrame(columns=['STID', 'air_temp', 'relative_humidity', 'precip_accum', 'fuel_moisture', 'wind_speed', 'solar_radiation'], index=range(0, len(meso_obss["STATION"])))
    # Loop through stations in returned data and add indicator of whether variable is present
    for i in range(0, station_df.shape[0]):
        station_df["STID"][i] = meso_obss["STATION"][i]["STID"]
        station_df["air_temp"][i] = int("air_temp" in meso_obss["STATION"][i]["SENSOR_VARIABLES"].keys())
        station_df["relative_humidity"][i] = int("relative_humidity" in meso_obss["STATION"][i]["SENSOR_VARIABLES"].keys())
        station_df["precip_accum"][i] = int("precip_accum" in meso_obss["STATION"][i]["SENSOR_VARIABLES"].keys())
        station_df["fuel_moisture"][i] = int("fuel_moisture" in meso_obss["STATION"][i]["SENSOR_VARIABLES"].keys())
        station_df["wind_speed"][i] = int("wind_speed" in meso_obss["STATION"][i]["SENSOR_VARIABLES"].keys())
        station_df["solar_radiation"][i] = int("solar_radiation" in meso_obss["STATION"][i]["SENSOR_VARIABLES"].keys())
    return station_df

def filter_stations(df):
    # Given dataframe of RAWS STIDs and data availability indicators, return STIDs of filters
    # Inputs: 
    # df: (DataFrame) output of check_data_availability
    # Return: list of strings, STIDs
    
    print("Filtering to Complete observations for all FMDA variables")
    
    # Filter to stations with complete observations over time period
    df = df[
        (df["fuel_moisture"]==1) & 
        (df["relative_humidity"]==1) &
        (df["precip_accum"]==1) &
        (df["air_temp"]==1) &
        (df["wind_speed"]==1) &
        (df["solar_radiation"]==1)
    ]
    # Extract station IDs
    ids = df['STID'].tolist()

    return ids


def format_raws(stn, fixnames = True):
    # Clean and format FMDA data dictionary 
    # Inputs:
    # stn: Station level dictionary, the part of the output of MesoPy timeseries
    # Outputs: formatted dictionary
    raws_dat = stn['OBSERVATIONS'].copy() # bug fix for in-place changing of dictionary outside of func call
    
    # Convert to Numpy arrays, check data type for floats
    for key in [*stn['OBSERVATIONS'].keys()]:
        if type(stn['OBSERVATIONS'][key][0]) is float:
            raws_dat[key] = np.array(stn['OBSERVATIONS'][key], dtype = 'float64')
        else:
            raws_dat[key] = np.array(stn['OBSERVATIONS'][key])
    
    # Transform Data
    raws_dat['air_temp_set_1'] = raws_dat['air_temp_set_1'] + 273.15 ## convert C to K
    if 'precip_accum_set_1' in raws_dat.keys():
        raws_dat['rain'] = format_precip(raws_dat['precip_accum_set_1']) ## format precip data, accumulated to hourly
    
    
    # Calculate Equilibrium Temps
    raws_dat['Ed'] = 0.924*raws_dat['relative_humidity_set_1']**0.679 + 0.000499*np.exp(0.1*raws_dat['relative_humidity_set_1']) + 0.18*(21.1 + 273.15 - raws_dat['air_temp_set_1'])*(1 - np.exp(-0.115*raws_dat['relative_humidity_set_1']))
    raws_dat['Ew'] = 0.618*raws_dat['relative_humidity_set_1']**0.753 + 0.000454*np.exp(0.1*raws_dat['relative_humidity_set_1']) + 0.18*(21.1 + 273.15 - raws_dat['air_temp_set_1'])*(1 - np.exp(-0.115*raws_dat['relative_humidity_set_1']))
    
    # Fix nan values
    for key in [*raws_dat.keys()]:
        if type(raws_dat[key][0]) is float:
            raws_dat[key] = fixnan(raws_dat[key], 2)
    
    # Add station id
    raws_dat['STID'] = stn['STID']
    
    # Add lat/lon
    raws_dat['LATITUDE'] = stn['LATITUDE']
    raws_dat['LONGITUDE'] = stn['LONGITUDE']
    
    # Simplify names 
    if fixnames:
        var_mapping = {
                'date_time': 'time', 'precip_accum': 'precipa', 'rain':'rain', 'solar_radiation': 'solar',
            'fuel_moisture': 'fm', 'relative_humidity': 'rh',
            'air_temp': 'temp', 'Ed': 'Ed', 'Ew': 'Ew', 'STID': 'STID',
            'LONGITUDE': 'lon', 'LATITUDE': 'lat'
            }
        old_keys = [*raws_dat.keys()]
        old_keys = [k.replace("_set_1", "") for k in old_keys]
        new_keys = []
        for key in old_keys:
            new_keys.append(var_mapping.get(key, key))
        old_keys = [*raws_dat.keys()]
        old_keys = [k.replace("_set_1", "") for k in old_keys]
        new_keys = []
        for key in old_keys:
            new_keys.append(var_mapping.get(key, key))
        raws_dat2 = dict(zip(new_keys, list(raws_dat.values())))
        return raws_dat2
    
    else: return raws_dat

def format_precip(precipa):
    rain=np.array(precipa, dtype = 'float64')
    rain = np.diff(rain) # first difference to convert accumulated to hourly
    rain = np.insert(rain, 0, [np.NaN]) # add NaN entry to account for diff
    # Highest ever recorded hourly rainfall in inches is about 16: https://www.weather.gov/owp/hdsc_world_record
    rain[rain > 100] = np.NaN # filter out erroneously high
    rain[rain < 0] = np.NaN # filter out negative, results from diff function after precipa goes to zero
    return rain

def build_dict(time_start, time_end, stids):
    # Given dates and station IDs, collect RAWS data
    # WARNING: MesoPy will only return data within the last year
    # Inputs: 
    # time_start: (str) format assumed to be %Y%m%d%H%M
    # time_end: (str)
    # stids: list of station ids, output of filter_stations()
    # Returns: formatted FMDA dictionary
    

    print('Gathering data from '+str(time_start)+' to '+str(time_end))
    meso_ts = m.timeseries(time_start, time_end, stid=stids, showemptystations = '0', vars=vars)   # ask the object for data
    # Dictionary to be saved for testing
    raws_dict = {}
    for i in range(0, len(meso_ts['STATION'])):
        raws1 = format_raws(meso_ts['STATION'][i])
        raws1['title']=f"{raws1['STID']}_{time_start}"
        raws1['descr']=f"Data collected from {time_start} to {time_end}. All data is RAWS sensor level data."
        # Filter out if less than 28 days of data, or greater than 20% None observations
        if len(raws1['fm']) < int(24*28) or np.mean(raws1['fm']==None)>.5:
            print(f"Excluding {raws1['STID']}, nobs = {len(raws1['fm'])}, nNone = {np.mean(raws1['fm'] == None)}")
        else:
            raws_dict[raws1['STID']+"_"+time_start] = raws1 # save to test dictionary
    print('Number of Stations: '+str(len(raws_dict)))
    return raws_dict

    

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print(f"Input length must be 5, received {len(sys.argv)}")
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <STATE_CODE> <target_file>' % sys.argv[0]))
        print("Example: python src/ingest/build_raws_dict.py 202305010000 202305302300 CO ~/data/rnn_data_dicts/raws_CO_202305.pickle")
        sys.exit(-1)
    start = sys.argv[1]
    end = sys.argv[2]
    state = sys.argv[3]
    target_file = sys.argv[4]

    df = check_data_availability(start, state) 
    stids = filter_stations(df)
    print(f"N STIDS: {len(stids)}")
    print(stids)

    raws1 = build_dict(start, end, stids)
    
    print(f"Writing output to {target_file}")

    # Write Output
    with open(target_file, 'wb') as handle:
        pickle.dump(raws1, handle, protocol=pickle.HIGHEST_PROTOCOL)


