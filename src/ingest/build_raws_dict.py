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

sys_cfg = Dict(json.load(open('etc/conf.json')))

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def build_raws(tstart_str, tend_str, stid, datapath, fmt = "%Y%m%d%H%M"):
    # Given times, stid, and path, return a dictionary of RAWS data fields (fm10, time, lat/lon, elevation) 
    
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
    ind2=np.where(df['d0'].eq(df['d0'][df['d0'] >= tend].min()))[0][0] + 1
    
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
    
    # Get Station data
    f = open(osp.join(datapath, 'stations.pkl'), 'br')
    st = pickle.load(f)
    st = st[st.index == stid]
    
    dict1 = {
        'time': raws1['datetime'].to_numpy(),
        'STID' : raws1['STID'].unique()[0],
        'fm' : raws1['fm10'].to_numpy(),
        'title' : 'RAWS Station '+stid,
        'hours':len(raws1['datetime']),
        'lat' : st['LATITUDE'].to_numpy()[0],
        'lon' : st['LONGITUDE'].to_numpy()[0],
        'other': {'elev': st['ELEVATION']}
    }
    return dict1

def build_hrrr(tstart_str, tend_str, lon, lat):
    # NOT IMPLEMENTED YET, HERE FOR FORMATTING

    hrrr1 = {
        'rain' : 'RAIN',
        'rh' : 'RH',
        'temp' : 'TEMP',
        'Ed' : 'ED',
        'Ew' : 'EW',
        'wind_speed' : 'WIND',
        'solar' : 'SOLAR',
    }
    
    return hrrr1

def build_dictionary(tstart_str, tend_str, stid, datapath, dict_name = "test_dict", fmt = "%Y%m%d%H%M"):
    dict1 = build_raws(tstart_str, tend_str, stid, datapath, fmt)
    
    hrrr1 = build_hrrr(tstart_str, tend_str, dict1['lon'], dict1['lat'])
    
    dict1['id'] = dict_name
    dict1['rain'] = hrrr1['rain'],
    dict1['rh'] = hrrr1['rh'],
    dict1['temp'] = hrrr1['temp'],
    dict1['Ed'] = hrrr1['Ed'],
    dict1['Ew'] = hrrr1['Ew'],
    dict1['wind'] = hrrr1['wind_speed'],
    dict1['solar'] = hrrr1['solar']

    return dict1

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 7:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <STID> <datapath> <target_directory> <filename>' % sys.argv[0]))
        sys.exit(-1)
    start = sys.argv[1]
    end = sys.argv[2]
    stid = sys.argv[3]
    datapath = sys.argv[4]
    outpath = sys.argv[5]

    print('Building RAWS Training Dictionary')
    print('Start: ' + str(start))
    print('End: ' + str(end))

    dict1 = build_dictionary(start, end, stid, datapath)
    
    print(dict1)

    print('Writing json to: ' + outpath)
    print('to_json('+sys.argv[6]+')')

