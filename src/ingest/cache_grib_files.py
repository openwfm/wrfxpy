from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd
import pickle
import copy
import csv
import os, sys, glob
import json
import time
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import timedelta, datetime
#add src directory to path json time
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from ingest.downloader import download_url
#from ngfs_dictionary import ngfs_dictionary  <-- use to force csv data columns into a type
import ngfs_helper as nh
import utils
#dictionary to convert between "CA" and "California", etc
import state_names as sn
import simple_forecast as sf
import ngfs_dictionary as nd
#import shapely
from shapely.geometry import Point, LineString, Polygon
try:
   from pyproj import Proj, transform, Transformer
except:
   from pyproj import Proj, transform # Transformer

print('Starting ingest code. Will try to cache all necessary grib files for forecasts')
current_utc = pd.Timestamp.now(tz='UTC')
print('\tCurrent UTC time: ',utils.utc_to_esmf(current_utc))
print('       supported GRIB sources: HRRR, HRRR_AK, NAM, NAM227, NAM196, NAM198, CFSR_P, CFSR_S, NARR, GFSA, GFSF_P, GFSF_S')
#only cache ALASKA during May to September
if current_utc.month in [5,6,7]:
    sats = ['NAM','NAM198','HRRR']
else:
    sats = ['NAM','HRRR']

#Sources whose cycle must be pinned explicitly.  NAM218 cycles every 6 hours, so
#letting the source choose "the most recent cycle" lands on 00/06/12/18 by itself and
#the cycle_start this loop computes is implicit in the time range.  **HRRR cycles every
#hour** (HRRR.cycle_hours = 1), so the same call would quietly pick whatever hourly
#cycle is newest -- 14z at 15:40Z -- and the 00/06/12/18 scheme this loop exists to
#implement would be silently lost, giving a cache of overlapping cycles that no longer
#lines up with the FMDA forecast cycles it is meant to share.
PINNED_CYCLE = {'HRRR'}

retrieve_cmd = './retrieve_gribs.sh {} {} {} ingest'
lookback_cycles = 2
for s in sats:
    print('\tWill download grib files for ',s)
    #loop through previous time steps to ensure cache cycle is complete
    for l in range(lookback_cycles):
        now = current_utc - timedelta(hours=6*l)
        #compute the cycle parameters
        cycle_hour = np.int8(np.int8(now.hour/6))*6
        cycle_start = pd.Timestamp(now.year,now.month,now.day,cycle_hour)
        print('\tStarting times: ',cycle_start)
        cycle_end = cycle_start + timedelta(hours = 33)
        print('\tEnding  times: ',cycle_end)
        cmd_str = retrieve_cmd.format(s,utils.utc_to_esmf(cycle_start),utils.utc_to_esmf(cycle_end))
        if s in PINNED_CYCLE:
            cmd_str += ' ' + utils.utc_to_esmf(cycle_start)
        print('\t',cmd_str)
        os.system(cmd_str)


