## NOTE: replace hrrr_conf json with argv. Using argv fails below

# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict
import json
import os.path as osp
# import xarray as xr
from datetime import date, timedelta, datetime
import pandas as pd
import subprocess
import sys

sys_cfg = Dict(json.load(open('etc/conf.json')))
cfg = Dict(json.load(open('wksp/hrrr_conf.json')))

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <forecast_hours> <target_directory>' % sys.argv[0]))
        print(('Example: %s 2023-08-10_09:00:00 2023-08-10_10:00:00 1 ./ingest/HRRR' % sys.argv[0]))
        sys.exit(-1)

    fmt = "%Y-%m-%d_%H:%M:%S" # Time format that pandas can recognize

    print('Gathering HRRR data')
    #starttime = datetime.strptime(sys.argv[1], fmt)
    #endtime = datetime.strptime(sys.argv[3], fmt)
    print('Start:'+str(cfg.start_time))
    print(sys.argv[1])
    print('End:'+str(cfg.end_time))
    print(sys.argv[3])


    # String to run with retrieve_gribs as subprocess
    base_str = "python src/ingest/retrieve_gribs.py HRRR "

    # Handle Date
    starttime = datetime.strptime(cfg.start_time, fmt)
    endtime = datetime.strptime(cfg.end_time, fmt)
    dates = pd.date_range(start=starttime,end=endtime, freq="1H") # Series of dates in 1 hour increments

    for t in range(0, 2): # replace 2 with dates.shape[0]
        print('-'*30)
        
        # Format time
        time=dates[t].strftime(fmt)
        tforecast = dates[t] + pd.DateOffset(hours = 1)
        tforecast = tforecast.strftime(fmt)

        temppath = osp.join(cfg.dest_dir, "temp")
        command = base_str + str(time) + " " + str(tforecast) + " " + temppath
        print(command)
        #subprocess.call(command,shell=True)


    # Temporarily download grib file at given time
    #t0 = datetime.strftime(dates[0], "%Y-%m-%d_%H:%M:%S") # UTC date format used by wrfxpy
    #tforecast = dates[0] + pd.DateOffset(hours = 1) # time to look forward to in forecast (confirm this is what is actually happening in retrieve_gribs)
    #tforecast = datetime.strftime(tforecast, fmt)
    #temppath = osp.join(cfg.dest_dir, "temp")
    #command = base_str + str(t0) + " " + str(tforecast) + " " + temppath
    
    #print(command)
    # subprocess.call(command,shell=True)
