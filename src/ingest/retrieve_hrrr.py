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
#cfg = Dict(json.load(open('wksp/hrrr_conf.json')))

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <target_directory>' % sys.argv[0]))
        print(('Example: %s 2023-08-10_09:00:00 2023-08-10_10:00:00 1 ./ingest/HRRR' % sys.argv[0]))
        sys.exit(-1)

    fmt = "%Y-%m-%d_%H:%M:%S" # Time format that pandas can recognize

    print('Gathering HRRR data')
    start=datetime.strptime(sys.argv[1], fmt)
    end=datetime.strptime(sys.argv[2], fmt)
    print("Start: "+str(start))
    print("End: "+str(end))
    outpath=str(sys.argv[3])
    print('Output Destination: '+outpath)


    # String to run with retrieve_gribs as subprocess
    base_str = "python src/ingest/retrieve_gribs.py HRRR "

    # Handle Date
    #starttime = datetime.strptime(cfg.start_time, fmt)
    #endtime = datetime.strptime(cfg.end_time, fmt)
    dates = pd.date_range(start=start,end=end, freq="1H") # Series of dates in 1 hour increments

    for t in range(0, dates.shape[0]): 
        print('-'*30)
        
        # Format time
        time=dates[t].strftime(fmt)
        tforecast = dates[t] + pd.DateOffset(hours = 1)
        tforecast = tforecast.strftime(fmt)

        #temppath = osp.join(cfg.dest_dir, "temp")
        
        command = base_str + str(time) + " " + str(tforecast) + " ~"
        print(command)
        subprocess.call(command,shell=True)


