## NOTE: replace hrrr_conf json with argv. Using argv fails below

# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict
import json
import os.path as osp
from datetime import date, timedelta, datetime
import pandas as pd
import subprocess
import sys

# Packages added manually to wrfx environment
import xarray as xr

sys_cfg = Dict(json.load(open('etc/conf.json')))

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def slice_hrrr(tempfile):

    ds1=xr.open_dataset(
            tempfile,
            filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 'instant'}
    )
    ds2=xr.open_dataset(
            tempfile,
            filter_by_keys={'typeOfLevel': 'surface', 'stepType': 2}
    )
        
    ds2=ds2.assign_coords({'heightAboveGround': np.float64(0)}) # Add height above ground field
    ds3=xr.open_dataset(
            tempfile,
            filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10}
        )


    return ds1, ds2, ds3

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <target_directory>' % sys.argv[0]))
        print(('Example: %s 2023-08-10_09:00:00 2023-08-10_10:00:00 ./ingest/HRRR' % sys.argv[0]))
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
    dates = pd.date_range(start=start,end=end, freq="1H") # Series of dates in 1 hour increments

    if dates.shape[0]>10: 
            print('Dont run this many times yet, test more')
            sys.exit(-1)

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

        # Slice HRRR grib file into surface, 2m, and 10m respectively
        hrrr_path="/home/hirschij/github/wrfxpy/ingest/HRRR/hrrr.20230113/conus/hrrr.t06z.wrfprsf00.grib2" # path to observed data, NOT THE RIGHT FILENAME CHANGE TO MATCH LOOP
        ds1, ds2, ds3 = slice_hrrr(hrrr_path)






