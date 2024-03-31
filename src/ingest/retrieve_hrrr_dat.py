# Set of functions and executable code for retrieving a range of data from HRRR used for FMDA. 
# Given a start date, end date, loop through the dates and
# execute retrieve_gribs from the main code repo, then extract and combine relevant fields to given output dir

# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict
import json
import os
import os.path as osp
from datetime import date, timedelta, datetime
import pandas as pd
import subprocess
import sys
from pathlib import Path
import numpy as np

# Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sys_cfg = Dict(json.load(open('etc/conf.json')))

geotiff_func_path = "/home/hirschij/github/notebooks/fmda/data/grib_to_geotiff.py" # path to grib_to_geotiff.py, used as subprocess

# DataFrame to track needed Bands for FMDA and metadata
## Note there is a surface temp, but no surface RH. 
## Sticking with 2m temp and rh for consistency since those are combined to EQ
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
              'surface Ground Heat Flux [W/m^2]'],
    'notes': "Precipitation data for HRRR is always 0 at time 0"
})

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def slice_gribs(grib_file, outpath):

    filename = osp.join(outpath, Path(osp.basename(grib_path)).stem)
    command = "python " + geotiff_func_path + " " + grib_file + " " + filename

    for i in range(0, band_df_hrrr.shape[0]):
        band = int(band_df_hrrr['Band'][i])
        subprocess.call(command + " " + str(band),shell=True)
    return

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

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print(('Usage: %s <esmf_from_utc> <esmf_to_utc> <forecast_hours> <target_directory>' % sys.argv[0]))
        print(('Example: %s 2024-01-01_00:00:00 2024-01-01_01:00:00 1 ./ingest/HRRR' % sys.argv[0]))
        sys.exit(-1)

    fmt = "%Y-%m-%d_%H:%M:%S" # Time format that pandas can recognize
    
    grib_source = "HRRR"
    print('Gathering HRRR data')
    start=datetime.strptime(sys.argv[1], fmt)
    end=datetime.strptime(sys.argv[2], fmt)
    print(f"Start Time: {start}")
    print(f"End Time: {end}")
    forecast_hours = sys.argv[3]
    print(f"HRRR Forecast Hours: {forecast_hours}")
    outpath=str(sys.argv[4])
    print('Output Destination: '+outpath)


    # String to run with retrieve_gribs as subprocess
    base_str = f"python src/ingest/retrieve_gribs.py {grib_source} "

    # Handle Date
    dates = pd.date_range(start=start,end=end, freq="1H") # Series of dates in 1 hour increments
    # Forecast hour format
    ft = ftimes(forecast_hours)

    print(f'Number of analysis hours: {dates.shape[0]}')

    for t in range(0, dates.shape[0]): 
        print('~'*30)
        
        # Format time
        time=dates[t].strftime(fmt) # start time for retrieve_gribs
        tforecast = dates[t] + pd.DateOffset(hours = int(forecast_hours)) # end time for retrieve_gribs
        tforecast = tforecast.strftime(fmt)
        
        # Call retrieve_gribs via subprocess
        command = base_str + str(time) + " " + str(tforecast) + " ~"
        print(command)
        subprocess.call(command,shell=True)
        # Set up destination file path
        temppath = osp.join(outpath, dates[t].strftime("%Y%m%d"))
        os.makedirs(temppath, exist_ok=True)

        # Loop over analysis/forecast hours and slice_gribs to get needed HRRR bands
        ft = ftimes(forecast_hours)
        for ts in ft:
            grib_path = osp.join(sys_cfg.sys_install_path, "ingest", grib_source, grib_source.lower()+'.'+dates[t].strftime("%Y%m%d"),'conus', f"{grib_source.lower()}.t{dates[t].strftime('%H')}z.wrfprsf{ts}.grib2")
            print("Grib Path: " + grib_path)
            if os.path.exists(grib_path):
                slice_gribs(grib_path, temppath)
            else:
                print(f"File does not exist: {grib_path}")

            # Remove gribs file for space cleanup
            ## NOTE: leaving associated .size there as paper trail 
            if os.path.exists(grib_path):
                os.remove(grib_path) # remove grib2 file

