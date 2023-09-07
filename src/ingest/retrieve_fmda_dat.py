# Set of functions and executable code for retrieving a range of data used for FMDA. 
# Given a start date, end date, and model (currently only tested with HRRR), loop through the dates abnd
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

# Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sys_cfg = Dict(json.load(open('etc/conf.json')))

geotiff_func_path = "/home/hirschij/github/notebooks/fmda/data/grib_to_geotiff.py" # path to grib_to_geotiff.py, used as subprocess

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

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def slice_gribs(grib_file, outpath, grib_source_name):
    
    if grib_source_name == "HRRR": band_df = band_df_hrrr

    filename = osp.join(outpath, Path(osp.basename(grib_path)).stem)
    command = "python " + geotiff_func_path + " " + grib_file + " " + filename

    for i in range(0, band_df.shape[0]):
        band = int(band_df['Band'][i])
        subprocess.call(command + " " + str(band),shell=True)
    return

# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print(('Usage: %s <grib_source_name> <esmf_from_utc> <esmf_to_utc> <target_directory>' % sys.argv[0]))
        print(('Example: %s HRRR 2023-08-10_09:00:00 2023-08-10_10:00:00 ./ingest/HRRR' % sys.argv[0]))
        sys.exit(-1)

    fmt = "%Y-%m-%d_%H:%M:%S" # Time format that pandas can recognize
    

    grib_source = sys.argv[1]
    print('Gathering '+str(grib_source)+' data')
    start=datetime.strptime(sys.argv[2], fmt)
    end=datetime.strptime(sys.argv[3], fmt)
    print("Start: "+str(start))
    print("End: "+str(end))
    outpath=str(sys.argv[4])
    print('Output Destination: '+outpath)


    # String to run with retrieve_gribs as subprocess
    base_str = "python src/ingest/retrieve_gribs.py " + str(sys.argv[1]) + " "

    # Handle Date
    dates = pd.date_range(start=start,end=end, freq="1H") # Series of dates in 1 hour increments

    print('Number of hours:'+str(dates.shape[0]))

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
        
        temppath = osp.join(outpath, dates[t].strftime("%Y%m%d"))
        os.makedirs(temppath, exist_ok=True)

        # Slice grib file to get needed bands
        grib_path = osp.join(sys_cfg.sys_install_path, "ingest", grib_source, grib_source.lower()+'.'+dates[t].strftime("%Y%m%d"),'conus', grib_source.lower()+'.t'+dates[t].strftime('%H')+'z.wrfprsf00.grib2')
        print("Grib Path: " + grib_path)
        slice_gribs(grib_path, temppath, grib_source)
        
        # Remove gribs file for space cleanup
        ## NOTE: leaving associated .size there as paper trail 
        if os.path.exists(grib_path):
            os.remove(grib_path) # remove grib2 file
            grib_path2 = grib_path.replace("wrfprsf00", "wrfprsf01") # get file name for 1 hr forecast
            os.remove(grib_path2) # remove 1 hr forecast file
        else:
            print("The file "+ grib_path +" does not exist")





