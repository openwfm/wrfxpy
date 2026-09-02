from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
#import geopandas as gpd
import pandas as pd
#import pickle
#import copy
import csv

from ngfs import config_manager
from ngfs.ngfs_day import ngfs_day
from ngfs.persistence import get_old_incidents

try:
   from pyproj import Proj, transform, Transformer
except:
   from pyproj import Proj, transform # Transformer
#from ingest.NWS_warnings import set_red_flags, red_flag_incident, subset_red_flag_data, subset_wfo_data

####### Functions  #######

#####
if __name__ == "__main__":
   print()
   now_local = pd.Timestamp.now()
   print(now_local)
   now_utc = pd.Timestamp.now('UTC')
   print('Starting ngfs script ',now_local)
   ### load or make configurations  ###load_cfgs()
   ngfs_cfg, wrfxpy_cfg = config_manager.load_cfgs()
   ngfs_directory = ngfs_cfg['ngfs_directory']

   #determine the number of forecastrs to possibly run and whether these will start automatically
   #auto_start, num_starts, force = config_manager.setup_auto_start(sys.argv,ngfs_cfg)  <<-------------------------- Remove?
   
   ### load previous incident information  ###
   started_inc_ids, old_incidents, old_ngfs_day  = get_old_incidents(ngfs_directory)
         
   #make ngfs_day object
   csv = ngfs_day(ngfs_cfg,started_inc_ids=started_inc_ids,old_incidents=old_incidents)

   #overirde any config settings with command line options
   csv.sys_args_override()

   #add detection and other data
   print('Adding data')
   csv.add_goes_data(df = old_ngfs_day.data.drop_duplicates() if len(old_ngfs_day.data) > 0 else None)
   #print(csv.data.keys())
   #add polar data for today forecasts
   #csv.add_polar_data()    <<-------------------------- Remove, this is being handled with add_viirs_data

   # add VIIRS data by new functions, if it works, remove that above
   csv.add_viirs_data()      ### maybe make this not return anything? 
   #get red flag warning data
   csv.add_red_flags()
   #adding population data
   csv.add_pop_data()

   #add the incidents from the data
   print('Adding incidents')
   csv.add_incidents()

   #process the incidents
   csv.process_incidents()

   #start forecasts
   csv.start_incidents()
   #csv.prioritize_incidents()

   #print map, save pickle file etc
   csv.save_outputs()

   print(pd.Timestamp.now()-now_local)


   

   


      
