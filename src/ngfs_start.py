# reads ngfs csv file and lists incidents and starts a simulation from user selection
# borrows heavily from simple_forecast.py
# usage for existing csv file:
#     python ngfs_start.py /path/to/ngsf_csv_file.csv
# usage for getting data from today:
#     python ngfs_start.py now
# if base configuration has already been created, then all job scripts may be
#    created without user input with:
#         python ngfs_start.py /path/to/ngsf_csv_file.csv auto
#         ngfs_start.py now auto
#        
# needs to be run under the correct python environment for wrfxpy --> conda activate  wrfx

#to do
#  1) get csv files automatically from:
#        https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-WEST/CONUS/NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv
#        /pub/volcat/fire_csv/NGFS_daily/GOES-WEST/CONUS
#        NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv


import pandas as pd
import csv
import sys
import json
import time
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import timedelta, datetime
#add src directory to path
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from downloader import download_url
import utils
import simple_forecast as sf

class ngfs_incident():
   def __init__(self,name):
      self.name = name
   #add a pandas dataframe subset from csv
   #def add_process_date(self,df):
   #   self.process_date = 
   def add_incident_id_string(self,id_string):
      self.incident_id_string = id_string
   def add_detections(self,df):
      self.data = df
   def add_incident_start_time(self):
      if self.data.shape[0] > 0:
         self.incident_start_time = self.data.iloc[0,3]
      else:
         self.incident_start_time = 'Unset start time'
   #def add_ignition(self):
      

if __name__ == '__main__':

   print()
   print('Starting ngfs script')
   
   #print(str(sys.argv))
   
   #setup of automatic download of csv file
   #stores csv file in the ingest/NGFS directory
   if 'now' in str(sys.argv):
      print('Attempting to download latest csv file:')
      #example csv name:
      #   NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv
      day_of_year = datetime.now().timetuple().tm_yday
      yyyy = datetime.now().timetuple().tm_year
      mm = datetime.now().timetuple().tm_mon
      dd = datetime.now().timetuple().tm_mday
      csv_str = 'NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year)+'.csv'
      csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-WEST/CONUS/'+csv_str
      print('    ',csv_str)
      #print('    ',csv_url)
      #configure download path, create directory if needed
      ngfs_dir = 'ingest/NGFS'
      print('     Will download to: ',utils.make_dir(ngfs_dir))
      
      csv_path = 'ingest/NGFS/'+csv_str
      download_url(csv_url,csv_path)
      
      csv_file = csv_path
   else:
      csv_file = sys.argv[1]
      print('Using existing csv file:')
      print('    ',csv_file) 
  
   
   
   #setup automatic start of simulations
   if 'auto' in str(sys.argv):
      sf.print_question('Autostart detected. Is this what you want? yes/no, default = [no]')
      auto_start = sf.read_boolean('no')
   else:
      auto_start = False
   if auto_start:
      print('Simulations will be started automatically. Make sure you have the resources.')
      time.sleep(2)
   else:
      print('Simulations will need to be started manually')
   
   #check to see if a  base ngfs configuration file exists
   if utils.file_exists('jobs/base_ngfs_cfg.json'):
      print('Configuration found')
      with open('jobs/base_ngfs_cfg.json','r') as openfile:
         temp_cfg = json.load(openfile)
      json_print = json.dumps(temp_cfg, indent=4)   
      print(json_print)
      sf.print_question('Use base configuration above? yes/no, default = [yes]')
      use_base_cfg = sf.read_boolean('yes')
      if use_base_cfg:
         base_cfg = temp_cfg
      else:
         base_cfg = sf.questionnaire()   
   else:
      # standard configuration to be used by all forecasts
      base_cfg = sf.questionnaire()
      print(base_cfg)
   # save base confg for next time
   json.dump(base_cfg, open('jobs/base_ngfs_cfg.json', 'w'), indent=4, separators=(',', ': '))
      
   
   #read the data and filter
   #columns with critical time information
   time_cols = ['incident_start_time','observation_time','initial_observation_time']
   data = pd.read_csv(csv_file, parse_dates=time_cols)
   #make sure the times are UTC
   data[time_cols[0]] = pd.DatetimeIndex(pd.to_datetime(data[time_cols[0]])).tz_localize('UTC')
   #assumes standard  filename for all detections in day
   #can be changed to use reg expression so that it's more general
   csv_date_str = csv_file[-18:-8]
   print('Detection data from ',csv_date_str)
   csv_year = int(csv_date_str[0:4])
   csv_month = int(csv_date_str[5:7])
   csv_day = int(csv_date_str[8:10])
   csv_timestamp = pd.Timestamp(year=csv_year,month=csv_month,day=csv_day)
   #maybe put information about csv file into a class?

   #filter out "null incidents"
   #maybe null incidents are what we are interested in?
   data = data.dropna(subset=['incident_name'])
   print(data)

   #print('Column headers from csv file')
   for i in data.columns:
      print(i)

   print('Incidents in csv file')
   incident_names = data['incident_name'].unique()
   num_incidents = incident_names.shape[0]
   print('Number of incidents in csv file: ',num_incidents)
   print(incident_names)


      
   #array of incident objects
   incidents  = [ngfs_incident(incident_names[i]) for i in range(num_incidents)]

   ign_lon = []
   ign_lat = []

   for i in range(num_incidents):
      print(incidents[i].name)
      #get only detections from individual nifc incidents
      incident_subset = data[data.incident_name == incident_names[i]]
      #sort by observation time	
      incident_subset = incident_subset.sort_values(by='observation_time')
      #incident_subset_idx = data[data.incident_name == incident_names[i]].index
      det_count = incident_subset.shape[0]
      print('      ',det_count,'  detections') 
      incidents[i].add_detections(incident_subset)
      #print(incidents[i].data)
      incidents[i].add_incident_start_time()
      incidents[i].add_incident_id_string(incident_subset['incident_id_string'][incident_subset.index[0]])
      print('       Incident start time : ', incidents[i].incident_start_time)
      print('       Incident start day : ', incidents[i].incident_start_time.day)
      #print('       Incident start month: ',incidents[i].incisdent_start_time.month)
      if (incidents[i].incident_start_time.day == csv_day and incidents[i].incident_start_time.month == csv_month):
         print('       New Incident From Today')
         print('       Finding start time and location')
         idx = incident_subset.index[0]
         #switch to terrain corrected locations
         ign_latlon = incident_subset.lat[idx], incident_subset.lon[idx]
         ign_utc = incident_subset.observation_time[idx]
         start_utc = utils.round_time_to_hour(ign_utc - timedelta(minutes=30))
         end_utc = start_utc + timedelta(hours=5)
         time_utc = utils.utc_to_esmf(ign_utc)
         print('       location: ',ign_latlon)
         print('       utc ignition time: ', time_utc)
         total_frp = max(incident_subset.total_frp)
         print('       total frp: ',total_frp)
         #append list of ignition pts
         ign_lon.append(ign_latlon[1])
         ign_lat.append(ign_latlon[0])
         #print(ign_lon)
         
         #change the base configuration file to match the individual incidents parameters
         cfg = base_cfg
         gc_string  = (incident_subset['incident_name'][idx]+'_'+utils.utc_to_esmf(start_utc)+
                                         '_'+incident_subset['incident_id_string'][idx][1:-1])
         cfg['grid_code'] = gc_string.replace(' ','_')
         cfg['domains']['1']['center_latlon'] = ign_latlon
         cfg['domains']['1']['truelats'] = (ign_latlon[0], ign_latlon[0])
         cfg['domains']['1']['stand_lon'] = ign_latlon[1]
         ignitions = cfg['ignitions']
         ign_dur = ignitions['1'][0]['duration_s']
         cfg['ignitions'] = { '1' : [ { 'time_utc' : time_utc,
                                   'duration_s' : ign_dur,
                                   'latlon' : ign_latlon } ] }
         cfg['start_utc'] = utils.utc_to_esmf(start_utc)
         cfg['end_utc'] = utils.utc_to_esmf(end_utc)
         #print(cfg)

         filename = 'jobs/' + cfg['grid_code'] + '.json'
         #json.dump(cfg, open(filename, 'w'), indent=4, separators=(',', ': '))
         print(('INT to start the simulation, execute ./forecast.sh %s' % filename))
       
   #setup mapping object
   #print('Starting to draw map. Keep fingers crossed')
   #m = Basemap(width=1200000,height=1200000,projection='lcc',resolution='l',lat_0=40.0,lon_0=-121.0)
   #m.latlon = True 
   #x, y = m(ign_lon,ign_lat)
   #m.scatter(x,y,s=15)
   #m.shadedrelief()
   #m.drawstates()
   #plt.show
