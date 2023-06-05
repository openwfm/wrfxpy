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

import numpy as np
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
      self.data = pd.DataFrame()
      self.new = False
   #add a pandas dataframe subset from csv
   #def add_process_date(self,df):
   #   self.process_date = 
   def add_incident_id_string(self,id_string):
      self.incident_id_string = id_string
   def add_detections(self,df):
      self.data = self.data.append(df)
   def add_incident_start_time(self):
      if self.data.shape[0] > 0:
         self.incident_start_time = self.data.iloc[0,3]
      else:
         self.incident_start_time = 'Unset start time'
   def new_incident(self):
      self.new = True
   #def add_ignition(self):
   
 #def autostart(job_json):
   #
   
def get_ngfs_csv(days_previous,goes):
  #downloads NGFS csv files from 0 to 4 days previous to now
  #   to download today's csv: get_ngfs_csv(0)
  #example csv name:
  #   NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv
   
  ngfs_dir = 'ingest/NGFS'
  print('    Will download to: ',utils.make_dir(ngfs_dir))
   
  #name the csv file
  csv_day = datetime.now() - timedelta(days = days_previous)
  day_of_year = csv_day.timetuple().tm_yday
  yyyy = csv_day.timetuple().tm_year
  mm = csv_day.timetuple().tm_mon
  dd = csv_day.timetuple().tm_mday
  
  
  
  if goes == 18:
   csv_str = 'NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year)+'.csv'
   csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-WEST/CONUS/'+csv_str
  else:
   csv_str = 'NGFS_FIRE_DETECTIONS_GOES-16_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year)+'.csv'
   csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-EAST/CONUS/'+csv_str
  
  print('    Downloading: ',csv_str)
  csv_path = ngfs_dir+'/'+csv_str
  #downloading
  download_url(csv_url,csv_path)
  return csv_str, csv_path
   
   
def make_job_configuration():
   #check to see if a  base ngfs configuration file exists
   if utils.file_exists('jobs/base_ngfs_cfg.json'):
     print('Configuration found')
     with open('jobs/base_ngfs_cfg.json','r') as openfile:
        temp_cfg = json.load(openfile)
     json_print = json.dumps(temp_cfg, indent=4)   
     #print(json_print)
     start_utc = utils.esmf_to_utc(temp_cfg['start_utc']) #esmf_to_utc
     end_utc = utils.esmf_to_utc(temp_cfg['end_utc'])
     forecast_length = end_utc-start_utc
     print('   Forecast length: ',forecast_length)
     print('   Ignition duration: ',temp_cfg['ignitions']['1'][0]['duration_s'], 'seconds')
     print('   Cell size: ', temp_cfg['domains']['1']['cell_size'])
     print('   Domain size: ', temp_cfg['domains']['1']['domain_size'])
     print('   Number of cores: ',temp_cfg['ppn'])
     print('   Number of nodes: ',temp_cfg['num_nodes'])
     print('   Wall time: ',temp_cfg['wall_time_hrs'], 'hours')
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
   return base_cfg

def read_csv_data(csv_file):
   #parameter csv_file can be a string with a path or a list of strings
   #reads csv file(s) and 
   
   #columns with critical time information
   time_cols = ['incident_start_time','observation_time','initial_observation_time']
   if type(csv_file) is list:
      #read the first csv file in the list
      data = pd.read_csv(csv_file[0], parse_dates=time_cols)
      print('Reading: ',csv_file[0])
      print('  Number of detections: ',data.shape[0])
      #read the rest csv files and merge
      for i in range(1,len(csv_file)):
         print('Merging with csv file: ', csv_file[i])
         data_read = pd.read_csv(csv_file[i], parse_dates=time_cols)
         data = pd.merge(data,data_read, how = 'outer')
         print('  Total number of detections: ',data.shape[0])
      #for naming purposes, get only first file name in list
      #assumes standard  filename for all detections in day
      csv_date_str = csv_file[0][-18:-8]
   else:
      data = pd.read_csv(csv_file, parse_dates=time_cols)
      csv_date_str = csv_file[-18:-8]
   
   #make sure the times are UTC
   data[time_cols[0]] = pd.DatetimeIndex(pd.to_datetime(data[time_cols[0]])).tz_localize('UTC')
   
   return data, csv_date_str
   
      

if __name__ == '__main__':

   print()
   print('Starting ngfs script')
   
   #print(str(sys.argv))
   
   #setup of automatic download of csv file
   #stores csv file in the ingest/NGFS directory
   if 'now' in str(sys.argv):
      
      print('Downloading latest csv files:')
      
      #configure download path, create directory if needed
      # this should be set elsewhere
      ngfs_dir = 'ingest/NGFS'
      
      
      #empy lists for storing path names
      csv_str = list()
      csv_path = list()
      
      #number of days of data
      days_to_get = 2
      
      
      for i in range(days_to_get):
         if i == 0:
            print('NGFS data from today')
         else:
            print('NGFS data from yesterday')
         #GOES-18
         s1, s2 = get_ngfs_csv(i,18)
         csv_str.append(s1)
         csv_path.append(s2)
         #GOES-16
         s1, s2 = get_ngfs_csv(i,16)
         csv_str.append(s1)
         csv_path.append(s2)
      
         
      
      '''
      #current data from GOES 18
      csv_str, csv_path = get_ngfs_csv(0,18)
      #yesterday's data from GOES 18
      csv_str_0, csv_path_0 = get_ngfs_csv(1,18)
      #current data from GOES 16
      csv_str_1, csv_path_1 = get_ngfs_csv(0,16)
      #yesterday's data from GOES 16
      csv_str_2, csv_path_2 = get_ngfs_csv(1,16)
     
      
      #merge two days of csv files
      with open(csv_path,'r') as fp:
         #cut headers from newest csv
         c = fp.readlines()[1:]
      with open(csv_path_0) as fp:
         c_0 = fp.readlines()[0:]
      with open(csv_path_1) as fp:
         c_1 = fp.readlines()[1:]
      with open(csv_path_2) as fp:
         c_2 = fp.readlines()[1:]
      #c_0 += '\n' #linebreak
      c_0 += c_1
      c_0 += c_2
      c_0 += c    # join the two csv files
      
      
      #read first csv_file, include headers
      with open(csv_path[0],'r') as fp:
         c_head = fp.readlines()[0:]
      
      #merge the other csv files
      for i in range(len(csv_path)-1):
         print('Merging ',csv_str[i+1])
         with open(csv_path[i+1],'r') as fp:
            c_new = fp.readlines()[1:]
         c_head += c_new
      
      '''
      
      #name of the new file
      #csv_file = ngfs_dir+'/merge_'+csv_str
      #csv_file = ngfs_dir+'/merged_'+csv_str[0]
      
      csv_file = csv_path
      #print('Will work with:')
      #print(csv_file)
      
      #write the results to the csv file
      #with open(csv_file,'w') as fp:
      #   fp.writelines(c_0[0:])
      #with open(csv_file,'w') as fp:
      #   fp.writelines(c_head[0:])
   else:
      # csv file passed as system argument
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
   
   ##################################################################
   # configure the job(s), uses questionaire from simple_forecast.p #
   ##################################################################
   base_cfg = make_job_configuration()
   
   ##################################################################
   # read the data from the csv file(s), get date string for naming #
   ##################################################################
   data, csv_date_str = read_csv_data(csv_file)
  
   
   
   
   #print('Detection data from ',csv_date_str)
   csv_year = int(csv_date_str[0:4])
   csv_month = int(csv_date_str[5:7])
   csv_day = int(csv_date_str[8:10])
   csv_timestamp = pd.Timestamp(year=csv_year,month=csv_month,day=csv_day)
   #maybe put information about csv file into a class?

   #filter out "null incidents"
   #maybe null incidents are what we are interested in?
   data_dirty = data
   data = data.dropna(subset=['incident_name'])
   #print(data)

   #print('Column headers from csv file')
   #for i in data.columns:
   #   print(i)

   #print('Incidents in csv file')
   incident_names = data['incident_name'].unique()
   num_incidents = incident_names.shape[0]
   
   #print(incident_names)

   #make an array of incident objects
   incidents  = [ngfs_incident(incident_names[i]) for i in range(num_incidents)]

   #list of ignition points for plotting later
   ign_lon = []
   ign_lat = []
   
   #counter
   new_incidents = 0

   for i in range(num_incidents):
      print(incidents[i].name)
      #get only subset of detections from individual nifc incidents
      incident_subset = data[data.incident_name == incident_names[i]]
      #sort by observation time	
      incident_subset = incident_subset.sort_values(by='observation_time')
      det_count = incident_subset.shape[0]
      print('      ',det_count,'  detections') 
      incidents[i].add_detections(incident_subset)
      #print(incidents[i].data)
      incidents[i].add_incident_start_time()
      incidents[i].add_incident_id_string(incident_subset['incident_id_string'][incident_subset.index[0]])
      print('       Incident start time : ', incidents[i].incident_start_time)
      #print('       Incident start day : ', incidents[i].incident_start_time.day)
      #print('       Incident start month: ',incidents[i].incisdent_start_time.month)
      
      #detect a new incident from today
      if (incidents[i].incident_start_time.day == csv_day and incidents[i].incident_start_time.month == csv_month):
      
         new_incidents += 1
         #change 'new' object attribute
         incidents[i].new_incident()
         print('       New Incident From Today')
         print('       Finding start time and location')
         
         idx = incident_subset.index[0] #gets the first row, having earliest time
         ign_latlon = incident_subset.lat_tc[idx], incident_subset.lon_tc[idx]
         unique_latlon = np.unique(ign_latlon)
# --> is first row really the earliest detection ??
         #search the entire csv for locations that have not been named yet
         dirty_subset = data_dirty[data_dirty['lon_tc'] == ign_latlon[1]]
         dirty_subset = dirty_subset[dirty_subset['lat_tc'] == ign_latlon[0]]
         if dirty_subset.shape[0] > incident_subset.shape[0]:
            print('       Found addtional detections')
            ign_utc = dirty_subset.observation_time[dirty_subset.index[0]]
         else:
            ign_utc = incident_subset.observation_time[idx]
         start_utc = utils.round_time_to_hour(ign_utc - timedelta(minutes=30))
         end_utc = start_utc + timedelta(hours=24)
         time_utc = utils.utc_to_esmf(ign_utc)
         print('       location: ',ign_latlon)
         print('       utc ignition time: ', time_utc)
         total_frp = max(incident_subset.total_frp)
         #print('       total frp: ',total_frp)
         #append list of ignition pts for plotting them spatially
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
         json.dump(cfg, open(filename, 'w'), indent=4, separators=(',', ': '))
         print(('INT to start the simulation, execute ./forecast.sh %s' % filename))
   print('Number of incidents in csv file: ',num_incidents)
   print('Number of new incidents: ',new_incidents)
       
   #setup mapping object
   #print('Starting to draw map. Keep fingers crossed')
   #m = Basemap(width=1200000,height=1200000,projection='lcc',resolution='l',lat_0=40.0,lon_0=-121.0)
   #m.latlon = True 
   #x, y = m(ign_lon,ign_lat)
   #m.scatter(x,y,s=15)
   #m.shadedrelief()
   #m.drawstates()
   #plt.show
