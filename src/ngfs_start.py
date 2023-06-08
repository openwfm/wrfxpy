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
#add src directory to path json time
#sys.path.insert(1, 'src/')
#sys.path.insert(1, 'src/ingest')
from ingest.downloader import download_url
import ingest.ngfs_helper 
import utils
#dictionary to convert between "CA" and "California", etc
import state_names as sn
import simple_forecast as sf

#To add or change

####### Classes  ######

#class that keeps atrributes about the day's data
class ngfs_day():
   #maybe add some attributes so that ongoing, new, and started incidents can be tracked
   def __init__(self,csv_date_str,today):
      self.date_str = csv_date_str
      self.timestamp = timestamp_from_string(csv_date_str)
      self.today = today  # <--- Boolean as to whether script was called with 'now', maybe change name to 'now'
      self.incidents = list()
      self.incident_names = list()
      #self.new_incident = list()
      #self.started_incident = list()
   def add_incidents(self,incident):
      self.incidents.append(incident)  ## <--- This will hold all the informations
      self.incident_names.append(incident.name)  ## needed?
   #def save_ngfs_day(self):
      #save with pickle
      #make script to run through all saved files to plot


class ngfs_incident():
   def __init__(self,name):
      self.name = name
      self.data = pd.DataFrame()
      self.new = False
      self.affected_population = 0
   def process_incident(self,data):
      #takes in DataFrame and fills in details about incident
      pass
   def new_incident(self):
      self.new = True
   #could be more than one location, use append?
   def add_incident_location(self,county,state):
      if len(state) == 2:
         state = sn.abbrev_to_us_state[state]
      #matches the first column of population data text file
      self.loc_str = '.'+county+', '+state
      self.state = state
      self.county = county
   def set_population(self,pop):
      self.affected_population = int(pop)
   def set_incident_id_string(self,id_string):
      self.incident_id_string = id_string
   #these are the columns of the csv file as pandas dataframe
   def add_detections(self,df):
      self.data = self.data.append(df)
   def set_incident_start_time(self):
      if self.data.shape[0] > 0:
         self.incident_start_time = self.data.iloc[0,3]
      else:
         self.incident_start_time = 'Unset start time'
   def set_json_start_code(self,json_file,grid_code):
      self.json_start_code = './forecast.sh '+json_file+' &> logs/'+grid_code.replace(' ','_')+'.log &'      
   #def add_ignition(self):
   '''
   def get_incident_bounding_box(self):
      #returns boundind box 
      lats = 
      lons = 
      bbox = (np.min(lats), np.min(lons), np.max(lats), np.max(lons))
   '''
   def make_incident_configuration(self,base_cfg):
      #change the base configuration file to match the individual incidents parameters
      #this is not done for older incidents, why? does it save time? makes code harder to read
      cfg = base_cfg
      gc_string  = self.name+'_'+utils.utc_to_esmf(self.start_utc)+'_'+self.incident_id_string[1:-1]
      cfg['grid_code'] = gc_string.replace(' ','_')
      cfg['domains']['1']['center_latlon'] = self.ign_latlon
      cfg['domains']['1']['truelats'] = (self.ign_latlon[0], self.ign_latlon[0])
      cfg['domains']['1']['stand_lon'] = self.ign_latlon[1]
      ignitions = cfg['ignitions']
      ign_dur = ignitions['1'][0]['duration_s']
      cfg['ignitions'] = { '1' : [ { 'time_utc' : self.time_utc,
                                   'duration_s' : ign_dur,
                                   'latlon' : self.ign_latlon } ] }
      cfg['start_utc'] = utils.utc_to_esmf(self.start_utc)
      cfg['end_utc'] = utils.utc_to_esmf(self.end_utc)
         #print(cfg)

      self.filename = 'jobs/' + cfg['grid_code'] + '.json'
      self.set_json_start_code(self.filename,cfg['grid_code']) 
      self.cfg = cfg

   def process_incident(self,incident_subset,data,data_full):

      #sort by observation time	
      incident_subset = incident_subset.sort_values(by='observation_time')

      self.det_count = incident_subset.shape[0]
      print('\t',self.det_count,'detections') 
      
      #demographic data for incident
      #could span more than one county, state --> use unique to generate a list?
      self.county = incident_subset['locale'][incident_subset.index[0]]
      self.state = incident_subset['state'][incident_subset.index[0]]
      if len(self.state) == 2:
         self.state = sn.abbrev_to_us_state[self.state]
      #matches the first column of population data text file
      self.loc_str = '.'+self.county+', '+self.state
      #print('\t',incidents[i].loc_str)
      
      #location and population data for incident
      loc_idx = pop_data[pop_data['Location'] == self.loc_str]
      
      if loc_idx.index.empty:
         self.affected_population = float('NaN')
         print('\tNo matching population data found')
      else:
         pop = loc_idx.iloc[0]['Population'] #this is a string with commas
         self.affected_population = float(pop.replace(',',''))  # convert string to float
         print('\tIncident county population is ',self.affected_population)
      
      #add detections to the incidnet
      self.add_detections(incident_subset)
      self.set_incident_start_time()
      self.set_incident_id_string(incident_subset['incident_id_string'][incident_subset.index[0]])
      print('\tIncident start time : ', self.incident_start_time)

      
      print('\tFinding start time and ignition location')
         
      idx = incident_subset.index[0] #gets the first row, having earliest time
      #not all csv files have terrain correction, so far
      try:
         self.ign_latlon = incident_subset.lat_tc[idx], incident_subset.lon_tc[idx]
      except:
         print('\tNo terrain corected lat/lon available')
         self.ign_latlon = incident_subset.lat[idx], incident_subset.lon[idx]
      #reduce to only unique lat/lon pairs
      self.unique_latlon = np.unique(self.ign_latlon)

      
      #search the entire csv for locations that have not been named yet
      try:
         dirty_subset = data_full[data_full['lon_tc'] == self.ign_latlon[1]]
         dirty_subset = dirty_subset[dirty_subset['lat_tc'] == self.ign_latlon[0]]
      except: #if no terrain corrected data present
         dirty_subset = data_full[data_full['lon'] == self.ign_latlon[1]]
         dirty_subset = dirty_subset[dirty_subset['lat'] == self.ign_latlon[0]]
      if dirty_subset.shape[0] > incident_subset.shape[0]:
         print('\tFound addtional, earlier detections')
         self.ign_utc = dirty_subset.observation_time[dirty_subset.index[0]]
      else:
         self.ign_utc = min(incident_subset.observation_time[idx],self.incident_start_time)
      
         
      
      #date to start the simulation, 30 minutes before the ignition
      self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=30))
      self.end_utc = self.start_utc + timedelta(hours=24)
      #time_utc goes into the namelist file
      self.time_utc = utils.utc_to_esmf(self.ign_utc)
      print('\tlocation: ',self.ign_latlon)
      print('\tutc ignition time: ', self.time_utc)
      #need to learn more about the FTP and how it's made
      self.total_frp = max(incident_subset.total_frp)
      #print('\ttotal frp: ',total_frp)

 #def autostart(job_json):
   #
        
      
####### Functions  #######
   
def make_base_configuration():
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
     print('\tForecast length: ',forecast_length)
     print('\tIgnition duration: ',temp_cfg['ignitions']['1'][0]['duration_s'], 'seconds')
     print('\tCell size: ', temp_cfg['domains']['1']['cell_size'])
     print('\tDomain size: ', temp_cfg['domains']['1']['domain_size'])
     print('\tNumber of cores: ',temp_cfg['ppn'])
     print('\tNumber of nodes: ',temp_cfg['num_nodes'])
     print('\tWall time: ',temp_cfg['wall_time_hrs'], 'hours')
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
   
def timestamp_from_string(csv_date_str):
   #print('Detection data from ',csv_date_str)
   csv_year = int(csv_date_str[0:4])
   csv_month = int(csv_date_str[5:7])
   csv_day = int(csv_date_str[8:10])
   return pd.Timestamp(year=csv_year,month=csv_month,day=csv_day)
   
    
def get_ngfs_csv(days_previous,goes):
  #downloads NGFS csv files from 0 to 4 days previous to now
  #   to download today's csv: get_ngfs_csv(0)
  #example csv name:
  #   NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv
   
  ngfs_dir = 'ingest/NGFS'
  print('\tWill download to: ',utils.make_dir(ngfs_dir))
   
  #name the csv file
  csv_day = datetime.now() - timedelta(days = days_previous)
  day_of_year = csv_day.timetuple().tm_yday
  yyyy = csv_day.timetuple().tm_year
  mm = csv_day.timetuple().tm_mon
  dd = csv_day.timetuple().tm_mday
  
  #join strings to have name and url of the csv file
  if goes == 18:
   csv_str = 'NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year)+'.csv'
   csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-WEST/CONUS/'+csv_str
  else:
   csv_str = 'NGFS_FIRE_DETECTIONS_GOES-16_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year)+'.csv'
   csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-EAST/CONUS/'+csv_str
  
  print('\tDownloading: ',csv_str)
  csv_path = ngfs_dir+'/'+csv_str
  #downloading
  download_url(csv_url,csv_path)
  return csv_str, csv_path
   
def download_csv_data(days_to_get):
   #downloads csv data by ftp, returns list of csv files obtained
   #empy lists for storing path names
   csv_str = list()
   csv_path = list()
      
   #number of days of data
   #days_to_get = 2
       
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
   return csv_str, csv_path
   

def read_NGFS_csv_data(csv_file):
   #parameter csv_file can be a string with a path or a list of strings
   #reads csv file(s) and merges them. Assigns a date to them too
   
   #columns with critical time information
   time_cols = ['incident_start_time','observation_time','initial_observation_time']
   if type(csv_file) is list:
      #read the first csv file in the list
      data = pd.read_csv(csv_file[0], parse_dates=time_cols)
      print('Reading: ',csv_file[0])
      print('\tNumber of detections: ',data.shape[0])
      #read the rest csv files and merge
      for i in range(1,len(csv_file)):
         print('Merging with csv file: ', csv_file[i])
         data_read = pd.read_csv(csv_file[i], parse_dates=time_cols)
         data = pd.merge(data,data_read, how = 'outer')
         print('\tTotal number of detections: ',data.shape[0])
      #for naming purposes, get only first file name in list
      #assumes standard  filename for all detections in day
      csv_date_str = csv_file[0][-18:-8]
   else:
      data = pd.read_csv(csv_file, parse_dates=time_cols)
      csv_date_str = csv_file[-18:-8]
   
   #make sure the times are UTC
   data[time_cols[0]] = pd.DatetimeIndex(pd.to_datetime(data[time_cols[0]])).tz_localize('UTC')
   
   return data, csv_date_str

def incident_demographics(data, pop_data):
   pass
   
####### indexing here is bad
def prioritize_incidents(incidents,new_idx,num_starts):
   #prioritizes the new incidents to be run by county population
   #new_dx is list of the new incidents
   #auto_start is number of the new simulations to start
   n = len(new_idx)
   started = np.zeros(n,dtype=int)
   pop = np.zeros(n)
   for i in range(n):
      pop[i] = incidents[i].affected_population
   #to get decreasing list, sort its negative
   sort_idx = np.argsort(-pop)
   print('New incidents by priority')
   start_count = 0
   for i in sort_idx:
      if incidents[i].new:
         print(incidents[i].name)
         print('\t Population affected: ',incidents[i].affected_population)
         print('\t',incidents[i].county,', ',incidents[i].state)
         #print('\t',incidents[i].json_start_code)
         if start_count < num_starts:
            print('\t Automatically starting this simulation')
            started[i] = 1
            start_count += 1
 
   return sort_idx, started

if __name__ == '__main__':

   print()
   print('Starting ngfs script')
   
   print('Reading in county population data')
   pop_data = pd.read_csv('ingest/NGFS/Population_by_US_County_July_2022.txt',sep='\t',encoding = "ISO-8859-1")
   
   #setup of automatic download of csv file
   #stores csv file(s) in the ingest/NGFS directory
   if 'now' in str(sys.argv):   
      print('Downloading latest csv files:')
      #configure download path, create directory if needed

      # this should be set elsewhere
      ngfs_dir = 'ingest/NGFS'
      
      #download the data
      days_to_get = 2
      csv_str, csv_path = download_csv_data(days_to_get)
      csv_file = csv_path ## <---- This can be a list of paths
      today = True
   else:
      # csv file passed as system argument
      csv_file = sys.argv[1]  ## <--- A single file path
      print('Using existing csv file:')
      print('\t',csv_file) 
      today = False
   
   #setup automatic start of simulations
   if 'auto' in str(sys.argv):
      sf.print_question('Autostart detected. Is this what you want? yes/no, default = [no]')
      auto_start = sf.read_boolean('no')
      if auto_start:
         sf.print_question('How many simulations to run at once? default = [4]')
         num_starts = sf.read_integer(4)
      else:
         num_starts = -1
   else:
      auto_start = False
      num_starts = -1
   if auto_start:
      print('Simulations will be started automatically. Make sure you have the resources.')
      time.sleep(2)
   else:
      print('Simulations will need to be started manually')
   
   ##################################################################
   # configure the job(s), uses questionaire from simple_forecast.p #
   ##################################################################
   base_cfg = make_base_configuration()
   
   ##################################################################
   # read the data from the csv file(s), get date string for naming #
   ##################################################################
   data, csv_date_str = read_NGFS_csv_data(csv_file)
   
   #initialize an object of the ngfs_day
   csv = ngfs_day(csv_date_str,today)

   '''
   if csv.today:
      print('Class object created, date is today')
      print(csv)
   else:
      print('Older csv file')
   '''

   #filter out "null incidents"
   #maybe null incidents are what we are interested in?
   data_full = data
   data = data.dropna(subset=['incident_name'])
   

   #print('Incidents in csv file')
   incident_names = data['incident_name'].unique()
   num_incidents = incident_names.shape[0]

   #number of id_strings and incidents doesn't match
   #id_strings = data['incident_id_string'].unique()
   #print('Number of id_strings: ',id_strings.shape[0])

   #list of ignition points for plotting later
   #replace by plotting based on status in the ngfs_day class
   ign_lon = []
   ign_lat = []
   
   #counter
   #new_incidents = 0
   '''
   idx = np.arange(num_incidents)
   print(idx,type(idx))
   '''
   new_idx = np.zeros((num_incidents,), dtype=int)

   ### for testing of new incident making class functions
   test_incident = ngfs_incident('Testing')
   print('test inc is a ',test_incident)
   
   #make an array of incident objects
   incidents  = [ngfs_incident(incident_names[i]) for i in range(num_incidents)]

   for i in range(num_incidents):
      print(incidents[i].name)

      #get only subset of detections from individual nifc incidents
      incident_subset = data[data.incident_name == incident_names[i]]

      ### Delete the stuff below and replace with incident_names[i].process_incident(incident_subset,data,data_full)
      incidents[i].process_incident(incident_subset,data,data_full)
      '''
      #sort by observation time	
      incident_subset = incident_subset.sort_values(by='observation_time')
      #try to make a function that takes incident_subset as argument and then does all of this below
      

      det_count = incident_subset.shape[0]
      print('\t',det_count,'detections') 
      
      #demographic data for incident
      #could span more than one county, state --> use unique to generate a list?
      county = incident_subset['locale'][incident_subset.index[0]]
      state = incident_subset['state'][incident_subset.index[0]]
      incidents[i].add_incident_location(county,state)
      loc_str = '\t'+county+', '+state
      print(loc_str)
      #print('\t',incidents[i].loc_str)
      
      #location and population data for incident
      loc_idx = pop_data[pop_data['Location'] == incidents[i].loc_str]
      
      if loc_idx.index.empty:
         pop = float('NaN')
         incidents[i].set_population(pop)
         print('\tNo matching population data found')
      else:
         pop = loc_idx.iloc[0]['Population'] #this is a string with commas
         incidents[i].set_population(float(pop.replace(',','')))  # convert string to float
         print('\tIncident county population is ',incidents[i].affected_population)
      
      #add detections to the incidnet
      incidents[i].add_detections(incident_subset)
      incidents[i].set_incident_start_time()
      incidents[i].set_incident_id_string(incident_subset['incident_id_string'][incident_subset.index[0]])
      print('\tIncident start time : ', incidents[i].incident_start_time)

      
      print('\tFinding start time and ignition location')
         
      idx = incident_subset.index[0] #gets the first row, having earliest time
      #not all csv files have terrain correction, so far
      try:
         ign_latlon = incident_subset.lat_tc[idx], incident_subset.lon_tc[idx]
      except:
         print('\tNo terrain corected lat/lon available')
         ign_latlon = incident_subset.lat[idx], incident_subset.lon[idx]
      #reduce to only unique lat/lon pairs
      unique_latlon = np.unique(ign_latlon)

      #search the entire csv for locations that have not been named yet
      try:
         dirty_subset = data_full[data_full['lon_tc'] == ign_latlon[1]]
         dirty_subset = dirty_subset[dirty_subset['lat_tc'] == ign_latlon[0]]
      except:
         dirty_subset = data_full[data_full['lon'] == ign_latlon[1]]
         dirty_subset = dirty_subset[dirty_subset['lat'] == ign_latlon[0]]
      if dirty_subset.shape[0] > incident_subset.shape[0]:
         print('\tFound addtional, earlier detections')
         ign_utc = dirty_subset.observation_time[dirty_subset.index[0]]
      else:
         ign_utc = min(incident_subset.observation_time[idx],incidents[i].incident_start_time)
         
         
      #date to start the simulation, 30 minutes before the ignition
      start_utc = utils.round_time_to_hour(ign_utc - timedelta(minutes=30))
      end_utc = start_utc + timedelta(hours=24)
      #time_utc goes into the namelist file
      time_utc = utils.utc_to_esmf(ign_utc)
      print('\tlocation: ',ign_latlon)
      print('\tutc ignition time: ', time_utc)
      #need to learn more about the FTP and how it's made
      total_frp = max(incident_subset.total_frp)
      #print('\ttotal frp: ',total_frp)
      

      
      #change the base configuration file to match the individual incidents parameters
      #this is not done for older incidents, why? does it save time? makes code harder to read
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
      incidents[i].set_json_start_code(filename,cfg['grid_code'])
      
      '''
      
      #append list of ignition pts for plotting them spatially
      ign_lon.append(incidents[i].ign_latlon[1])
      ign_lat.append(incidents[i].ign_latlon[0])
      #print(ign_lon)


      #detect a new incident from today
      if (incidents[i].incident_start_time.day == csv.timestamp.day and incidents[i].incident_start_time.month == csv.timestamp.month):
         #new_incidents += 1
         new_idx[i] = 1
         #change 'new' object attribute
         incidents[i].new_incident()
         incidents[i].make_incident_configuration(base_cfg)
         print('\tNew Incident From Today')
         json.dump(incidents[i].cfg, open(incidents[i].filename, 'w'), indent=4, separators=(',', ': '))
         print(('INT to start the simulation, execute ./forecast.sh %s' % incidents[i].filename))
         

      
   print('Number of incidents in csv file: ',num_incidents)
   print('Number of new incidents: ',sum(new_idx))
   
   
   #prioritize and autostart the new incidents
   #this functionallity should be put in the ngfs_day class
   sort_idx, started = prioritize_incidents(incidents,new_idx,num_starts)

   #m = Basemap(width=1200000,height=1200000,projection='lcc',resolution='l',lat_0=40.0,lon_0=-121.0)
   m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
        projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
   m.latlon = True 
   #m.bluemarble()
   m.shadedrelief()
   m.drawstates()
   m.drawcountries()
   m.drawcoastlines()
   #m.drawmapboundary(fill_color='blue')
   #m.fillcontinents(lake_color='blue')
   #print(ign_lon,ign_lat)
   ign_lon = np.array(ign_lon)
   ign_lat = np.array(ign_lat)
   x, y = m(ign_lon,ign_lat)
   #for new incidents
   new_ign_lon = ign_lon[new_idx==1]
   new_ign_lat = ign_lat[new_idx==1]
   print('new_idx',new_idx)
   x_new,y_new = m(new_ign_lon,new_ign_lat)

   if not len(ign_lon):
      print('No new incidents to add to map')
   #for started sims
   print('started',started)
   started_lon = ign_lon[started==1]  
   started_lat = ign_lat[started==1]
   x_started,y_started = m(started_lon,started_lat)
   m.scatter(x,y,s=15,label='Ongoing Incident')
   m.scatter(x_new,y_new,s=30,label='New Incident')
   m.scatter(x_started,y_started,s=30,label='New, Forecast Started')
   t_str = 'NIFC Fire Incidents \n' + csv_date_str
   plt.title(t_str)
   plt.legend(loc ="lower left")
   #save the figure
   save_str = 'NGFS_map_'+csv_date_str+'.png'
   plt.savefig(save_str,bbox_inches='tight')
   #plt.show
   
