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
'''
to do list:
1) 
2) 
3) 
4) 
5) Script to get older files from VIIRS and modis
6) try to get standalone working for a case

'''
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd
import pickle
import csv
import os, sys
import json
import time
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import timedelta, datetime
#add src directory to path json time
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from downloader import download_url
#from ngfs_dictionary import ngfs_dictionary  <-- use to force csv data columns into a type
import ngfs_helper as nh
import utils
#dictionary to convert between "CA" and "California", etc
import state_names as sn
import simple_forecast as sf
import ngfs_dictionary as nd

####### Classes  ######

#class that keeps atrributes about the day's data
class ngfs_day():
   #maybe add some attributes so that ongoing, new, and started incidents can be tracked
   def __init__(self,csv_date_str,today):
      self.date_str = csv_date_str # <<---------------------------------------- rename ????
      self.timestamp = timestamp_from_string(csv_date_str)
      self.today = today  # <--- Boolean as to whether script was called with 'now', maybe change name to 'now'
      #self.set_new = list()
      #self.started_incident = list()

   def add_incidents(self,incidents):
      self.incidents = incidents  ## <--- This will hold all the information
   def add_data(self,data):
      self.data = data # this is the csv file stripped of null incidents
      self.sats = self.data['satellite_name'].unique()
      print('Data from: ',self.sats) 
   def add_full_data(self,full_data):
      self.full_data = full_data # this is the csv file will entries in place
   def num_incidents(self):
      return len(self.incidents)
   #boolen arrays that show status of incidents
   def set_ongoing(self):
      #print('started',self.started,type(self.started))
      #print('new',self.new,type(self.new))
      self.ongoing = np.logical_not(self.new + self.started)
   def set_new(self,new):
      self.new = new
   def set_started(self,started):
      self.started = started
   def set_save_name(self):
      if self.today:
         self.sat_name = 'GOES 16 & 18'
         self.map_save_str = 'ngfs/NGFS_'+csv_date_str+'.png'
         self.pickle_save_str = 'ngfs/pkl_ngfs_day_'+csv_date_str+'.pkl'
      else:
         self.sat_name = self.sats[0].replace('-','_')
         self.map_save_str = 'ngfs/NGFS_'+csv_date_str+'_'+self.sat_name+'.png'
         self.pickle_save_str = 'ngfs/pkl_ngfs_day_'+csv_date_str+'_'+self.sat_name+'.pkl'

      
   def incident_ign_latlons(self):
      ign_latlons = np.zeros([self.num_incidents(),2])
      for i in range(self.num_incidents()):
         ign_latlons[i,:] = self.incidents[i].ign_latlon
      return ign_latlons
   
   def save_pickle(self):
      print('Saving as ',self.pickle_save_str)
      with open(self.pickle_save_str,'wb') as f:
         pickle.dump(self,f)
      #make script to run through all saved files to plot


   def print_base_map(self):
      print('Starting map making')
      m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
         projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
      m.latlon = True 
      #m.bluemarble()
      m.shadedrelief()
      m.drawstates()
      m.drawcountries()
      m.drawcoastlines()
 
      #for i in range(self.num_incidents()):
      #   print(i,self.incidents[i].name,self.incidents[i].started)

      #all incident ignition locations
      latlons = self.incident_ign_latlons()

      #get ongoing incidents
      ongoing_latlons = latlons[self.ongoing,:]
      ongoing_lon = ongoing_latlons[:,1]
      ongoing_lat= ongoing_latlons[:,0]
      x_ongoing, y_ongoing = m(ongoing_lon,ongoing_lat)
      
      #new incidents
      new_latlons = latlons[self.new,:]
      new_lon = new_latlons[:,1]
      new_lat= new_latlons[:,0]
      x_new, y_new = m(new_lon,new_lat)

      #started incidents
      started_latlons = latlons[self.started,:]
      started_lon = started_latlons[:,1]
      started_lat= started_latlons[:,0]
      x_started, y_started = m(started_lon,started_lat)

      #optionally scatter the ongoing, new, and started incident locations
      if any(self.ongoing):
         m.scatter(x_ongoing,y_ongoing,s=15,label='Ongoing Incident',edgecolors='black')
      if any(self.new):
         m.scatter(x_new,y_new,s=30,label='New Incident',edgecolors='black')
      if any(self.started):
         m.scatter(x_started,y_started,s=30,label='New, Forecast Started',edgecolors='black')
      
      
      plt.legend(loc ="lower left")

      #title and save the figure
      t_str = 'NIFC Fire Incidents \n' + csv_date_str + '\n ' + self.sat_name
      plt.title(t_str)

      print('Saving ',self.map_save_str)
      plt.savefig(self.map_save_str,bbox_inches='tight')
      #plt.show
      



class ngfs_incident():
   def __init__(self,name):
      self.name = name
      self.data = pd.DataFrame()
      self.new = False
      self.started = False
      self.affected_population = 0
      self.force = False
      self.auto_start = False
   # def process_incident(self,data):
   #    #takes in DataFrame and fills in details about incident   <<<----- done below already ??
   #    pass
   def set_new(self):
      self.new = True
   def set_started(self):
      self.started = True
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
         #self.incident_start_time = self.data.iloc[0,3]
         self.incident_start_time = min(self.data['incident_start_time'])
      else:
         self.incident_start_time = 'Unset start time'
   def set_json_start_code(self,json_file,grid_code):
      self.json_start_code = './forecast.sh '+json_file+' &> logs/'+grid_code.replace(' ','_')+'.log &'      
   #def add_ignition(self):
   # These attributes are set when viirs data finds detection pixels within the goes ignition pixel
   def set_new_ign_latlon(self,new_ign_latlon):
      self.new_ign_latlon = new_ign_latlon
   def set_new_ign_utc(self,new_gn_utc):
      self.new_ign_utc = new_ign_utc
      if new_ign_utc < self.ign_utc:

         print('new',new_ign_utc)
         print('old',self.ign_utc)

         print('Changing the time paramaters of the job for ealier detection information')
         self.ign_utc = new_ign_utc
         #the start and finish times of the simulation get adjusted too
         self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=30))
         self.end_utc = self.start_utc + timedelta(hours=24)

   
   
   def set_incident_bounding_box(self):
      #returns bounding box 
      try:
         min_lat = min(min(self.data['lat_tc_c1']),min(self.data['lat_tc_c2']),min(self.data['lat_tc_c3']),min(self.data['lat_tc_c4']))
         max_lat = max(max(self.data['lat_tc_c1']),max(self.data['lat_tc_c2']),max(self.data['lat_tc_c3']),max(self.data['lat_tc_c4']))
         min_lon = min(min(self.data['lon_tc_c1']),min(self.data['lon_tc_c2']),min(self.data['lon_tc_c3']),min(self.data['lon_tc_c4']))
         max_lon = max(max(self.data['lon_tc_c1']),max(self.data['lon_tc_c2']),max(self.data['lon_tc_c3']),max(self.data['lon_tc_c4'])) 
      except:
         print('\tNo terrain corrected corners')
         min_lat = min(min(self.data['lat_c1']),min(self.data['lat_c2']),min(self.data['lat_c3']),min(self.data['lat_c4']))
         max_lat = max(max(self.data['lat_c1']),max(self.data['lat_c2']),max(self.data['lat_c3']),max(self.data['lat_c4']))
         min_lon = min(min(self.data['lon_c1']),min(self.data['lon_c2']),min(self.data['lon_c3']),min(self.data['lon_c4']))
         max_lon = max(max(self.data['lon_c1']),max(self.data['lon_c2']),max(self.data['lon_c3']),max(self.data['lon_c4']))
      self.bbox = min_lon,min_lat,max_lon,max_lat
      print('\tBounding box: ' , self.bbox)
   
   def make_incident_configuration(self,base_cfg):
      #change the base configuration file to match the individual incidents parameters
      
      cfg = base_cfg
      gc_string  = self.name+'_'+utils.utc_to_esmf(self.start_utc)+'_'+self.incident_id_string[1:-1]
      cfg['grid_code'] = gc_string.replace(' ','_')
      cfg['domains']['1']['center_latlon'] = self.ign_latlon
      cfg['domains']['1']['truelats'] = (self.ign_latlon[0], self.ign_latlon[0])
      cfg['domains']['1']['stand_lon'] = self.ign_latlon[1]
      ignitions = cfg['ignitions']
      ign_dur = ignitions['1'][0]['duration_s']
      #change if there is a viirs detection to work with
      #the domain stays the same otherwise so comparision between viirs and goes ignitions may be examined
      if not hasattr(self,'new_ign_latlon'):
         cfg['ignitions'] = { '1' : [ { 'time_utc' : self.time_utc,
                                    'duration_s' : ign_dur,
                                    'latlon' : self.ign_latlon } ] }
      else:
         cfg['ignitions'] = { '1' : [ { 'time_utc' : self.time_utc,
                                   'duration_s' : ign_dur,
                                   'latlon' : self.new_ign_latlon } ] }  # <<---- new 
      cfg['start_utc'] = utils.utc_to_esmf(self.start_utc)
      cfg['end_utc'] = utils.utc_to_esmf(self.end_utc)
         #print(cfg)
      descrip_string = self.name
      cfg['postproc']['description'] = descrip_string
      self.filename = 'jobs/' + cfg['grid_code'] + '.json'
      self.set_json_start_code(self.filename,cfg['grid_code']) 
      self.cfg = cfg

   def process_incident(self,incident_data,data,full_data):

      #sort data by observation time	
      incident_data = incident_data.sort_values(by='observation_time')

      self.det_count = incident_data.shape[0]
      print('\tNumber of detections: ',self.det_count) 

      #incident ignition pixel has corners of the pixel, etc
      self.ignition_pixel = incident_data.iloc[0]
      #print(self.ignition_pixel)
      
      #demographic data for incident
      #could span more than one county, state --> use unique to generate a list?
      self.county = incident_data['locale'][incident_data.index[0]]
      self.state = incident_data['state'][incident_data.index[0]]
      if len(self.state) == 2:
         self.state = sn.abbrev_to_us_state[self.state]
      #matches the first column of population data text file
      self.loc_str = '.'+self.county+', '+self.state
      print('\tLocation: ',self.county,', ',self.state)
      
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
      self.add_detections(incident_data)
      self.set_incident_start_time()
      self.set_incident_id_string(incident_data['incident_id_string'][incident_data.index[0]])
      print('\tNIFC incident start time : ', self.incident_start_time)

      #get the bounding box
      self.set_incident_bounding_box()


      print('\tFinding start time and ignition location')
         
      idx = incident_data.index[0] #gets the first row, having earliest time
      #not all csv files have terrain correction, so far
      try:
         self.ign_latlon = [incident_data.lat_tc[idx], incident_data.lon_tc[idx]]
      except:
         print('\tNo terrain corected lat/lon available')
         self.ign_latlon = [incident_data.lat[idx], incident_data.lon[idx]]
      #reduce to only unique lat/lon pairs
      self.unique_latlon = np.unique(self.ign_latlon)
      
      #search the entire csv for pixel locations that have not been assigned to the incident yet
      try:
         full_subset = full_data[full_data['lon_tc'] == self.ign_latlon[1]]
         full_subset = full_subset[full_subset['lat_tc'] == self.ign_latlon[0]]
      except: #if no terrain corrected data present
         full_subset = full_data[full_data['lon'] == self.ign_latlon[1]]
         full_subset = full_subset[full_subset['lat'] == self.ign_latlon[0]]
      if full_subset.shape[0] > incident_data.shape[0]:
         print('\tFound addtional, earlier detections')
         self.ign_utc = full_subset.observation_time[full_subset.index[0]]
         print('\tNew earliest pixel time: ',self.ign_utc)

      #minimum of the earliest detection pixel and stated NIFC inceident start time
      self.ign_utc = min(incident_data.observation_time[idx],self.incident_start_time,incident_data.initial_observation_time[idx])
      
      #date to start the simulation, 30 minutes before the ignition
      self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=30))
      self.end_utc = self.start_utc + timedelta(hours=24)
      #time_utc goes into the namelist file
      self.time_utc = utils.utc_to_esmf(self.ign_utc)
      print('\tLocation of earliest GOES pixel in csv: ',self.ign_latlon)
      print('\tUTC ignition time: ', self.time_utc)
      #need to learn more about the FRP and how it's made
      self.total_frp = max(incident_data.total_frp)
      #print('\ttotal frp: ',total_frp)

   #used for creating string used for system call to start simulation
   def set_cmd_str(self,cmd_str):
      self.cmd_str = cmd_str


      
####### Functions  #######
   
def make_base_configuration(force):
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
      

      if force:
         print('Using base configuration')
         use_base_cfg = True
      else:
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
      #how to interpret the columns with all NULL ?
      null_columns = {
         'incident_name':'string',
         'incident_conf':'string',
         'incident_type':'string'
      }
      data = pd.read_csv(csv_file[0], parse_dates=time_cols)
      data = data.astype(nd.ngfs_dictionary) # <--- force data to have specified dtypes
      #data = data.astype(null_columns)
      print('Reading: ',csv_file[0])
      print('\tNumber of detections: ',data.shape[0])
      #read the rest csv files and merge
      for i in range(1,len(csv_file)):
         print('Merging with csv file: ', csv_file[i])
         data_read = pd.read_csv(csv_file[i], parse_dates=time_cols)
         data = data.astype(nd.ngfs_dictionary)
         #data_read = data_read.astype(null_columns)
         print('\tNumber of detections in csv: ',data_read.shape[0])
         try:
            data = pd.merge(data,data_read, how = 'outer') 
         except Exception as e:
            print('Fail to merge with ',csv_file[i])
            print(f"An exception occurred: {str(e)}")
            #data = pd.concat([data,data_read],axis=1,join ='outer')
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
   
####### put inside of the ngfs_day object
def prioritize_incidents(incidents,new_idx,num_starts):
   #prioritizes the new incidents to be run by county population
   #new_idx is boolean mask
   #nums_start is number of the new simulations to start set tpo be -1 
   n = len(new_idx)
   started = np.zeros(n,dtype=int)
   pop = np.zeros(n)
   for i in range(n):
      pop[i] = incidents[i].affected_population
   #to get decreasing list, sort its negative
   sort_idx = np.argsort(-pop)
   if num_starts > 0:
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
               os.system(incidents[i].cmd_str)
               #print(incidents[i].cmd_str)
               os.system('jobs')
               time.sleep(3)
               incidents[i].set_started()
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
      #working with today's data
      today = True

      print('Downloading latest csv files:')
      #configure download path, create directory if needed
      # this should be set elsewhere
      ngfs_dir = 'ingest/NGFS'
      
      #download the data
      days_to_get = 2
      csv_str, csv_path = download_csv_data(days_to_get)
      csv_file = csv_path ## <---- This can be a list of paths
      
   else:
      # csv file passed as system argument
      csv_file = sys.argv[1]  ## <--- A single file path
      print('Using existing csv file:')
      print('\t',csv_file) 
      today = False
   
   # this belongs in the ngfs_day class
   #setup automatic start of simulations
   if 'auto' in str(sys.argv):
      sf.print_question('Autostart detected. Is this what you want? yes/no, default = [no]')
      force = False
      auto_start = sf.read_boolean('no')
      if auto_start:
         sf.print_question('How many simulations to run at once? default = [10]')
         num_starts = sf.read_integer(10)
      else:
         num_starts = -1
   elif 'force' in str(sys.argv):
      #this forces auto_start and avoids questions
      force = True
      auto_start = True
      num_starts = 10    # <<<------------------------------- set to -1 for testing, change to 10
   else:
      force = False
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
   #if force = false, confirm base_cfg otherwise load it

   base_cfg = make_base_configuration(force)
   
   ##################################################################
   # read the data from the csv file(s), get date string for naming #
   ##################################################################
   data, csv_date_str = read_NGFS_csv_data(csv_file)
   
   ##################################################################
   #               initialize an object of the ngfs_day class       #
   ##################################################################
   csv = ngfs_day(csv_date_str,today)


   #filter out "null incidents"
   #maybe null incidents are what we are interested in?
   full_data = data
   data = data.dropna(subset=['incident_name'])

   #print('Incidents in csv file')
   incident_names = data['incident_name'].unique()
   num_incidents = incident_names.shape[0]
   
   #should belong to the ngfs_day class object
   new_idx = np.zeros((num_incidents,), dtype=int)
   
   #should belong to the ngfs_day class object
   #make an array of incident objects
   incidents  = [ngfs_incident(incident_names[i]) for i in range(num_incidents)]

   print('Acquiring Polar data')
   polar = nh.polar_data(csv.timestamp)
   if csv.today:   ### <<<----------------------------------------------------------- Maybe download older data too?
      print('\tGetting the polar data for the previous 24 hours')
      polar.add_firms_24(sat='noaa_20',csv_timestamp = csv.timestamp)
      polar.add_firms_24(sat='suomi',csv_timestamp = csv.timestamp)
      polar.add_firms_dates(sat='noaa_20',csv_timestamp = csv.timestamp,days_to_get = 2)
      polar.add_firms_dates(sat='suomi',csv_timestamp = csv.timestamp,days_to_get = 2)
   else:
      print('Getting the polar data for ',csv_date_str,csv.timestamp.day_of_year)
      polar.add_firms_dates(sat='noaa_20',csv_timestamp = csv.timestamp,days_to_get = 2)
      polar.add_firms_dates(sat='suomi',csv_timestamp = csv.timestamp,days_to_get = 2)

      
   #polar.add_modis() #<<----------- different columns than the VIIIRS dat sets 
   

   for i in range(num_incidents):
      print(incidents[i].name)

      #get only subset of detections from individual nifc incidents
      incident_subset = data[data.incident_name == incident_names[i]]

      #process all the data for an individual incident
      incidents[i].process_incident(incident_subset,data,full_data)

      #get polar data ignition estimate, if polar data exists
      if polar:
         new_ign_latlon, new_ign_utc = polar.best_ign_estimate(incidents[i].ignition_pixel,incidents[i].bbox)
         #print(type(new_ign_latlon))
         if not np.any(np.isnan(new_ign_latlon)):
            print('\tChanging ignition location to VIIRS data location')
            # print(type(new_ign_latlon))
            # print(type(incidents[i].ign_latlon))
            incidents[i].set_new_ign_latlon(new_ign_latlon)
            incidents[i].set_new_ign_utc(new_ign_utc)
      
      #detect a new incident from today Maybe this should go in the process stage
      if (incidents[i].incident_start_time.day == csv.timestamp.day and incidents[i].incident_start_time.month == csv.timestamp.month):
         new_idx[i] = 1
         #change 'new' object attribute
         incidents[i].set_new()
         incidents[i].make_incident_configuration(base_cfg)
         print('\tNew Incident From Today')
         json.dump(incidents[i].cfg, open(incidents[i].filename, 'w'), indent=4, separators=(',', ': '))
         print('\tTo start the simulation, execute: ')
         cmd_str = './forecast.sh ' + incidents[i].filename + ' &> logs/' + incidents[i].filename[5:-5] +'.log &'
         incidents[i].set_cmd_str(cmd_str)   
         #print(\n\t\t ./forecast.sh %s' % incidents[i].filename)
         print(cmd_str)

      
         

      
   print('Number of incidents in csv file: ',num_incidents)
   print('Number of new incidents: ',sum(new_idx))

   #prioritize and autostart the new incidents
   #this functionallity should be put in the ngfs_day class
   sort_idx, started = prioritize_incidents(incidents,new_idx,num_starts)

   #update the ngfs_day object
   #maybe some of these could be called inside each other??
   csv.add_data(data)
   csv.add_full_data(full_data)
   csv.add_incidents(incidents)
   csv.set_save_name()
   csv.set_new(new_idx)
   #print('new',csv.new)
   csv.set_started(started)
   #print('started',csv.started)
   #detrmine ongoing events
   csv.set_ongoing()
   #print('ongoing',csv.ongoing)

   #save data and pictures
   csv.save_pickle()
   csv.print_base_map()

   
   
   