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
import os, sys, glob
import json
import time
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import timedelta, datetime
#add src directory to path json time
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from ingest.downloader import download_url
#from ngfs_dictionary import ngfs_dictionary  <-- use to force csv data columns into a type
import ngfs_helper as nh
import utils
#dictionary to convert between "CA" and "California", etc
import state_names as sn
import simple_forecast as sf
import ngfs_dictionary as nd
#import v2_to_v1_dict as v2d

####### Classes  ######


#class that keeps atrributes about the day's data
class ngfs_day():
   #maybe add some attributes so that ongoing, new, and started incidents can be tracked
   def __init__(self,csv_date_str,today):
      self.date_str = csv_date_str # <<---------------------------------------- rename ????
      if today:
         self.timestamp = pd.Timestamp.now(tz='UTC')
         print('CSV timestamp: ',self.timestamp)
      else:
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

   def save_incident_text(self):
      #saves the ignition point information for all csv events as a csv file. Needs to add 
      # additional fields to the csv file such as viirs pixel location and assumed ignition time
      # use pd.to_csv() function
      ign_pix = pd.DataFrame()
        #print(self.ngfs_days)
      new_ign_lat = np.array([])
      new_ign_lon = np.array([])
      new_ign_time = list()
      viirs_pixel = list()
      
         #print(ngfs_day)
      for inc in self.incidents:
         #print(inc.name,inc.new)
         #time.sleep(3)
         if inc.new:
            #print(inc.name,inc.new)
            #print(inc.ignition_pixel)
            ign_pix = ign_pix.append(inc.ignition_pixel)
            #for adding columns to the csv file
            try:
               new_ign_lat = np.append(new_ign_lat,inc.new_ign_latlon[0])
               new_ign_lon = np.append(new_ign_lon,inc.new_ign_latlon[1])
               new_ign_time.append(inc.ign_utc)
               viirs_pixel.append(True)
            except:
               new_ign_lat = np.append(new_ign_lat,inc.ign_latlon[0])
               new_ign_lon = np.append(new_ign_lon,inc.ign_latlon[1])
               new_ign_time.append(inc.ign_utc)
               viirs_pixel.append(False)
      #print('length of new vars',len(new_ign_lat),len(new_ign_lon),len(new_ign_time))
      #print('dataframe shape',ign_pix.shape)
      ign_pix['forecast_ign_lat'] = new_ign_lat
      ign_pix['forecast_ign_lon'] = new_ign_lon
      ign_pix['forecast_ign_UTC'] = new_ign_time
      ign_pix['viirs_pixel_ign'] = viirs_pixel
      print(ign_pix)
      csv_save_str = 'ngfs/forecast_ignition_pixels_'+self.date_str+'.csv'
      ign_pix.to_csv(csv_save_str,index=False)

   def print_base_map(self):
      print('Starting map making')
      m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
         projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
      m.latlon = True 
      #m.bluemarble()
      m.shadedrelief()
      m.drawstates()
      m.drawcountries()
      
 
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
      
   #these are the columns of the csv file as pandas dataframe
   def add_data(self,df):
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
   def set_new_ign_utc(self,new_ign_utc):
      self.new_ign_utc = new_ign_utc
      if new_ign_utc < self.ign_utc:
         print('\tChanging the time paramaters of the job for ealier detection information')
         print('\tNew UTC time:',new_ign_utc)
         print('\tOld UTC time:',self.ign_utc)
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
      except:  #not in the v2 csv at all
         print('\tNo terrain corrected corners')
         min_lat = min(min(self.data['lat_c1']),min(self.data['lat_c2']),min(self.data['lat_c3']),min(self.data['lat_c4']))
         max_lat = max(max(self.data['lat_c1']),max(self.data['lat_c2']),max(self.data['lat_c3']),max(self.data['lat_c4']))
         min_lon = min(min(self.data['lon_c1']),min(self.data['lon_c2']),min(self.data['lon_c3']),min(self.data['lon_c4']))
         max_lon = max(max(self.data['lon_c1']),max(self.data['lon_c2']),max(self.data['lon_c3']),max(self.data['lon_c4']))
      self.bbox = min_lon,max_lon,min_lat,max_lat
      print('\tBounding box: ' , self.bbox)
   
   def make_incident_configuration(self,base_cfg):
      #change the base configuration file to match the individual incidents parameters
      
      cfg = base_cfg
      try:
         print('\tGrib source: ',cfg['grib_source'])
      except:
         print('\tGrib source is unset')
      gc_string  = self.incident_name+'_'+utils.utc_to_esmf(self.start_utc)+'_'+self.incident_id_string[1:-1]
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

      #HRR needs special data handling because 48 hour forecasts are only issued 
      #   at t00z, t06z, t12z, and t18z
      cycle_start = self.start_utc
      if cfg['grib_source'] == 'HRRR':
         print('\tComputing start of grib cycle')
         cycle_hour = np.int8(np.trunc(self.start_utc.hour/6))*6
         cycle_start = cycle_start.replace(hour = cycle_hour)
      cfg['cycle_start_utc'] = utils.utc_to_esmf(cycle_start)
      cfg['download_whole_cycle'] = 'true'

         #print(cfg)
      descrip_string = self.incident_name
      cfg['postproc']['description'] = descrip_string
      self.filename = 'jobs/' + cfg['grid_code'] + '.json'
      self.set_json_start_code(self.filename,cfg['grid_code']) 
      self.cfg = cfg
      #viewprint(cfg)

   def process_incident(self,incident_data,data,full_data):

      #sort data by observation time	
      incident_data = incident_data.sort_values(by='actual_image_time')

      id_strings = incident_data['incident_id_string'].unique()
      if id_strings.shape[0] > 1:
         print('\tFound multiple id strings for named incident')
         time.sleep(10)

      #take the first id string and name from the list
      self.incident_id_string = incident_data['incident_id_string'][incident_data.index[0]]
      self.incident_name  = incident_data['incident_name'][incident_data.index[0]]
      print(self.incident_name)

      self.det_count = incident_data.shape[0]
      print('\tNumber of detections: ',self.det_count) 

      #get unique locations
      self.unique_latlons = np.asarray(incident_data.loc[:,['lat_tc','lon_tc']].drop_duplicates())
      print('\tNumber of unique detection locations: ',self.unique_latlons.shape[0])
      #print(self.unique_latlons)

      #incident ignition pixel has corners of the pixel, etc
      self.ignition_pixel = incident_data.iloc[0]
      #print('\tIgnition pixel: ')
      #pixel_corners(self.ignition_pixel)
      
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
      self.add_data(incident_data)

      self.set_incident_start_time()
      print('\tNIFC incident start time : ', self.incident_start_time)

      #get the bounding box
      self.set_incident_bounding_box()


      print('\tFinding start time and ignition location')
         
      idx = incident_data.index[0] #gets the first row, having earliest time
      #not all csv files have terrain correction, so far
      try:
         tc_ign_latlon = [incident_data.lat_tc[idx], incident_data.lon_tc[idx]]
         swir_ign_latlon = [incident_data.lat_tc_swir[idx], incident_data.lon_tc_swir[idx]]
         #average of all detections within three hours of first detection
         time_msk = incident_data.actual_image_time - incident_data.actual_image_time[idx] < timedelta(hours = 3.0)
         swir_msk = abs(incident_data.lon_tc_swir) < 180.0

         if sum(time_msk & swir_msk) > 1:
            mean_ign_latlon = [np.mean(incident_data.lat_tc_swir[time_msk & swir_msk]), np.mean(incident_data.lon_tc_swir[time_msk & swir_msk])]
            print('\tAveraging SWIR pixels from first half hour')
         else:
            mean_ign_latlon = [np.mean(incident_data.lat_tc), np.mean(incident_data.lon_tc)]
         print('\tStandard ign_latlon: ',tc_ign_latlon)
         print('\tSWIR ign_latlon : ',swir_ign_latlon)
         print('\tMean ign_latlon : ',mean_ign_latlon)

         # Calculate the absolute differences
         lat_diff = abs(tc_ign_latlon[0] - swir_ign_latlon[0])
         lon_diff = abs(tc_ign_latlon[1] - swir_ign_latlon[1])

         if max(lat_diff, lon_diff) < 0.01:
            self.ign_latlon = swir_ign_latlon
            print('\tUsing SWIR ignition location')
         else:
            self.ign_latlon = tc_ign_latlon

      except:
         print('\tNo terrain corected lat/lon available')
         self.ign_latlon = [incident_data.lat[idx], incident_data.lon[idx]]
         mean_ign_latlon = [np.mean(incident_data.lat_tc), np.mean(incident_data.lon_tc)]
      #reduce to only unique lat/lon pairs
      
      #if mean_ign is true, the ignition point will be estimated from average location of initial detections
      mean_ign = True
      if mean_ign:
         self.ign_latlon = mean_ign_latlon
         print('\tUsing mean ignition point')
   
      
      
      self.unique_latlon = np.unique(self.ign_latlon)
      
      #search the entire csv for pixel locations that have not been assigned to the incident yet
      try:
         full_subset = full_data[full_data['lon_tc'] == self.ign_latlon[1]]
         full_subset = full_subset[full_subset['lat_tc'] == self.ign_latlon[0]]
      except: #if no terrain corrected data present
         full_subset = full_data[full_data['lon'] == self.ign_latlon[1]]
         full_subset = full_subset[full_subset['lat'] == self.ign_latlon[0]]
      if full_subset.shape[0] > incident_data.shape[0]:
         print('\tFound additional, earlier detections')
         self.ign_utc = full_subset.observation_time[full_subset.index[0]]
         print('\tNew earliest pixel time: ',self.ign_utc)

      #search also for the case where t multiple incident names may exist
      id_subset = data[data['incident_id_string'] == self.incident_id_string]
      if id_subset.shape[0] > self.data.shape[0]:
         print('More data exists. Appending..')
         time.sleep(10)

      #minimum of the earliest detection pixel and stated NIFC incident start time
      self.ign_utc = min(incident_data.observation_time[idx],self.incident_start_time,incident_data.initial_observation_time[idx])
      
      #date to start the simulation, 30 minutes before the ignition
      self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=30))
      self.end_utc = self.start_utc + timedelta(hours=24)
      #time_utc goes into the namelist file
      self.time_utc = utils.utc_to_esmf(self.ign_utc)
      print('\tLocation of earliest GOES pixel in csv: ',self.ign_latlon)
      print('\tUTC ignition time: ', self.time_utc)
      #need to learn more about the FRP and how it's made
      try:
         self.total_frp = max(incident_data.frp)
      except:
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
      print('\tTime step: ',temp_cfg['domains']['1']['time_step'], 'seconds')
      
      #handling of fmda if the json file  has path to fmda "fmda_geogrid_path"
      # /data/WRFXPY/wksp_fmda/CONUS/202307/
      if 'fmda_geogrid_path' in temp_cfg:
         print('Will use FMDA for fuel moisture')
         fmda_year = temp_cfg['start_utc'][:4]
         fmda_month = temp_cfg['start_utc'][5:7]
         fmda_day = temp_cfg['start_utc'][8:10]
         fmda_hour = fmda_day = temp_cfg['start_utc'][11:13]
         base_folder = '/data/WRFXPY/wksp_fmda/CONUS/' + fmda_year + fmda_month + '/'
         date_folder = 'fmda-CONUS-'+fmda_year+fmda_month+fmda_day+'-'+fmda_hour+'.geo'
         temp_cfg['fmda_geogrid_path'] = base_folder + date_folder

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
   return pd.Timestamp(year=csv_year,month=csv_month,day=csv_day,tz='UTC')
   
    
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
   csv_str = 'NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year).zfill(3)+'.csv'
   csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-WEST/CONUS/'+csv_str
  else:
   csv_str = 'NGFS_FIRE_DETECTIONS_GOES-16_ABI_CONUS_'+str(yyyy)+'_'+str(mm).zfill(2)+'_'+str(dd).zfill(2)+'_'+str(day_of_year).zfill(3)+'.csv'
   csv_url = 'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-EAST/CONUS/'+csv_str
  
  print('\tDownloading: ',csv_str)
  csv_path = ngfs_dir+'/'+csv_str
  #downloading
  if days_previous < 2:
     download_url(csv_url,csv_path)
  else:
     if os.path.exists(csv_path):
        print('\tUsing older data in ingest path')
     else:
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
         print('NGFS data from yesterday or before')
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
   
   #columns with critical time information from v1 and v2 csv files
   time_cols =  ['incident_start_time','observation_time','initial_observation_time']
   time_cols2 = ['acq_date_time','pixel_date_time'] #for v2 csv files
   if type(csv_file) is list:
      #read the first csv file in the list
      #how to interpret the columns with all NULL, prevent reading of null data into type float when it should be a string?
      null_columns = {
         'incident_name':'string',
         'incident_conf':'string',
         'incident_type':'string'
      }
      #initialize an empty DataFrame, to be appended with NGFS csv file(s)
      data = pd.DataFrame()
      
      for i in range(len(csv_file)): 
         try: #read v2 csv
         #if True:
            data_read = pd.read_csv(csv_file[i], parse_dates=time_cols2)
            #rename v2 variables as v1 counterparts as much as possible
            data_read['initial_observation_time'] = data_read[time_cols2[1]]
            data_read['incident_start_time'] = data_read[time_cols2[1]]
            #force csv data to adhere to a specific type
            data_read.rename(columns=nd.v2_to_v1, inplace=True)
            data_read = data_read.astype(nd.v2_dict)
             
         except:
            data_read = pd.read_csv(csv_file[i], parse_dates=time_cols)
            #data_read[time_cols[1]] = pd.DatetimeIndex(pd.to_datetime(data_read[time_cols[1]])).tz_localize('UTC')
            #force csv data to adhere to a specific type
            data_read['actual_image_time'] = data_read['observation_time']
            data_read = data_read.astype(nd.ngfs_dictionary)
         #make sure the times are UTCread
         
         #data_read = data_read.astype(null_columns)
         if len(data) == 0:
            print('Reading first csv: ',csv_file[i])
            data = data_read
         else:
            print('Merging with csv file: ', csv_file[i])
            try:
               data = pd.merge(data,data_read, how = 'outer')
            except Exception as e:
               print('Fail to merge with ',csv_file[i])
               print(f"An exception occurred: {str(e)}")
               #data = pd.concat([data,data_read],axis=1,join ='outer')
         print('\tNumber of detections in csv: ',data_read.shape[0])
         print('\tTotal number of detections: ',data.shape[0])
      #for naming purposes, get only first file name in list
      #assumes standard  filename for all detections in day
      csv_date_str = csv_file[0][-18:-8]
   else:
      #print('\tReading :',csv_file)
      #time.sleep(30)
      try:
         data = pd.read_csv(csv_file, parse_dates=time_cols)
         data = data.astype(nd.ngfs_dictionary)
      except:
         data = pd.read_csv(csv_file, parse_dates=time_cols2)
         data.rename(columns=nd.v2_to_v1, inplace=True)
         data['initial_observation_time'] = data[time_cols[2]]
         data['incident_start_time'] = data[time_cols[2]]
         data = data.astype(nd.v2_dict)
      csv_date_str = csv_file[-18:-8]
   
   #make sure the times are UTC, seems to be take care of already, above
   try:
      data[time_cols[1]] = pd.DatetimeIndex(pd.to_datetime(data[time_cols[1]])).tz_localize('UTC')
   except:
      print('\tData column has correct timezone already')
   
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
            print(incidents[i].incident_name)
            print('\t',incidents[i].ign_latlon)
            print('\t',incidents[i].ign_utc)
            print('\t Population affected: ',incidents[i].affected_population)
            print('\t',incidents[i].county,', ',incidents[i].state)
         #print('\t',incidents[i].json_start_code)
            if start_count < num_starts:
               print('\t Automatically starting this simulation')
               os.system(incidents[i].cmd_str)
               #print(incidents[i].cmd_str)
               os.system('jobs')
               #pause between jobs to help avoid crash in metgrid. Working ? <------------------
               time.sleep(120)
               incidents[i].set_started()
               started[i] = 1
               start_count += 1
            else:
               print('\tNew incident, but unstarted simulation')
 
   return sort_idx, started

def pixel_corners(csv_row):
   #read the corners of the pixel and return an array
   #csv_row is the entire rown from a csv data frame
   corners = np.zeros([4,2])
   corners[0,1] = csv_row.loc['lat_c1']
   corners[1,1] = csv_row.loc['lat_c2']
   corners[2,1] = csv_row.loc['lat_c3']
   corners[3,1] = csv_row.loc['lat_c4']
   corners[0,0] = csv_row.loc['lon_c1']
   corners[1,0] = csv_row.loc['lon_c2']
   corners[2,0] = csv_row.loc['lon_c3']
   corners[3,0] = csv_row.loc['lon_c4']
   print('\t',corners)
   


if __name__ == '__main__':

   print()
   print('Starting ngfs script')
   
   print('Reading in county population data')
   pop_data = pd.read_csv('ingest/NGFS/Population_by_US_County_July_2022.txt',sep='\t',encoding = "ISO-8859-1")

   print('Reading previous 3 (?) pickle files')
   full_pick_list = glob.glob('ngfs/*.pkl')
   full_pick_list.sort(key=os.path.getmtime)
   pick_list = list()

   for i in full_pick_list:
      #print(i)
      if 'GOES' not in i:
         pick_list.append(i)
   #print(pick_list)
   #time.sleep(10)
   old_ngfs_incidents = list()
   for i in range(-1,-10,-1):
      #print(pick_list[i])
      with open(pick_list[i],'rb') as f:
         df = pickle.load(f)
      old_ngfs_incidents.extend(df.incidents)

   #for incs in old_ngfs_incidents:
   #   print(incs.name,incs.incident_name,incs.started)


   #print(old_ngfs_incidents)
   #time.sleep(10)
   
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
      days_to_get = 3
      csv_str, csv_path = download_csv_data(days_to_get)
      csv_file = csv_path ## <---- This can be a list of paths
      
   else:
      # csv file passed as system argument
      csv_file = list()
      csv_file.append(sys.argv[1])  ## <--- A single file path
      print('Using existing csv file:')
      print('\t',csv_file[0]) 
      today = False
   
   # this belongs in the ngfs_day class
   #setup automatic start of simulations
   if 'auto' in str(sys.argv):
      sf.print_question('Autostart detected. Is this what you want? yes/no, default = [no]')
      force = False
      auto_start = sf.read_boolean('no')
      if auto_start:
         sf.print_question('How many simulations to run at once? default = [25]')
         num_starts = sf.read_integer(25)
      else:
         num_starts = -1
   elif 'force' in str(sys.argv):
      #this forces auto_start and avoids questions
      force = True
      auto_start = True
      num_starts = 25   # <<<------------------------------- set to -1 for testing, change to 10 or more for operations
   else:
      force = True  # <<<------------------------------- change this to require confirmation of forecast configuration
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
   print('\tFull data shape: ',data.shape)
   
   ##################################################################
   #               initialize an object of the ngfs_day class       #
   ##################################################################
   csv = ngfs_day(csv_date_str,today)
   

   #filter out "null incidents"
   #maybe null incidents are what we are interested in?
   full_data = data
   #filter out the possible artifacts 'possible_instrument_artifact'
   # idx = data.possible_instrument_artifact == 'Y'
   # print('Found', sum(idx), ' artifact pixels')
   #try:
      # data = data[data['possible_instrument_artifact'] == 'N']
      #new data shape
      #print('\tNew data shape after filtering possible artifacts: ',data.shape)
   # except:
      #print('\tNo filter for possible artifact applied')

   data = data.dropna(subset=['incident_name'])  #  <------ change to use id_string?

   #print('Incidents in csv file')
   incident_names = data['incident_name'].unique()
   incident_id_strings = data['incident_id_string'].unique()
   num_names = incident_names.shape[0]
   num_id_strings = incident_id_strings.shape[0]
   print('Found ',num_names, 'unique incident names')
   print('Found ',num_id_strings, 'unique id strings')
   time.sleep(5)


   initialize_by_names  = False

   if initialize_by_names:
      print('\tInitializing incident object by incident names')
      num_incidents = num_names
   else:
      print('\tInitializing incident object by incident id string')
      num_incidents = num_id_strings


   
   #should belong to the ngfs_day class object
   #array for recording which incidents are new
   new_idx = np.zeros((num_incidents,), dtype=int)
   
   
   #incidents_names = [ngfs_incident(incident_names[i]) for i in range(num_names)] #  <------ change to use id_string?
   #incidents_id_strings = [ngfs_incident(incident_id_strings[i]) for i in range(num_id_strings)]

   #make an array of incident objects
   if initialize_by_names:
      incidents = [ngfs_incident(incident_names[i]) for i in range(num_names)]
   else:
      incidents = [ngfs_incident(incident_id_strings[i]) for i in range(num_id_strings)]

   print('Acquiring Polar data')
   polar = nh.polar_data(csv.timestamp)

   if csv.today:
      print('\tGetting the polar data for the previous 24 hours')
      ###current data
      try:
         polar.add_firms_24(sat='noaa_20',csv_timestamp = csv.timestamp)
      except:
         print('Error getting NOAA_20 24-hour data')
      try:
         polar.add_firms_24(sat='suomi',csv_timestamp = csv.timestamp)
      except:
         print('Error getting Suomi 24-hour data')
      try:
         polar.add_firms_24(sat='landsat',csv_timestamp = csv.timestamp)
      except:
         print('Error getting Landsat 24-hour data')
      ###archived data
      try:
         polar.add_firms_dates(sat='noaa_20',csv_timestamp = csv.timestamp,days_to_get = 3)
      except:
         print('Error getting NOAA_20 date data')
      try:
         polar.add_firms_dates(sat='suomi',csv_timestamp = csv.timestamp,days_to_get = 3)
      except:
         print('Error getting Suomi date data')
   #archived data
   else:
      print('Getting the polar data for ',csv_date_str,csv.timestamp.day_of_year)
      try:
         polar.add_firms_dates(sat='noaa_20',csv_timestamp = csv.timestamp,days_to_get = 3)
      except:
         print('Error getting NOAA_20 date data')
      try:
         polar.add_firms_dates(sat='suomi',csv_timestamp = csv.timestamp,days_to_get = 3)
      except:
         print('Error getting Suomi date data')

      
   #polar.add_modis() #<<----------- different columns than the VIIIRS dat sets 
   

   for i in range(num_incidents):
      print('')
      print(incidents[i].name)
      

      #get only subset of detections from individual nifc incidents
      if initialize_by_names:
         incident_subset = data[data.incident_name == incident_names[i]]
      else:
         incident_subset = data[data.incident_id_string == incident_id_strings[i]]

      #process all the data for an individual incident
      #print('Subsetting data for incident. Number of detections = ',incident_subset.shape[0])
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
      inc_started = False
      for old_inc in old_ngfs_incidents:
         if incidents[i].name == old_inc.name:
            print('\tKnown incident from before')
            print('\t\t',old_inc.ign_latlon)
            print('\t\t',old_inc.ign_utc)
            print('\t\t Incident started: ',old_inc.started)
            incidents[i].incident_start_time = min(incidents[i].incident_start_time,old_inc.ign_utc)
            if old_inc.started:
               inc_started = True
            #time.sleep(2)
      #start new incidents to run  change the hours back to 163
      if (csv.timestamp - incidents[i].incident_start_time) < timedelta(hours = 96) and not inc_started:
         #(incidents[i].incident_start_time.day == csv.timestamp.day and incidents[i].incident_start_time.month == csv.timestamp.month):
         new_idx[i] = 1
         #change 'new' object attribute
         incidents[i].set_new()
         incidents[i].make_incident_configuration(base_cfg)
         print('\tNew Incident From previous 24 hours')
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
   csv.add_data(data) #only csv rows with an incident name or id
   csv.add_full_data(full_data) #all csv data rows
   csv.add_incidents(incidents)
   csv.set_save_name()
   csv.set_new(new_idx.astype(bool))
   csv.set_started(started.astype(bool))
   #print('started',csv.started)
   #detrmine ongoing events
   csv.set_ongoing()
   #print('ongoing',csv.ongoing)

   #save data and pictures
   if csv.today:
      csv.save_incident_text()
   csv.save_pickle()
   csv.print_base_map()

   
   
   
