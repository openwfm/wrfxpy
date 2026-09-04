"""
One run of the NGFS forecasting cycle.

The ngfs_day class owns everything belonging to a single execution: the
detection data, the incidents built from it, and the ledger of incident ids
already forecast. The entry point calls its methods in order -- acquire data,
add_incidents, process_incidents, start_incidents, save_outputs.

Two detection streams are kept separate. self.data holds GOES; self.viirs_data
holds NGFS VIIRS plus a NASA FIRMS URT top-up for detections NGFS has not
published yet. Incident ids appearing only in the VIIRS stream are picked up by
add_incidents as VIIRS-only incidents, which the GOES-driven production
monolith cannot see.

self.data is carried forward from the previous run's saved state rather than
re-downloaded from scratch, and deduplicated on
(latitude, longitude, acq_date_time), so an incident's detection history
survives across cycles. save_pickle trims it to the most recent 48 hours before
writing.

add_incidents is where old and new state meet: for each incident id in the data
it either revives the ngfs_incident object from the previous run and appends new
detections to it, or builds a fresh one. An id already in started_inc_ids is
built fresh but immediately marked started, which is what prevents a second
forecast for a fire that dropped out of the data and came back.

cluster_data, add_polar_data and prioritize_incidents have no callers; they are
retained from earlier versions. See README section 9.
"""
from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
sys.path.insert(1, 'src/ngfs')
import utils, logging, traceback
import json
import utils
import pandas as pd
import numpy as np
import pickle
import time
from mpl_toolkits.basemap import Basemap
import config_manager
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.cluster import DBSCAN
try:
   from pyproj import Proj, transform, Transformer
except:
   from pyproj import Proj, transform # Transformer
from PIL import Image
from datetime import timedelta, datetime
from ngfs import constants as cons
from ngfs import ngfs_api as ngfs_api
from ngfs import ngfs_ftp as ngfs_ftp
from ngfs import firms_data
from ngfs.ngfs_incident import ngfs_incident
import ngfs_dictionary as nd
import ngfs_helper as nh  #<<<<------------------------  fix this messy module
from ingest.NWS_warnings import set_red_flags
from ngfs.persistence import print_base_map as print_base_map_import
from ngfs.persistence import detection_summary as detection_summary_import
from ngfs.persistence import save_incident_text as save_incident_text_import

#to do
#make utilities.py to han to handle times, etc

#### change the download structures ####
### add_goes_data() will call download_goes_data(args), which returns the dataframe


def make_csv_date_str(sys_args):
      for sa in sys_args:
         if 'csv' in sa:
            print('CSV file passed as argument')
            #file with name like ingest/NGFS/NGFS_FIRE_DETECTIONS_GOES-19_ABI_CONUS_2026_01_12_012.csv
            return sa[-18:-8]
      return pd.Timestamp.now('UTC').date().isoformat()

def timestamp_from_string(csv_date_str):
      #string of form '/NGFS_FIRE_DETECTIONS_GOES-19_ABI_CONUS_YYYY-MM-DD_xxx.csv'
      return pd.Timestamp(year=int(csv_date_str[:4]), month=int(csv_date_str[5:7]), day=int(csv_date_str[8:10]), tz='UTC')
   
def today_is_now(ngfs_cfg):
      if 'now' in sys.argv or ngfs_cfg['run_cfg']['today_forecasts']:
         return True
      return False

def data_to_v2(df):
   #reverse the v2_to_v1 dictionary in nd
   for c in df.columns:
      if 'lon_tc' in c: #detects a typical version 1 column name
         v1_to_v2 = {v: k for k, v in nd.v2_to_v1.items()}
         df.rename(columns=v1_to_v2,inplace=True)
         return df
   return df #no version 1 column detected

def incident_distance(inc1,inc2):
   lat1 = inc1.data.latitude.mean()
   lat2 = inc2.data.latitude.mean()
   lon1 = inc1.data.longitude.mean()
   lon2 = inc2.data.longitude.mean()
   d = np.sqrt((lon2-lon1)**2+(lat2-lat1)**2)
   return d

      
def cluster_data(df,params=None):
   ##sort the data by image time
   df = df.sort_values(by = 'ign_utc')
   print('Number of points to cluster: ',len(df))
   print(df)
   #uses DBSCAN to cluster detection data wfo
   clust = pd.DataFrame()

   x,y = df.longitude,df.latitude

   #project the coordinates into someting with meters as unit
   if 'Transformer' not in dir():
      inProj = Proj(init='epsg:4326')
      outProj = Proj(init='epsg:3857')
      xp,yp = transform(inProj,outProj,np.array(x),np.array(y))
   else:
   #if correct pyproj version is available and has Transformer
      tf = Transformer.from_crs("EPSG:4326", "EPSG:3857")
      xp,yp = tf.transform(y,x)
   lonlat = np.transpose(np.array([xp,yp]))
   #dbscan parameters
   if params is None:
      min_pts = 1 
      db_eps = 2*np.sqrt(np.mean(df.pixel_area))*1000 #about 1 times average pixel resolution of 
   else:
      min_pts = params['min_pts']
      db_eps = params['db_eps']
   print('DBSCAN parameters: ',db_eps,min_pts)

   #clustering
   db = DBSCAN(eps = db_eps,min_samples = min_pts).fit(lonlat)
   db_labels = db.labels_c_cloc
   print(f'Found  {len(set(db_labels))}  clusters in the data of length {len(df)}')

   #rename the incident data for the first entry in the cluster is for the clusters and join to dataframe to return
   for i in set(db_labels):
      #noise label is -1
      if i != -1:
         tmp = df[db_labels == i]
         id = tmp.incident_id_string.unique()[0]
         name = tmp.incident_name.unique()[0]
         #assign the new names
         tmp.loc[:,'incident_id_string'] = id
         tmp.loc[:,'incident_name'] = name
         clust = clust.append(tmp)
         del tmp

   return clust


def incident_locations(incident_list):
   #returns a datframe with lat/lon of the incidents average locations
   lats = []
   lons = []
   incident_id_string = []
   incident_name = []
   ign_utc = []
   for inc in incident_list:
      lats.append(inc.data.latitude.mean())
      lons.append(inc.data.longitude.mean())
      incident_id_string.append(inc.incident_id_string)
      incident_name.append(inc.incident_name)
      ign_utc.append(inc.ign_utc)
   d = {
      'incident_name' : incident_name,
      'incident_id_string' : incident_id_string,
      'latitude' : lats,
      'longitude' : lons,
      'ign_utc' :ign_utc
   }
   return pd.DataFrame(d,index=range(len(lats))).reset_index(drop=True)

def download_viirs_data(ngfs_cfg=None,df=pd.DataFrame(),end_time = pd.Timestamp.now('UTC'),start_time = None):
      #downloads NGFS VIIRS data and then gets update via NASA FIRMS URT
      '''
      if not ngfs_cfg:
         ngfs_cfg = self.ngfs_cfg
      print(type(ngfs_cfg),ngfs_cfg)
      days_to_get = ngfs_cfg['viirs_cfg']['days_to_get']
      '''
      if not start_time:
         days_to_get = 2
         start_time = pd.Timestamp.now('UTC')-timedelta(hours=days_to_get*24)
      

      #get the NGFS detections via FTP ar API
      viirs_data = ngfs_ftp.add_ngfs_scene(ngfs_cfg,parse_times = True,sat='viirs',start_time=start_time,end_time=end_time) ### prepare for API access, if possible
      print(f'Loaded {len(viirs_data)} VIIRS detections from NGFS')

      #get the NASA FIRMS URT data, keep only new stuff not in NGFS yet
      fs =  ["noaa_20", "noaa_21", "suomi","landsat","noaa_20_Alaska", "noaa_21_Alaska","suomi_Alaska"]
      firms_days_to_get = ngfs_cfg['firms_cfg'].get('days_to_get',3)
      firms_satellites = ngfs_cfg['firms_cfg'].get('sats',fs) 
      t_max = viirs_data['acq_date'].max()
      firms_viirs = pd.DataFrame()
      for satellite in firms_satellites:
         data_read = firms_data.add_firms_urt(ngfs_cfg,satellite,csv_timestamp=None,df = pd.DataFrame())
         firms_viirs = pd.concat([firms_viirs,data_read],ignore_index = True)
         firms_viirs = firms_viirs[firms_viirs['acq_date_time']>=t_max]
      
      viirs_data = viirs_data.drop_duplicates(subset=['latitude','longitude'],keep = 'first')
      viirs_data = viirs_data.reset_index(drop=True)

      if hasattr(viirs_data,'known_incident_id'):
         known_viirs_incidents = viirs_data['known_incident_id'].unique()
         print(f'\tThere are {sum(~pd.isna(known_viirs_incidents))} incidents in the NGFS polar data, len =  {len(viirs_data)}')

      if len(df) > 0:
         viirs_data = pd.concat([df,viirs_data],ignore_index=True)

      return viirs_data 


#class that keeps atrributes about the day's data
class ngfs_day():
   '''
   The ngfs_day class provides functionality to download and organize data found within the csv files provided by the NGFS system.
   '''
   #maybe add some attributes so that ongoing, new, and started incidents can be tracked
   def __init__(self,ngfs_cfg,started_inc_ids = None,old_incidents = None,old_data = None):
      #where to store output like maps, data, etc.
      self.sys_args = sys.argv
      self.ngfs_cfg = ngfs_cfg
      self.base_cfg = config_manager.make_base_configuration(True,ngfs_cfg)
      self.ngfs_directory = ngfs_cfg['ngfs_directory']
      self.date_str = make_csv_date_str(self.sys_args)
      self.today = today_is_now(ngfs_cfg)  # <--- Boolean as to whether script was called with 'now', maybe change name to 'now'
      if self.today:
         self.timestamp = pd.Timestamp.now(tz='UTC')
      else:
         self.timestamp = timestamp_from_string(self.date_str)
      print('CSV timestamp: ',self.timestamp)

      #data source
      self.data_source = self.ngfs_cfg['goes_cfg']['data_source']

      #detection data
      self.full_data = pd.DataFrame()
      if not old_data is None:
         self.data = old_data
      else:
         self.data = pd.DataFrame()
      self.viirs_data = pd.DataFrame()
      self.sats = []
      self.known = []
      self.pop_data = None

      #incidents
      self.incidents = []
      if not old_incidents is None:
         self.old_incidents = old_incidents
      else:
         self.old_incidents = []
      if not started_inc_ids is None:
         self.started_inc_ids = started_inc_ids
      else:
         self.started_inc_ids = []
      self.start_count = 0
      self.rf_warnings = []
      self.rf_zones = []

   '''
   def __setstate__(self,state):
      defaults = {'sys_args': None,
         'ngfs_cfg': None,
         'base_cfg': None,
         'ngfs_directory': None,
         'data_source': None,
         'viirs_data': None,
         'known': None,
         'pop_data': None,
         'old_incidents': None,
         'started_inc_ids': None,
         'rf_warnings': None,
         'rf_zones': None,
         'polar': None}
      
      for key, value in defaults.items():
         if key not in self.__dict__:
            setattr(self,key,value)
   '''

      
   
   def sys_args_override(self):
      #override for command line argument 
      #add override for burn model, forecast length, 
      sys_args = self.sys_args
      if 'ftp' in sys_args:
         self.ngfs_cfg['goes_cfg']['data_source'] = 'ftp'
         self.data_source = 'ftp'
      elif 'api' in sys_args:
         self.ngfs_cfg['goes_cfg']['data_source'] = 'api'
         self.data_source = 'api'
      for sa in self.sys_args:
         if '.csv' in sa:
            self.data_source = sa
            self.today = False
      if 'now' in sys_args:
         self.ngfs_cfg['run_cfg']['today_forecasts'] = True
         self.today = True
      if 'behave' in sys_args:
         self.ngfs_cfg["fire_namelist_path"] = "etc/nlists/default.fire_behave_13"
      if 'cawfe' in sys_args:
         self.ngfs_cfg["fire_namelist_path"] = "etc/nlists/default.fire_cawfe_13"



   def add_pop_data(self):
      print('Reading in county population data')
      if 'pop_data' in self.ngfs_cfg.keys():
         pop_file = self.ngfs_cfg['pop_data']
      else:
         pop_file = 'ingest/NGFS/Population_by_US_County_July_2024.txt'
      try:
         self.pop_data = pd.read_csv(pop_file,sep='\t',encoding = "ISO-8859-1")
         #print(self.pop_data)
      except:
         print('Error reading population data')
         self.pop_data = pd.DataFrame(columns=['Location','Population'])

   def add_red_flags(self):
      #get red flag warning data
      rf_keys = ['nws_fire_wx_code', 'event_type']
      if (rf_keys[0] in self.data.keys() or rf_keys[1] in self.data.keys()):
         #empty lists are initialized for self.rf_warnings, self.rf_zones
         print('Using NGFS Fire Weather data')
      else:
         print('Using NWS Fire Weather data')
         try:
            self.rf_warnings, self.rf_zones = set_red_flags(self.date_str)
         except:
            #already empty lists
            print('Error in finding the red flag zones')

   def add_goes_data(self,df = pd.DataFrame()):
      ### loads new detection data  ###
      if len(df) > 0:
         df = data_to_v2(df) #assure its comptatible
         print(f'Starting with data of length {len(df)}')
      #find if a csv file was used in system call
      if self.data_source == 'ftp':
         print('Will download data via FTP')
         #data_read, csv_date_str = ngfs_ftp.get_ngfs_data(ngfs_cfg=self.ngfs_cfg,data=df)   # <<------ change this to the same syntax as for api download
         data_read = ngfs_ftp.add_ngfs_scene(ngfs_cfg=self.ngfs_cfg,data=df,parse_times=True)
      elif self.data_source == 'api':
         print('Will download data via API')
         data_read = ngfs_api.get_ngfs_data(ngfs_cfg=self.ngfs_cfg,data=df)
      elif 'csv' in self.data_source: #passing a csv as command line arugment
         print(f'Will use {self.data_source} for data')
         data_read, csv_date_str = ngfs_ftp.read_NGFS_csv_data([self.data_source])
      if not self.data_source:
         print('No valid data source specified')
         data_read = pd.DataFrame()

      #join data_read to existing df
      df = pd.concat([df,data_read],ignore_index=True)

      #add the data to the object
      if len(df):
         self.data = pd.concat([self.data,df],ignore_index=True)
         self.data = self.data.drop_duplicates(subset =['latitude','longitude','acq_date_time'],keep='first',ignore_index=True)
         self.sats = self.data['satellite'].unique()
         print('Added data from: ',self.sats)
      else:
         print(f'No data added via {self.data_source}')
   # this is the csv file with all entries in place.ngfs_day
         
   def add_polar_data(self,empty=False):
   
      def add_firms_data(satellite, csv_timestamp, days_to_get):   #<<------ move into the NGFS_helper module?
         try:
            self.polar.add_firms_24(sat=satellite, csv_timestamp=csv_timestamp)
            ##nrt = firms_data.add_firms_nrt(ngfs_cfg,sat,csv_timestamp=None,df = pd.DataFrame()):
         except:
            print(f'Error getting {satellite} NRT data')
         try:
            self.polar.add_firms_urt(sat=satellite, csv_timestamp=csv_timestamp)
            ##urt = firms_data.add_firms_urt(ngfs_cfg,sat,csv_timestamp=None,df = pd.DataFrame()):
         except Exception as e:
            logging.error(f'Error getting {satellite} URT data %s' % repr(e))
            traceback.print_exc()

      print('Acquiring Polar data')
      self.polar = nh.polar_data(self.timestamp,cfg=self.ngfs_cfg)   #  <<<<---------- move to init?
      if empty:
         return
      fs = self.polar.satlist #fs = self.ngfs_cfg['firms_cfg']['sats']                        <<<<<-----for removing polar
      firms_days_to_get = self.ngfs_cfg['firms_cfg'].get('days_to_get',3)
      firms_satellites = self.ngfs_cfg['firms_cfg'].get('sats',fs) 
      if self.today:
         print('\tGetting the polar data for the previous 48 hours')
      else:
         print(f'Getting the polar data for {self.date_str}, {self.timestamp.day_of_year}')
      #add dat to ngfs_day object
      for satellite in firms_satellites:
         add_firms_data(satellite, self.timestamp, firms_days_to_get)
      print(f'Size of FIRMS data: {len(self.polar.data)}')
      
      #add ngfs_viirs detections to it
      self.polar.add_ngfs_viirs(ngfs_cfg = self.ngfs_cfg)                             ### <<<<<----- this will be either ngfs_ftp or ngfs_api moudule
      '''
      end_time = pd.Timestamp.now('UTC')
      start_time = pd.Timestamp.now('UTC')-timedelta(hours=2*24)
      df = ngfs_ftp.add_ngfs_scene(ngfs_cfg,parse_times = False,sat='viirs',start_time=start_time,end_time=end_time)
      self.viirs_data = df

      '''
      print('Viirs size before sort and remove: ',len(self.polar.data))
      self.polar.data = self.polar.data.sort_values(by='acq_date')
      self.polar.data = self.polar.data.drop_duplicates(subset =['latitude','longitude','acq_date'],keep='first',ignore_index=True)
      self.polar.data = self.polar.data.reset_index(drop=True)
      print('Viirs size after sort and remove: ',len(self.polar.data))
      #print(self.polar.data.keys())
      if hasattr(self.polar.data,'known_incident_id'):
         known_viirs_incidents = self.polar.data['known_incident_id'].unique()
         print(f'\tThere are {sum(~pd.isna(known_viirs_incidents))} incidents in the NGFS polar data, polar  len =  {len(self.polar.data)}')


   def add_viirs_data(self,ngfs_cfg=None,df=pd.DataFrame(),end_time = pd.Timestamp.now('UTC'),start_time = None):
      if not ngfs_cfg:
         ngfs_cfg = self.ngfs_cfg
      #print(type(ngfs_cfg),ngfs_cfg)
      days_to_get = ngfs_cfg['viirs_cfg']['days_to_get']
      if not start_time:
         start_time = pd.Timestamp.now('UTC')-timedelta(hours=days_to_get*24)

      #download_viirs_data(ngfs_cfg=None,df=pd.DataFrame(),end_time = pd.Timestamp.now('UTC'),start_time = None):
      self.viirs_data = download_viirs_data(ngfs_cfg,df=df,end_time = end_time,start_time = start_time)
      
   
   
   def known_data(self):
      #returns a subset of the data containing only known incidents
      for k in self.data.keys():
         if 'incident_id' in k:
            incident_id_string = k
         if 'incident_name' in k:
            incident_name = k
      return self.data.dropna(subset=[incident_id_string])
   
   def unknown_incidents(self):
      #aasigns incident id and incident name for unknown "possible wildland fire" in specified list of WFO regions
      '''
      translation between'type_description' and '[type]'
      Known Wildland Fire Incident [1]
      Possible Wildland Fire Near a Solar Farm [7]
      Possible Wildland Fire Near a Persistent Emitter [8]
      Industrial [2]
      Possible Wildland Fire [0]
      Likely an Urban Source [5]
      Oil/Gas [2]
      Possible Solar Farm [4]
      Solar panel (database + spectral) [3]
      VIIRS Static Sources [2]
      Likely a Solar Farm [2]
      Volcano [2]
      '''

      if self.ngfs_cfg['unknown']['run']:
         description_types = ['Possible Wildland Fire Near a Solar Farm',
                             'Possible Wildland Fire Near a Persistent Emitter',
                             'Possible Wildland Fire']
         types = [7,8,0]
         wfo_list = self.ngfs_cfg['unknown']['wfo_list']
         print("Checking for hotspots of unknown fires")
         print(f"Checking in WFO: {wfo_list}")

         # Remove detections already associated with known incidents
         known_df = self.known_data()
         known_features = known_df.feature_tracking_id.unique()
         
         
         #subset of possible wildland fires, make boolean mask
         try:
            msk_type = self.data['type_description'].isin(description_types)
            #temp_data = self.data[self.data.type_description == 'Possible Wildland Fire']
         except:
            msk_type = self.data['type'].isin(types)
            #temp_data = self.data[self.data.type == 0]
         #subset further by specified WFO regions
         
         unknown_count = 0
         #find the feature tracking id for the posible wildland fires in weah WFO and assign a name and incident ID to each
         for w in wfo_list:
            print(f'\tLook at feature tracking ids in WFO {w}')
            w_msk = (self.data['nws_wfo_code'] == w) & msk_type
            tracking_id = self.data['feature_tracking_id'][w_msk].unique()
            skipped = 0
            unknown_count += len(tracking_id)
            print(f'\tFound {len(tracking_id)} tracking ids of possible wildland fires in WFO {w}')
            for ti in tracking_id:
               if ti in known_features:
                  #print(f'Skipping {ti}, it belongs to a known incident')
                  skipped += 1
                  continue
               print(f'\t\t{ti}')
               '''
                  Assign an knwon_incident_id, based on WFO & feature tracking id
                  tracking id: 'ID-2026-03-22T15:26:17.000Z_0019'
                  temp_id: 'I6-2026-03-22T15:26:17.000Z_0019'
                  temp_name: 'KOHX_2026_0322_617.000Z0019'
               '''
               fid = ti.replace(':', '').replace('-', '').replace('_', '')
               temp_id = f'{{{w}{fid[-4:]}-{fid[2:6]}-{fid[6:10]}-{fid[11:15]}-{fid[-12:]}}}'
               temp_name = f"{w}_{fid[2:6]}_{fid[6:10]}_{fid[-12:]}"
               print(f'\t\tincident id: {temp_id}')
               print(f'\t\tincident name: {temp_name}')
               #insert these tempoarary identifiers into the data
               self.data.loc[self.data['feature_tracking_id'] == ti,'known_incident_id'] = temp_id
               self.data.loc[self.data['feature_tracking_id'] == ti,'known_incident_name'] = temp_name
         print(f'\tFound {unknown_count-skipped} posible wildland fire with the list of WFO regions')
      else:
         print("Ignoring potential wildland fires that are unknown")

   
   def clean_data(self):
      #functions to better organzize the data/remove anything weird
      #index of known incidents
      for k in self.data.keys():
         if 'incident_id' in k:
            incident_id_string = k
         if 'incident_name' in k:
            incident_name = k
      self.known = self.data[incident_id_string].notna()
         
   def add_incidents(self):
      #older incidents
      print()
      print(f'Adding incidents. Starting with {len(self.old_incidents)} older incidents and {len(self.started_inc_ids)} incident id strings')

      #check for forecasting of unknown fire hotpots, this will add incident ids and names for possible wildland fires within WFO regions
      self.unknown_incidents()

      #find IRWIN incident id strings for new data
      try:
         inc_ids = self.data['known_incident_id'].unique()
      except:
         inc_ids = self.data['incident_id_string'].unique()  #older NGFS version]
      inc_ids = [x for x in inc_ids if type(x) == str] #keeps only string objects

      #check for incidents in the NGFS viirs data too.
      viirs_count = 0
      try:
         #viirs_inc_ids = self.polar.data['known_incident_id'].unique()       ### <<<<< remove ----------------------------------
         viirs_inc_ids = self.viirs_data['known_incident_id'].unique()
         viirs_inc_ids = [x for x in viirs_inc_ids if type(x) == str] #keeps only string objects+
         for vi in viirs_inc_ids:
            if vi not in inc_ids:
               print('\tFound VIIRS-only incident',vi)
               inc_ids.append(vi)
               viirs_count += 1
         print(f'Found {viirs_count} VIIRS-only incidents')
      except:
         print('Error reading VIIRS incident ids')

      if len(inc_ids):
         print(f'\tFound {len(inc_ids)} incident id strings in the data')
         for ii in inc_ids:
            added = False
            #get data for new or possibly unstarted incident
            data_subset = self.data[self.data['known_incident_id'] == ii]
            if len(data_subset) == 0:
               #print('\tAdding VIIRS-only incident',ii)
               #data_subset = self.polar.data[self.polar.data['known_incident_id'] == ii]    ### <<<<< remove ----------------------------------
               data_subset = self.viirs_data[self.viirs_data['known_incident_id'] == ii]
            #initialize new ngfs_incident for the data
            new_inc = ngfs_incident(name=ii,data=data_subset,base_cfg = self.base_cfg, ngfs_cfg=self.ngfs_cfg )
            #see if the incident is not new, but unstarted
            for inc in self.old_incidents:
               if (ii == inc.incident_id_string): #adds the incident object already processed
                  print(f'\tAdding incident {inc.incident_id_string} {inc.incident_name} from older pickle data')
                  if inc.started:
                     #print(f'\t\tIncident {ii} has been started')
                     new_inc.started = True
                     new_inc.ign_latlon = inc.ign_latlon #keep this for map making
                  else:
                     #print(f'\t\tIncident {ii} is unstarted')
                     #append with newst data
                     print(f'\tBefore update, incident has {len(inc.data)} data length')
                     inc.add_data(data_subset)
                     print(f'\tAfter update, incident has {len(inc.data)} data length')
                  if 'full_process' in sys.argv: # <<-------- testing remove
                     self.incidents.append(new_inc)
                  else:
                     self.incidents.append(inc)
                  added = True 
            if not added:
               if ii not in self.started_inc_ids:
                  print(f'\tAdding new incident {new_inc.incident_id_string} {new_inc.incident_name}')
               else: #this is an incident previously started but not in lastest batch of data from previous run
                  print(f'\tAdding started incident {new_inc.incident_id_string} {new_inc.incident_name}')
                  new_inc.started = True
               self.incidents.append(new_inc)
      else:
         print('No incident id strings in the data')
      print(f'\tThere are {len(self.incidents)} incidents in the data')

   def process_incidents(self):
      #runs through the list of incidents and process the new ones
      #print information about old incidents
      lookback_time = self.ngfs_cfg['run_cfg'].get('lookback_time', 24)
      print(f'Processing {len(self.incidents)} incidents')
      for inc in self.incidents:
         if not inc.started and inc.ign_utc > (self.timestamp-timedelta(hours=lookback_time)):
            print('Processing incident')
            print(f'{inc.incident_id_string},{inc.incident_name}')
            try:
               #inc.process_incident(self.data, viirs_data = self.polar.data, rf_zones = self.rf_zones, pop_data = self.pop_data)   #self.polar.data   ### <<<<< remove ----------------------------------
               inc.process_incident(self.data, viirs_data = self.viirs_data, rf_zones = self.rf_zones, pop_data = self.pop_data)   #self.viirs_data
               inc.make_incident_configuration(self.base_cfg,self.ngfs_cfg)
            except:
               print('Error processing the incident, marking it started')
               inc.started = True
         else:
            if len(self.incidents) < 200:
                  print('Skipping started or old incident')
                  try:
                     inc.print_incident()
                  except:
                     print(inc.incident_name)
                     print()
         print()

   ####### put inside of the ngfs_day object
   def prioritize_incidents(self):
      #prioritizes the new incidents to be run by county population or by total population
      #new_idx is boolean mask
      #nums_start is number of the new simulations to start set tpo be -1 
      num_starts = self.ngfs_cfg['run_cfg']['num_starts']
      incidents = self.incidents
      self.set_new()
      new_idx = self.new


      priority_by_population = False
      job_sleep = 150 # time to pause beteen jobs
      
      if sum(new_idx) > num_starts:
         frp_cutoff = 5e3 # will filter out low-frp incidents
      else:
         frp_cutoff = 0

      n = len(new_idx)
      started = np.zeros(n,dtype=int)
      pop = np.zeros(n)
      frp = np.zeros(n)
      domain_size = 31*np.ones(n) #default domain size
      for i in range(n):
         pop[i] = incidents[i].affected_population
         frp[i] = np.sum(incidents[i].data.frp)  # <<--------------------------------- Could be double counting if both GOES see it?
         if hasattr(incidents[i],'cfg') and not incidents[i].cfg is None:
            domain_size[i] = incidents[i].cfg['domains']['1']['domain_size'][0]
      
      
      #to get decreasing list, sort its negative
      if priority_by_population:
         sort_array = np.array([-pop,domain_size])
         #sort_idx = np.argsort(-pop)
      else:
         sort_array = np.array([-frp,domain_size])
         #sort_idx = np.argsort(-frp)
      
      #sort jobs by decreasing priority, but move large domain jobs to the end of the queue
      sort_array = np.transpose(sort_array)
      sort_idx = np.lexsort((sort_array[:,0],sort_array[:,1]))


      if num_starts > 0:
         print('New incidents by priority')
         start_count = 0
         for i in sort_idx:
            if incidents[i].new and not incidents[i].started:
               print(incidents[i].incident_name)
               print('\t',incidents[i].ign_latlon)
               print('\t',incidents[i].ign_utc)
               print('\t Total FRP: ',np.sum(incidents[i].data.frp))
               print('\t Population affected: ',incidents[i].affected_population)
               print('\t',incidents[i].county,', ',incidents[i].state)
               
               #detect prescribed burn
               rx = 'RX' in incidents[i].incident_name

               #don't start RX incidents or low FRP incidents under certain circumstances, frp cutoff is zero when there are not many incidents
               frp_filtered = (not rx or np.sum(incidents[i].data.frp) > frp_cutoff) 
               #start anything with viirs detections, but don't wait more than 4 hours for them
               viirs_filtered  = (len(incidents[i].viirs_data ) > 0 ) or (incidents[i].ign_utc + timedelta(hours=4) < pd.Timestamp.now(tz='UTC'))
            #print('\t',incidents[i].json_start_code)
               if start_count < num_starts and frp_filtered and viirs_filtered:

                  print('\t Automatically starting this simulation')
                  #os.system('jobs')
                  #os.system(incidents[i].cmd_str)
                  #pause between jobs to help avoid crash in metgrid. <------------------
                  print(incidents[i].cmd_str)
                  time.sleep(job_sleep) # now 
                  incidents[i].started = True
                  started[i] = 1
                  start_count += 1
               else:
                  print('\t New incident, but unstarted simulation')
                  if rx:
                     print('\t Skipping low-FRP RX incident')
                  if not viirs_filtered:
                     print('\t Waiting for VIIRS detection data')
               # note started incident that is new within time frame but previously started
               #keep track of what is started
               
            if incidents[i].started:
                  started[i] = 1
   
      return sort_idx, started

   def start_incidents(self):
      #starts the incidents
      '''
      for renaiming incidents to use in initial forecast
      for inc in df.incidents:
         print(inc.incident_id_string)
         print(inc.incident_id_string[:-13] +'InitialFcast}')

         {3955EA36-346F-4238-9976-39156D1F171A}
         {3955EA36-346F-4238-9976-InitialFcast}
         {A5F72521-C511-4FA9-AFDA-C400E719496D}
         {A5F72521-C511-4FA9-AFDA-InitialFcast}
         {CAEE60E9-9828-4FCE-B32C-309CAA44461D}
         {CAEE60E9-9828-4FCE-B32C-InitialFcast}

      '''
      num_starts = self.ngfs_cfg['run_cfg'].get('num_starts', 30)
      lookback_time = self.ngfs_cfg['run_cfg'].get('lookback_time', 24)
      job_sleep = self.ngfs_cfg['run_cfg'].get('job_sleep',150)
      start_count = 0
      wait_count = 0
      print('Starting the incidents')
      now = pd.Timestamp.now('UTC')
      #loop through the unstarted incidents, start all with VIIRS ignition estimate or older that 4 hours
      for inc in self.incidents:
         print(f'{inc.incident_name} -- {inc.incident_id_string}')
         hours_old = now - inc.ign_utc
         if start_count >= num_starts:
            inc.started = False
            print(f'\tNot starting {inc.incident_name}. Exceeded maximum start count.')
            continue
         if inc.started:
            print(f'\t{inc.incident_name}, {inc.incident_id_string} is already started.')
            #inc.new = False
            continue
         if hours_old >  timedelta(hours = lookback_time):
            print(f'\tNot starting {inc.incident_name}. Ignition time older than {lookback_time} hours')
            inc.new = False
            #inc.started = False
            continue
         inc.new = True
         if (hours_old > timedelta(hours = 4)) or (len(inc.viirs_data) > 0):
            print(f'\tStarting the forecast of {inc.incident_name}, {inc.incident_id_string}')
            #print(f"hours old: {hours_old}, number of viirs detections: {len(inc.viirs_data)}")
            #print(inc.viirs_data)
            #print(f'Now:{now},ign_utc:{inc.ign_utc},viirs detections: {len(inc.viirs_data)}')
            try:
               inc.start_forecast(sleep_time = start_count*job_sleep)
            except: #older incidents may not have start methodself.ongoing
               inc.filename = 'jobs/' + inc.cfg['grid_code'] + '.json'
               inc.set_json_start_code(inc.filename,inc.cfg['grid_code'])
               #os.system(inc.json_start_code)
               print(inc.json_start_code)
            #time.sleep(job_sleep)
            inc.started = True
            start_count += 1  
         else:
            print(f'\tWaiting for VIIRS data for {inc.incident_name}, {inc.incident_id_string}')
            #check to see if an initial forcast has been made already
            old_id = inc.incident_id_string
            new_id = old_id[:-13]+'InitialFcast}'
            if new_id not in self.started_inc_ids:
               print('\tWill start an initial forecast using NIFC location or GOES location')
               inc.initial_forecast()
               self.started_inc_ids.append(inc.incident_id_string)
               inc.start_forecast(sleep_time = start_count*job_sleep)
               start_count += 1
            else:
               print('\tInitial forecast already made, awaiting VIIRS data for full forecast')
            wait_count +=1
         print()
      #loop through the incidents and mark those started
      self.start_count = start_count
      self.set_started()

      print(f'Started {start_count} of {len(self.incidents)} incidents in the data. {wait_count} incidents are awaiting VIIRS data.')
      print(f'Of the {len(self.incidents)} incidents in the data, {sum(self.started)} have been started.')


         
   def add_full_data(self,full_data):      # <<<<----------------------------------------------------- remve eventually
      self.full_data = pd.concat([self.full_data,full_data],ignore_index = True)
   def num_incidents(self):
      return len(self.incidents)
   
   #boolen arrays that show status of incidents
   def set_new(self):
      new= []
      for inc in self.incidents:
         if inc.new:
            new.append(True)
         else:
            new.append(False)
      self.new = np.array(new)

   def set_started(self):
      started = []
      for inc in self.incidents:
         if inc.started:
            started.append(True)
         else:
            started.append(False)
      self.started = np.array(started)
   def set_ongoing(self):
      self.ongoing = ~self.new

      
   def set_save_name(self):
      if self.today:
         self.sat_name = 'GOES 18 & 19'
         self.map_save_str = self.ngfs_directory+'/NGFS__testing'+self.date_str+'.png'   ## <<<<<--------------------------------------- remove
         # Current UTC time string: YYYYMMDD_HH_MM
         now = self.timestamp
         time_str = (
            f"{now.year}"
            f"{str(now.month).zfill(2)}"
            f"{str(now.day).zfill(2)}_"
            f"{str(now.hour).zfill(2)}_"
            f"{str(now.minute).zfill(2)}"
         )
         self.pickle_save_str = f'{self.ngfs_directory}/pkl_ngfs_day_{self.date_str}_{time_str}_testing{cons.PICKLE_SUFFIX}' ## <<<<<--------------------------------------- remove
      else:
         try:
            self.sat_name = self.sats[0].replace('-','_')
         except:
            self.sat_name = 'No_Incidents'
         self.map_save_str = f'ngfs/NGFS_{self.date_str}_{self.sat_name}_testing.png'                 ## <<<<<--------------------------------------- remove
         self.pickle_save_str = f'ngfs/pkl_ngfs_day_{self.date_str}_{self.sat_name}_testing{cons.PICKLE_SUFFIX}'

   def incident_ign_latlons(self):
      return np.array([incident.ign_latlon for incident in self.incidents])

   def save_pickle(self,):
      #only keep data back 48 hours
      #find time key for older dataframes
      for k in self.data.keys():
         if 'time' in k:
            time_key = k
            break
      cutoff_time  = pd.Timestamp.now('UTC') - timedelta(hours = 48)
      self.data = self.data[self.data[time_key] > cutoff_time]
      print('Saving as ',self.pickle_save_str)
      #write to a private temporary name and rename into place, so an interrupted
      #save cannot leave a truncated pickle where the next run looks for state
      tmp_save_str = f'{self.pickle_save_str}.{os.getpid()}.tmp'
      pd.to_pickle(self,tmp_save_str,compression=cons.PICKLE_COMPRESSION,protocol=3)
      os.replace(tmp_save_str,self.pickle_save_str)


   def save_incident_text(self):
      save_incident_text_import(self)

   def detection_summary(self):
      detection_summary_import(self)
   
   def print_base_map(self):
      print_base_map_import(self)

      #return det_summary
   def save_outputs(self):
      #print map, save pickle file etc
      self.detection_summary()
      self.set_save_name()
      self.print_base_map()
      self.save_incident_text()
      if self.start_count > 0: #only save pickle file if an incident has been forecast
         self.save_pickle()

   

   
         
if __name__ == '__main__':
   print('Testing ngfs_day')
   ngfs_cfg, wrfxpy_cfg = config_manager.load_cfgs() #in __main__, for testing
   csv_date_str = '2026_02_06'
   print(sys.argv)
   nd = ngfs_day(ngfs_cfg,sys_args=sys.argv)
   df = pd.read_csv('ingest/NGFS/NGFS_FIRE_DETECTIONS_GOES-19_ABI_CONUS_2026_02_06_037.csv')
   nd.add_data(df)
   print(nd.data)
   print(len(nd.data))
   print(len(nd.incidents))
   nd.add_incidents()
   
