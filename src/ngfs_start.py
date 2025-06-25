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
'''
to do list:
1) 
2) 
3) 
4) 
5) 
6) 

'''
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd
import pickle
import copy
import csv
import os, sys, glob
import json
import time
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.cluster import DBSCAN
from PIL import Image
from datetime import timedelta, datetime
#add src directory to path json time
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from ingest.downloader import download_url
#from ngfs_dictionary import ngfs_dictionary  <-- use to force csv data columns into a type
import ngfs_helper as nh
import utils, logging, traceback
#dictionary to convert between "CA" and "California", etc
import state_names as sn
import simple_forecast as sf
import ngfs_dictionary as nd
#import shapely
from shapely.geometry import Point, LineString, Polygon

#depending on the pyproj available
try:
   from pyproj import Proj, transform, Transformer
except:
   from pyproj import Proj, transform # Transformer
from ingest.NWS_warnings import set_red_flags, red_flag_incident, subset_red_flag_data, subset_wfo_data
#import geopandas as gpd
####### Classes  ######

#class that keeps atrributes about the day's data
class ngfs_day():
   '''
   The ngfs_day class provides functionality to download and organize data found within the csv files provided by the NGFS system.
   '''
   #maybe add some attributes so that ongoing, new, and started incidents can be tracked
   def __init__(self,ngfs_cfg,csv_date_str,today):
      #where to store output like maps, data, etc.
      self.ngfs_directory = ngfs_cfg['ngfs_directory']
      self.ngfs_cfg = ngfs_cfg
      self.date_str = csv_date_str # <<---------------------------------------- rename ????
      self.today = today  # <--- Boolean as to whether script was called with 'now', maybe change name to 'now'
      if today:
         self.timestamp = pd.Timestamp.now(tz='UTC')
         print('CSV timestamp: ',self.timestamp)
      else:
         self.timestamp = timestamp_from_string(csv_date_str)

      #detection data
      self.full_data = pd.DataFrame()
      self.data = pd.DataFrame()
      self.sats = []

      #incidents
      self.incidents = []
      
      #self.set_new = list()
      #self.started_incident = list()
      self.red_flag_warnings = list()
      self.red_flag_zones = list()
      self.started_inc_ids = list()

   def add_incidents(self,incidents):
      self.incidents =  incidents## <--- This will hold all the information
   def add_data(self,df):
      self.data = self.data.append(df,ignore_index = True)
      self.sats = self.data['satellite_name'].unique()
      print('Data from: ',self.sats) 
   # this is the csv file with all entries in place
   def add_full_data(self,full_data):
      self.full_data = self.full_data.append(full_data,ignore_index = True)
   def num_incidents(self):
      return len(self.incidents)
   #boolen arrays that show status of incidents
   def set_ongoing(self):
      #print('started',self.started,type(self.started))
      #print('new',self.new,type(self.new))
      self.ongoing = np.logical_not(self.new) # + self.started)
   def set_new(self,new):
      self.new = new
   def set_started(self,started):
      self.started = started
   def set_save_name(self):
      if self.today:
         self.sat_name = 'GOES 18 & 19'
         self.map_save_str = self.ngfs_directory+'/NGFS_'+csv_date_str+'.png'
         time_str = str(time.gmtime().tm_year)+str(time.gmtime().tm_mon).zfill(2)+str(time.gmtime().tm_mday).zfill(2)+'_'+str(time.gmtime().tm_hour).zfill(2) + '_' + str(time.gmtime().tm_min).zfill(2)
         self.pickle_save_str = self.ngfs_directory+'/pkl_ngfs_day_'+csv_date_str+'_'+time_str+'.pkl'
      else:
         try:
            self.sat_name = self.sats[0].replace('-','_')
         except:
            self.sat_name = 'No_Incidents'
         self.map_save_str = 'ngfs/NGFS_'+csv_date_str+'_'+self.sat_name+'.png'
         self.pickle_save_str = 'ngfs/pkl_ngfs_day_'+csv_date_str+'_'+self.sat_name+'.pkl'

   def incident_ign_latlons(self):
      return np.array([incident.ign_latlon for incident in self.incidents])

   def save_pickle(self,):
      print('Saving as ',self.pickle_save_str)
      with open(self.pickle_save_str,'wb') as f:
         pickle.dump(self,f)
 
   def save_incident_text(self):
      """
      Saves the estimated ignition point information for all new incidents as a CSV file.
      """
      
      #Needs to add 
      # additional fields to the csv file such as viirs pixel location and assumed ignition time
      # use pd.to_csv() function
      ign_pix = pd.DataFrame()
      new_ign_lat = []
      new_ign_lon = []
      new_ign_time = []
      viirs_pixel = []
      
      for inc in self.incidents:
         #print(inc.name,inc.new)
         if inc.new:
            #print(inc.ignition_pixel)
            ign_pix = ign_pix.append(inc.ignition_pixel)
            #for adding columns to the csv 
            try:
               new_ign_lat = np.append(new_ign_lat,inc.new_ign_latlon[0])
               new_ign_lon = np.append(new_ign_lon,inc.new_ign_latlon[1])
               new_ign_time.append(inc.ign_utc)
               viirs_pixel.append(True)
            except AttributeError:
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
      #only save non-empty dataframe
      if not ign_pix.empty:
         time_str = time.strftime("%H_%M", time.localtime())
         csv_save_str = f'ngfs/forecast_ignition_pixels_{self.date_str}_{time_str}.csv'
         ign_pix.to_csv(csv_save_str, index=False)

   def print_base_map(self):

      print('Starting map making')

      def get_coords(latlons, mask, map):
         #map is a Basemap object, this returns projected coordinates
         filtered = latlons[mask]
         return map(filtered[:, 1], filtered[:, 0])
      
      def scatter_if_any(x, y, size, label, map):
         #scatters projected coordinates and adds legend text to plot
         if len(x) > 0:
            map.scatter(x, y, s=size, label=label, edgecolors='black')

      #CONUS map
      m = Basemap(llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64, urcrnrlat=49, projection='lcc', lat_1=33, lat_2=45, lon_0=-95)
      m.latlon = True 
      m.shadedrelief()
      m.drawstates()
      m.drawcountries()
      m.drawcoastlines()

      latlons = self.incident_ign_latlons()

      #tuple with incident status, icon size, legend text label for plotting on map
      incident_data = [
         (self.ongoing, 15, 'Ongoing Incident'),
         (self.new, 30, 'New Incident'),
         (self.started * self.new, 30, 'New, Forecast Started')
      ]

      #scatter the CONUS incidents
      for data, size, label in incident_data:
         x, y = get_coords(latlons, data, m)
         scatter_if_any(x, y, size, label, m)

      plt.legend(loc="lower left")
      plt.title(f'NIFC Fire Incidents \n{self.date_str}\n {self.sat_name}')

      print('Saving ', self.map_save_str)
      plt.savefig(self.map_save_str, bbox_inches='tight')
      plt.cla()

      #if there are fires in Alaska, make a map and join it to the existing CONUS map
      if max(latlons[:, 0]) > 54.0:

         n = Basemap(llcrnrlon=-164, llcrnrlat=54, urcrnrlon=-130, urcrnrlat=73, projection='lcc', lat_1=63, lat_2=68, lon_0=-151)
         n.latlon = True 
         n.shadedrelief()
         n.drawstates()
         n.drawcountries()
         n.drawcoastlines()

         #scatter the Alaska incidents on Alaskan map
         for data, size, label in incident_data:
            x, y = get_coords(latlons, data, n)
            scatter_if_any(x, y, size, label, n)

         plt.legend(loc="lower right")
         plt.title(f'NIFC Alaska Fire Incidents \n{self.date_str}\n GOES-18')

         sv_str = self.map_save_str.replace('ngfs/', 'ngfs/Alaska_')
         print('Saving ', sv_str)
         plt.savefig(sv_str, bbox_inches='tight')

         conus_str = self.map_save_str.replace('ngfs/', 'ngfs/CONUS_')
         os.system(f'cp {self.map_save_str} {conus_str}')

         time.sleep(20)
         img0 = Image.open(self.map_save_str)
         img1 = Image.open(sv_str)

         height = img0.size[1]
         width = int(np.round(height * img1.size[0] / img1.size[1]))
         img1 = img1.resize((width, height), Image.LANCZOS)
         img1.save(sv_str)

         image_new = Image.new("RGB", (img0.size[0] + img1.size[0], height), "white")
         image_new.paste(img1, (0, 0))
         image_new.paste(img0, (width, 0))
         image_new.save(self.map_save_str)

class ngfs_incident():
   '''
   The ngfs_incident class provides functionality to organize and configure forecasts for data associated with a named incident
   '''
   def __init__(self,name):
      self.name = name
      self.data = pd.DataFrame()
      self.viirs_data = pd.DataFrame()
      self.feature_data = pd.DataFrame()
      self.new = False
      self.started = False
      self.affected_population = 0
      self.force = False
      self.auto_start = False
      self.feature_tracking_id = list()
      self.red_flag = False
      #set current time as estimate ignition time
      self.ign_utc = pd.Timestamp.now(tz='UTC') 
      self.start_utc = self.ign_utc - timedelta(hours = 2)
      self.end_utc = self.ign_utc + timedelta(hours = 4)
      self.ign_latlon = [40.0,-120.0]

   #could be more than one location, use append?
   def add_incident_location(self,county,state):
      if len(state) == 2:
         state = sn.abbrev_to_us_state[state]
      #matches the first column of population data text file
      self.loc_str = '.'+county+', '+state
      self.state = state
      self.county = county
      
   #these are the columns of the csv file as pandas dataframe
   def add_data(self,df):
      self.data = self.data.append(df,ignore_index = True)

   def add_viirs_data(self,viirs_df):
      self.viirs_data = self.viirs_data.append(viirs_df,ignore_index = True)

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
         print('\tChanging the time parameters of the job for ealier detection information')
         print('\tNew UTC time:',new_ign_utc)
         print('\tOld UTC time:',self.ign_utc)
         self.ign_utc = new_ign_utc
         #the start and finish times of the simulation get adjusted too
         self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=60))
         self.end_utc = self.start_utc + timedelta(hours=24)

   
   
   def set_incident_bounding_box(self):
      #returns bounding box 

      try:
         lat_columns = [f'lat_tc_c{i}' for i in range(1, 5)]
         lon_columns = [f'lon_tc_c{i}' for i in range(1, 5)]
         min_lat = min(min(self.data[col]) for col in lat_columns)
         max_lat = max(max(self.data[col]) for col in lat_columns)
         min_lon = min(min(self.data[col]) for col in lon_columns)
         max_lon = max(max(self.data[col]) for col in lon_columns)
      except:
         print('\tNo terrain corrected corners')
         lat_columns = [f'lat_c{i}' for i in range(1, 5)]
         lon_columns = [f'lon_c{i}' for i in range(1, 5)]
         min_lat = min(min(self.data[col]) for col in lat_columns)
         max_lat = max(max(self.data[col]) for col in lat_columns)
         min_lon = min(min(self.data[col]) for col in lon_columns)
         max_lon = max(max(self.data[col]) for col in lon_columns)

      self.bbox = min_lon, max_lon, min_lat, max_lat
      print('\tBounding box:', self.bbox)
 
   def make_incident_configuration(self,base_cfg,ngfs_cfg):
      #change the base configuration file to match the individual incidents parameters
      
      cfg = copy.deepcopy(base_cfg)
      #select the burn model via namelist.fire
      if 'fire_namelist_path' in ngfs_cfg['run_cfg'].keys():
         cfg['fire_namelist_path'] = ngfs_cfg['run_cfg']['fire_namelist_path']
         print('Using ',cfg['fire_namelist_path'], 'for namelist.fire')
         if 'behave_13' in cfg['fire_namelist_path']: 
            cfg['domains']['1']['time_step'] = 6
      
      if 'region_cfg' in ngfs_cfg.keys():
         for r in ngfs_cfg['region_cfg']:
            if self.data.state[0] in ngfs_cfg['region_cfg'][r]['state']:
               cfg['grib_source'] = ngfs_cfg['region_cfg'][r]['grib_source']
               cfg['geo_vars_path'] = ngfs_cfg['region_cfg'][r]['geo_vars_path']
               print(ngfs_cfg['region_cfg'][r]['msg'])
      else:
         update_states = ['CA','AZ','NV','UT', 'NM'] #['OR','WA','ID','MT','WY','CO']
         if self.data.state[0] in update_states:
               cfg['geo_vars_path'] = 'etc/vtables/geo_vars.json_2024'
               print('\tUsing updated Landfire maps')
         #look for Alaska/Hawaii and updated regions to use latest Landfire and appropriate weather products
         if any(self.data.state == 'AK'):
            cfg['grib_source'] = 'NAM198'
            print('\tAlaska incident detected, using NAM198 and Alaska Landfire data')
            cfg['geo_vars_path'] = 'etc/vtables/geo_vars.json_alaska'
         if any(self.data.state == 'HI'):
            cfg['grib_source'] = 'NAM196'
            print('\tHawaii incident detected, using NAM196 and Hawaii Landfire data')
            cfg['geo_vars_path'] = 'etc/vtables/geo_vars.json_hawaii'
            #maybe use the old adrjrw here because there is so much ocean in the domain?
            #cfg['wrf_namelist_path']  = "etc/nlists/default.input_adjrw"
            #cfg['fire_namelist_path'] = "etc/nlists/default.fire_adjrw"
      try:
         print('\tGrib source: ',cfg['grib_source'])
      except:
         print('\tGrib source is unset')
      #remove the { and } characters at the end of the incident id strings
      gc_string = self.incident_name+'_'+utils.utc_to_esmf(self.start_utc)+'_'+self.incident_id_string[1:-1]
      
      replace_chars = ['#','(',')',':']
      for rc in replace_chars:
         gc_string = gc_string.replace(rc,'_')
      #gc_string = gc_string.replace(')','_')
      cfg['grid_code'] = gc_string.replace(' ','_')
      cfg['domains']['1']['center_latlon'] = self.ign_latlon
      cfg['domains']['1']['truelats'] = (self.ign_latlon[0], self.ign_latlon[0])
      cfg['domains']['1']['stand_lon'] = self.ign_latlon[1]
      ignitions = cfg['ignitions']
      ign_dur = ignitions['1'][0]['duration_s']
      #change if there is a viirs detection to work with
      #the domain stays the same otherwise so comparision between viirs and goes ignitions may be examined
      if not hasattr(self,'new_ign_latlon'):
         cfg['ignitions'] = { '1' : [ { 'time_utc' : utils.utc_to_esmf(self.ign_utc),
                                    'duration_s' : ign_dur,
                                    'latlon' : self.ign_latlon } ] }
      else:
         t_utc = min(self.ign_utc,self.new_ign_utc)
         cfg['ignitions'] = { '1' : [ { 'time_utc' : utils.utc_to_esmf(t_utc),
                                   'duration_s' : ign_dur,
                                   'latlon' : self.new_ign_latlon } ] }  # <<---- new 
         cfg['domains']['1']['center_latlon'] = self.new_ign_latlon
         cfg['domains']['1']['truelats'] = (self.new_ign_latlon[0], self.new_ign_latlon[0])
         cfg['domains']['1']['stand_lon'] = self.new_ign_latlon[1]

      #check to see if domain will cover the the detection footprint
      #make a geoseries  of the ingition point and construct a buffer around it
      #see if buffer contains the bounding box, increas the size of the domain if not
      def make_detection_buffer(x,y,r):
      #make a buffer region arounf the ignition pixel
         #ctr_pt = Point(x,y)
         #convert the x,y to meter units
         if 'Transformer' not in dir():
            
            inProj = Proj(init='epsg:4326')
            outProj = Proj(init='epsg:3857')
            xp,yp = transform(inProj,outProj,x,y)
            
            print('\tProjected coords: ',xp,yp)
            #see https://pyproj4.github.io/pyproj/stable/gotchas.html#upgrading-to-pyproj-2-from-pyproj-1
            #make a Point gemoetry and buffer around it
            pt = Point(xp,yp)
            b = pt.buffer(r)
            #transform the buffer back to lat/lot in degree
            xb,yb = transform(outProj,inProj,b.exterior.coords.xy[0],b.exterior.coords.xy[1])
            
         else:
            tf = Transformer.from_crs("EPSG:4326", "EPSG:3857")
            yp,xp = tf.transform(y,x)
            pt = Point(xp,yp)
            b = pt.buffer(r)
            tg = Transformer.from_crs("EPSG:3857","EPSG:4326")
            yb,xb = tg.transform(b.exterior.coords.xy[1],b.exterior.coords.xy[0])
            
         buff = Polygon(zip(xb,yb))
         return buff
      
      x = self.ign_latlon[1]
      y = self.ign_latlon[0]
      r = 15*1000  #15 km <-- change this to take parameters from the configuration file
      size_increase = 0 # number of time the size of the domain has doubled
      buff = make_detection_buffer(x,y,r)

      #look through the other detections and resize the domain if detections are not within the buffer region
      #dets = gpd.points_from_xy(self.unique_latlons[:,1],self.unique_latlons[:,0])
      dets = tuple(zip(self.unique_latlons[:,1],self.unique_latlons[:,0]))
      for d in dets:
         if not Point(d).within(buff):
            r = 2*r
            buff = make_detection_buffer(x,y,r)
            size_increase += 1
      #double the domain size in each dimension and quadruple the proccessor count
      #example domain size [31,31] --> [61,61]
      if size_increase > 0 and size_increase < 5:
         sz = cfg['domains']['1']['domain_size'][0] - 1
         new_sz = sz*2**size_increase + 1
         cfg['domains']['1']['domain_size'] = [new_sz,new_sz]
         sr = cfg['domains']['1']['subgrid_ratio'][0]
         new_sr = int(sr*(1/2)**(size_increase-1))
         cfg['domains']['1']['subgrid_ratio'] = [new_sr,new_sr]
         ppn = cfg['ppn']
         cfg['ppn'] = min(400,ppn*2**(size_increase*2))
         print('Resizing domain, doublings = ',size_increase)

      #if size_increase > 2:
      #   print('Domain too large to forecast')
      


      cfg['start_utc'] = utils.utc_to_esmf(self.start_utc)
      cfg['end_utc'] = utils.utc_to_esmf(utils.round_time_to_hour(self.ign_utc + timedelta(hours =  25.5))) #utils.round_time_to_hour(self.ign_utc - timedelta(minutes=30))

      #HRR and others need special data handling because 48 hour forecasts are only issued 
      #   at t00z, t06z, t12z, and t18z
      cycle_start = self.start_utc
      cycle_list = ['HRRR','HRRR_AK','NAM198','NAM196']
      if cfg['grib_source'] in cycle_list:
         print('\tComputing start of grib cycle')
         cycle_hour = np.int8(np.trunc(self.start_utc.hour/6))*6
         cycle_start = cycle_start.replace(hour = cycle_hour)
         cfg['cycle_start_utc'] = utils.utc_to_esmf(cycle_start)
         #cfg['download_whole_cycle'] = 'true' 

      ##how the FMC will be handled
      non_conus = ['Alaska','Hawaii']
      #handling of fmda if the json file  has path to fmda "fmda_geogrid_path"
      # /data/WRFXPY/wksp_fmda/CONUS/202307/
      if self.state not in non_conus: #'fmda_geogrid_path' in cfg:
         print('\tWill use FMDA for fuel moisture')
         cfg['fmda_geogrid_path'] = get_fmda_path(cfg)# base_folder + date_folder
      else:
         print('\tNon-CONUS fire, will use eqilibrium FMC')
         #removes the key from dictionary
         cfg.pop('fmda_geogrid_path',None)


         #print(cfg)
      descrip_string = self.incident_name
      cfg['postproc']['description'] = descrip_string
      self.filename = 'jobs/' + cfg['grid_code'] + '.json'
      self.set_json_start_code(self.filename,cfg['grid_code']) 
      self.cfg = cfg
      del cfg
      #viewprint(cfg)
         
   def process_incident(self, incident_data, data, full_data, rf_zones):
    """
    Processes an incident by analyzing the incident data, extracting unique feature tracking IDs, 
    updating detection data, estimating start time, and checking for red flag warnings.
    """
    # Sort data by observation time
    incident_data = incident_data.sort_values(by='actual_image_time')

    # Initialize incident attributes
    self.incident_id_string = incident_data['incident_id_string'].iloc[0]
    self.incident_name = incident_data['incident_name'].iloc[0]
    self.feature_tracking_id = list(incident_data.feature_tracking_id.unique())
    print(self.incident_name)
    print(f'\tIncident has {len(self.feature_tracking_id)} feature tracking ids')

    # Search for feature tracking IDs in full data
    for f in self.feature_tracking_id:
        feature_subset = full_data[full_data.feature_tracking_id == f]
        feature_subset = feature_subset.drop(feature_subset[feature_subset.incident_id_string == self.incident_id_string].index)
        feature_subset = feature_subset[feature_subset.actual_image_time < min(incident_data.actual_image_time)]

        if not feature_subset.empty:
            print(f'\tFound {len(feature_subset)} earlier hotspots with feature tracking id: {f}')
            print(f'\tMean lat / lon: {np.mean(feature_subset.lat_tc)} / {np.mean(feature_subset.lon_tc)}')
            self.feature_data = self.feature_data.append(feature_subset, ignore_index=True)

            feature_subset = feature_subset[
                (feature_subset.lon_tc >= min(incident_data.lon_tc)) &
                (feature_subset.lon_tc <= max(incident_data.lon_tc)) &
                (feature_subset.lat_tc >= min(incident_data.lat_tc)) &
                (feature_subset.lat_tc <= max(incident_data.lat_tc))
            ]

            if not feature_subset.empty:
                feature_subset['incident_id_string'] = self.incident_id_string
                feature_subset['incident_name'] = self.incident_name
                incident_data = incident_data.append(feature_subset, ignore_index=True)
                incident_data = incident_data.sort_values(by='actual_image_time')

    self.det_count = incident_data.shape[0]
    print(f'\tNumber of detections: {self.det_count}')

    self.unique_latlons = np.asarray(incident_data[['lat_tc', 'lon_tc']].drop_duplicates())
    print(f'\tNumber of unique detection locations: {self.unique_latlons.shape[0]}')

    self.ignition_pixel = incident_data.iloc[0]

    self.county = incident_data.locale.iloc[0]
    self.state = incident_data.state.iloc[0]
    if len(self.state) == 2:
        try:
            self.state = sn.abbrev_to_us_state[self.state]
            self.loc_str = f'{self.county}, {self.state}'
        except KeyError:
            print('Error translating the state abbreviation(s)')
            self.state = 'Unknown'
            self.loc_str = 'Unknown'
    try:
       print(f'\tLocation: {self.loc_str}')
       loc_idx = pop_data[pop_data['Location'] == self.loc_str]
    except:
       print('\tUnknown location')
       loc_idx = pd.DataFrame()

    
    if loc_idx.empty:
        self.affected_population = float('NaN')
        print('\tNo matching population data found')
        #self.affected_population = 0.0
    else:
        pop = loc_idx.iloc[0]['Population']
        self.affected_population = float(pop.replace(',', '')) if isinstance(pop, str) else pop
        print(f'\tIncident county population is {self.affected_population}')

    self.add_data(incident_data)
    self.set_incident_bounding_box()

    print('\tFinding start time and ignition location')
    self.set_incident_start_time()
    print(f'\tEstimated incident start time: {self.incident_start_time}')

    idx = incident_data.index[0]
    try:
        tc_ign_latlon = [incident_data.lat_tc[idx], incident_data.lon_tc[idx]]
        swir_ign_latlon = [incident_data.lat_tc_swir[idx], incident_data.lon_tc_swir[idx]]

        time_msk = (incident_data.actual_image_time - incident_data.actual_image_time[idx] < timedelta(hours=3.0))
        #-999 used as fill values when lon_tc_swir not determined
        swir_msk = abs(incident_data.lon_tc_swir) < 180.0

        if sum(time_msk & swir_msk) > 1:
            mean_ign_latlon = [np.mean(incident_data.lat_tc_swir[time_msk & swir_msk]), np.mean(incident_data.lon_tc_swir[time_msk & swir_msk])]
            print('\tAveraging SWIR pixels from first three hours')
        else:
            mean_ign_latlon = [np.mean(incident_data.lat_tc[time_msk]), np.mean(incident_data.lon_tc[time_msk])]

        print(f'\tStandard ign_latlon: {tc_ign_latlon}')
        print(f'\tSWIR ign_latlon: {swir_ign_latlon}')
        print(f'\tMean ign_latlon: {mean_ign_latlon}')

        lat_diff = abs(tc_ign_latlon[0] - swir_ign_latlon[0])
        lon_diff = abs(tc_ign_latlon[1] - swir_ign_latlon[1])

        self.ign_latlon = swir_ign_latlon if max(lat_diff, lon_diff) < 0.04 else tc_ign_latlon

    except KeyError:
        print('\tNo terrain corrected lat/lon available')
        self.ign_latlon = [incident_data.lat[idx], incident_data.lon[idx]]
        mean_ign_latlon = [np.mean(incident_data.lat_tc), np.mean(incident_data.lon_tc)]

    self.ign_latlon = mean_ign_latlon
    print('\tUsing mean ignition point')

    self.unique_latlon = np.unique(self.ign_latlon)

    try:
        full_subset = full_data[
            (full_data['lon_tc'] == self.ign_latlon[1]) &
            (full_data['lat_tc'] == self.ign_latlon[0])
        ]
    except KeyError:
        full_subset = full_data[
            (full_data['lon'] == self.ign_latlon[1]) &
            (full_data['lat'] == self.ign_latlon[0])
        ]

    if full_subset.shape[0] > incident_data.shape[0]:
        print('\tFound additional, earlier detections')
        self.ign_utc = full_subset.observation_time.iloc[0]
        print(f'\tNew earliest pixel time: {self.ign_utc}')

    id_subset = data[data['incident_id_string'] == self.incident_id_string]
    if id_subset.shape[0] > self.data.shape[0]:
        print('More data exists. Appending...')
        time.sleep(10)

    self.ign_utc = min(
        incident_data.observation_time[idx],
        self.incident_start_time,
        incident_data.initial_observation_time[idx]
    )

    self.red_flag = red_flag_incident(self.ign_latlon[1], self.ign_latlon[0], self.ign_utc, rf_zones)
    if self.red_flag:
        print('\tIncident in a red flag warning zone')

    self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=60))
    self.end_utc = self.start_utc + timedelta(hours=24)
    self.time_utc = utils.utc_to_esmf(self.ign_utc)
    print(f'\tLocation of earliest GOES pixel in csv: {self.ign_latlon}')
    print(f'\tUTC ignition time: {self.time_utc}')

    try:
        self.total_frp = max(incident_data.frp)
    except KeyError:
        self.total_frp = max(incident_data.total_frp)


   #used for creating string used for system call to start simulation
   def set_cmd_str(self,cmd_str):
      self.cmd_str = cmd_str


      
####### Functions  #######
      
'''
def get_old_incidents(ngfs_directory):
   print('Reading previous pickle file(s)')
   full_pick_list = glob.glob(ngfs_directory + '/*.pkl')
   full_pick_list.sort(key=os.path.getmtime)
   pick_list = list()
   #filters  out pickle files generated by running only a specific csv file, keeps only newest
   current_time = time.time()
   for i in full_pick_list:
      #print(i)
      #look back one week
      if 'GOES' not in i and ((current_time - os.path.getmtime(i))/3600 < 24*7):
         pick_list.append(i)
   #print(pick_list)
   #time.sleep(10)
   #find old incidents in reverse chronological order, newest is first
   print('Number of pickle files to possibly look at: ',len(pick_list))
   pick_list.reverse()
   old_ngfs_incidents = []
   started_inc_ids = []

   for i in pick_list:
      print('\tReading ',i)
      df = pd.read_pickle(i)
      old_ngfs_incidents.extend(df.incidents)
      #will exit loop if the previous pickle file has statrted_inc_ids attribute
      if hasattr(df,'started_inc_ids'):
         started_inc_ids.extend(df.started_inc_ids)
         print('\tTracking ',len(started_inc_ids),' old incident ID strings')
         break
      else:
         print('\tPickle file does not have started_inc_ids')
      print(len(started_inc_ids))
      
      
   #print the old incidents, if not from the list, check if they are started
   print('Found the previous incidents:')
   for incs in old_ngfs_incidents:
      if incs.incident_id_string in started_inc_ids:
         incs.started = True
      elif incs.started:
         started_inc_ids.append(incs.incident_id_string)
      print(incs.incident_id_string,incs.incident_name,'Started = ',incs.started)
   
   return list(set(started_inc_ids)),old_ngfs_incidents

'''
def get_old_incidents(ngfs_directory):
    """
    Looks at older, saved ngfs_days objects and finds incidents that have already been processed or forecasted.
    """
    print('Reading previous pickle file(s)')
    
    # Get and sort pickle files by modification time
    full_pick_list = sorted(glob.glob(f'{ngfs_directory}/*.pkl'), key=os.path.getmtime)
    current_time = time.time()

    # Filter out old pickle files not generated by specific CSV files
    pick_list = [
        i for i in full_pick_list
        if 'GOES' not in i and (current_time - os.path.getmtime(i)) / 3600 < 24 * 7
    ]
    
    print(f'Number of pickle files to possibly look at: {len(pick_list)}')
    
    # Reverse the list to get the newest files first
    pick_list.reverse()
    
    old_ngfs_incidents = []
    started_inc_ids = []

    for i in pick_list:
        print(f'\tReading {i}')
        df = pd.read_pickle(i)
        old_ngfs_incidents.extend(df.incidents)

        if hasattr(df, 'started_inc_ids'):
            started_inc_ids.extend(df.started_inc_ids)
            print(f'\tTracking {len(started_inc_ids)} old incident ID strings')
            break
        else:
            print('\tPickle file does not have started_inc_ids')
        print(len(started_inc_ids))
    
    # Process and print old incidents
    print('Found the previous incidents:')
    for inc in old_ngfs_incidents:
        if inc.incident_id_string in started_inc_ids or inc.started:
            inc.started = True
            if inc.incident_id_string not in started_inc_ids:
                started_inc_ids.append(inc.incident_id_string)
        print(f'{inc.incident_id_string} {inc.incident_name} Started = {inc.started}')
    
    return list(set(started_inc_ids)), old_ngfs_incidents

   
def make_base_configuration(force,ngfs_cfg):
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
      if 'fmda_geogrid_path' in temp_cfg or ngfs_cfg['fmda_cfg']['use_fmda']:
         print('Will use FMDA for fuel moisture')
         temp_cfg['fmda_geogrid_path'] = get_fmda_path(temp_cfg)#base_folder + date_folder

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
   #base configuration file not found --> make one, using simple_forecast (sf)
   else:
      print('No base forecast configuration found, will run simple_forecast to make one...')
      # standard configuration to be used by all forecasts
      base_cfg = sf.questionnaire()
      print(base_cfg)
      # save base confg for next time
      json.dump(base_cfg, open('jobs/base_ngfs_cfg.json', 'w'), indent=4, separators=(',', ': '))

   
   return base_cfg

def timestamp_from_string(csv_date_str):
    return pd.Timestamp(year=int(csv_date_str[:4]), month=int(csv_date_str[5:7]), day=int(csv_date_str[8:10]), tz='UTC')

   
def cluster_data(df):
   ##sort the data by image time
   df = df.sort_values(by = 'actual_image_time')
   print('Number of points to cluster: ',len(df))
   print(df)
   #uses DBSCAN to cluster detection data wfo
   clust = pd.DataFrame()

   x,y = df.lon_tc,df.lat_tc

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
   min_pts = 10#could be one pixel that was seen for a while, or several pixels in different locations
   db_eps = 2*np.sqrt(np.mean(df.pixel_area))*1000 #about 1 times average pixel resolution of 
   print('DBSCAN parameters: ',db_eps,min_pts)

   #clustering
   db = DBSCAN(eps = db_eps,min_samples = min_pts).fit(lonlat)
   db_labels = db.labels_
   print('Found ',len(set(db_labels)),' clusters in the data')

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

def get_fmda_path(cfg):
   fmda_year = cfg['start_utc'][:4]
   fmda_month = cfg['start_utc'][5:7]
   fmda_day = cfg['start_utc'][8:10]
   fmda_hour = cfg['start_utc'][11:13]
   base_folder = '/data/WRFXPY/wksp_fmda/CONUS/' + fmda_year + fmda_month + '/'
   date_folder = 'fmda-CONUS-'+fmda_year+fmda_month+fmda_day+'-'+fmda_hour+'.geo'
   return base_folder + date_folder

def get_ngfs_csv(days_previous,sat,domain):
   #downloads NGFS csv file, returns file name and local paths
   #example csv name:
   #   NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv

   #name the csv file
   # Get the date components
   csv_day = datetime.now() - timedelta(days=days_previous)
   day_of_year = csv_day.timetuple().tm_yday
   yyyy, mm, dd = csv_day.year, csv_day.month, csv_day.day
   #domain = {'CONUS','Full-Disk'}

 
   #join strings to have name and url of the csv file
   # Define the base of the csv_str
   #base_str = 'NGFS_FIRE_DETECTIONS_GOES-{}_ABI_CONUS_{}_{}_{}_{}.csv'
   base_str = 'NGFS_FIRE_DETECTIONS_GOES-{}_ABI_{}_{}_{}_{}_{}.csv'

   # Determine the satellite
   satellite = str(sat)

   # Format the csv_str using satellite and other variables
   csv_str = base_str.format(satellite,domain,yyyy, str(mm).zfill(2), str(dd).zfill(2), str(day_of_year).zfill(3))

   # Define the csv_url using csv_str
   csv_url = f'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-{"WEST" if sat == 18 else "EAST"}/{domain}/{csv_str}'

  
   print('\tDownloading: ',csv_str)
   csv_path = 'ingest/NGFS/'+csv_str
   # today's and yesterday's
   if days_previous < 2:
      download_url(csv_url,csv_path)
   #older file should be the same as what is cached
   else:
      if os.path.exists(csv_path):
         print('\tUsing older data in ingest path')
      else:
         download_url(csv_url,csv_path)
         
   return csv_str, csv_path
   
def get_csv_lists(days_to_get,ngfs_cfg):
   #returns list of csv files to be obtained
   #empy lists for storing path names
   csv_str = list()
   csv_path = list()
      
   #number of days of data
   #days_to_get = 
   for i in range(days_to_get):
      print('NGFS data from today' if i == 0 else 'NGFS data from yesterday or before')
      '''
      for satellite in ngfs_cfg['goes_cfg']['sat_sectors']:
         print(f'GOES-{satellite}')
         for sector in ngfs_cfg['goes_cfg']['sat_sectors'][satellite]:
            print(sector)
            
            try:
               c_str, c_path = get_ngfs_csv(i, satellite, sector)
               csv_str.append(c_str)
               csv_path.append(c_path)
            except:
               print('Error getting the csv for GOES -',satellite,sector)
               #print(f'\tLooked for it at {c_path}')
      '''
      
      for satellite, sector in [(18, 'CONUS'), (19,'CONUS'), (18, 'Full-Disk')]: ### removed this (18, 'Full-Disk') on March 13 becuase of download error
         # fix it so download error doesn't crash, also this list should come from the ngfs.json config file
         try:
            c_str, c_path = get_ngfs_csv(i, satellite, sector)
            csv_str.append(c_str)
            csv_path.append(c_path)
         except:
            print('Error getting the csv for GOES-',satellite,sector)
      
   return csv_str, csv_path
   
def read_NGFS_csv_data(csv_file):
    """
    Reads and merges CSV file(s), assigning a date to them too.
    
    Parameters:
    - csv_file (str or list): A string path or a list of strings representing the CSV files.
    
    Returns:
    - data (pd.DataFrame): Merged data from CSV files.
    - csv_date_str (str): Date string extracted from the first file name.
    """
    # Define time columns for different CSV versions
    time_cols_v1 = ['incident_start_time', 'observation_time', 'initial_observation_time']
    time_cols_v2 = ['acq_date_time', 'pixel_date_time']
    
    # Convert csv_file to list if it's not already
    csv_file = [csv_file] if not isinstance(csv_file, list) else csv_file

    # Define dictionary for null columns
    null_columns = {
        'incident_name': 'string',
        'incident_conf': 'string',
        'incident_type': 'string'
    }

    def get_alaska_wfo(data):
      borough_to_wfo = {
      'Anchorage Municipality': 'KAFC',   # Anchorage Forecast Office
      'Fairbanks North Star Borough': 'KAFG',  # Fairbanks Forecast Office
      'Juneau City and Borough': 'KAJK',  # Juneau Forecast Office
      'Matanuska-Susitna Borough': 'KAFC',
      'Kenai Peninsula Borough': 'KAFC',
      'North Slope Borough': 'KAFG',
      'Nome Census Area': 'KAFG',
      'Ketchikan Gateway Borough': 'KAJK',
      'Yakutat City and Borough': 'KAJK',
      'Bethel Census Area': 'KAFG',
      'Dillingham Census Area': 'KAFC',
      'Kusilvak Census Area': 'KAFG',
      'Valdez-Cordova Census Area': 'KAFC',
      'Wrangell City and Borough': 'KAJK',
      'Haines Borough': 'KAJK',
      'Petersburg Borough': 'KAJK',
      'Sitka City and Borough': 'KAJK',
      'Prince of Wales-Hyder Census Area': 'KAJK',
      'Kodiak Island Borough': 'KAFC',
      'Bristol Bay Borough': 'KAFC',
      'Aleutians East Borough': 'KAFC',
      'Aleutians West Census Area': 'KAFC',
      'Denali Borough': 'KAFG',
      'Southeast Fairbanks Census Area': 'KAFG',
      'Yukon-Koyukuk Census Area': 'KAFG',
      'Lake and Peninsula Borough': 'KAFC',
      'Hoonah-Angoon Census Area': 'KAJK',
      'Skagway Municipality': 'KAJK',
      'Copper River Census Area': 'KAFC',
      'Northwest Arctic Borough': 'KAFG'
      # Add more as needed
      }
      data = data[data.state == 'AK']
      wfo_list = []
      for i in range(len(data)):
         try:
            wfo_list.append(borough_to_wfo[data['county'].iloc[i]])
         except:
            wfo_list.append('Unknown')
      data['nws_wfo_code'] = wfo_list
      return data
          

    # Initialize an empty DataFrame
    data = pd.DataFrame()

    # Loop through each CSV file
    for file_path in csv_file:
        # Modify feature_tracking_id string for GOES-16 or GOES-19 files
        if 'GOES-16' in file_path or 'GOES-19' in file_path:
            os.system(f"sed -i 's/,ID-20/,I6-20/g' {file_path}")

        try:
            # Attempt to read as v2 CSV
            data_read = pd.read_csv(file_path, parse_dates=time_cols_v2)
            if 'Full' in file_path:
                print(f'\nReading full-disk data: {file_path}')
                print(f'\tNumber of all Full-Disk detections: {len(data_read)}')
                data_read = data_read[data_read.country == 'United States']
                print(f'\tNumber of USA Full-Disk detections: {len(data_read)}')
                data_read = get_alaska_wfo(data_read)

            # Update fields to match legacy NGFS fields
            data_read['initial_observation_time'] = data_read[time_cols_v2[1]]
            data_read['incident_start_time'] = data_read[time_cols_v2[1]]
            data_read.rename(columns=nd.v2_to_v1, inplace=True)
            data_read = data_read.astype(nd.v2_dict)
        except Exception:
            # If reading as v2 fails, try reading as v1
            data_read = pd.read_csv(file_path, parse_dates=time_cols_v1)
            data_read['actual_image_time'] = data_read['observation_time']
            try:
                data_read = data_read.astype(nd.ngfs_dictionary)
            except Exception as e:
                print('Trouble parsing data types, probably early ngfs version')
                print(e)
                data_read = pd.DataFrame()

        # Merge the current data with existing data
        if data.empty:
            print(f'Reading first CSV: {file_path}')
            data = data_read
        else:
            print(f'Merging with CSV file: {file_path}')
            try:
                data = pd.merge(data, data_read, how='outer')
            except Exception as e:
                print(f"Failed to merge with {file_path}: {str(e)}")

        print(f'\tNumber of detections in CSV: {data_read.shape[0]}')
        print(f'\tTotal number of detections: {data.shape[0]}')

    # Extract date string from the first file name
    csv_date_str = csv_file[0][-18:-8]

    # Ensure times are in UTC
    try:
        data[time_cols_v1[1]] = pd.to_datetime(data[time_cols_v1[1]], utc=True)
    except Exception:
        print('\tData column has correct timezone already')

    return data, csv_date_str

   
####### put inside of the ngfs_day object
def prioritize_incidents(incidents,new_idx,num_starts):
   #prioritizes the new incidents to be run by county population or by total population
   #new_idx is boolean mask
   #nums_start is number of the new simulations to start set tpo be -1 
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
   domain_size = 31*np.ones(n)
   for i in range(n):
      pop[i] = incidents[i].affected_population
      frp[i] = np.sum(incidents[i].data.frp)  # <<--------------------------------- Could be double counting if both GOES see it?
      if hasattr(incidents[i],'cfg'):
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
               os.system(incidents[i].cmd_str)
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

def pixel_corners(csv_row):
   #read the corners of the pixel and return an array
   #csv_row is the entire rown from a csv data frame
   corners = np.zeros([4,2])

   # Fill corners array with latitude and longitude values
   for i in range(4):
        corners[i, 1] = csv_row[f'lat_c{i+1}']
        corners[i, 0] = csv_row[f'lon_c{i+1}']

   print('\t',corners)

# Setup for automatic download of CSV files
def setup_csv_download(sys_args,ngfs_cfg):
    # Determine if working with today's data
    today = 'now' in str(sys_args) or ngfs_cfg['run_cfg'].get('today_forecasts', False)
    
    if today:
        print('Downloading latest CSV files:')
        # Configure download path and create directory if needed
        ngfs_ingest = ngfs_cfg['goes_cfg'].get('ingest_directory', 'ingest/NGFS')
        days_to_get = ngfs_cfg['goes_cfg'].get('days_to_get', 2)
        print(f'\tWill download to: {utils.make_dir(ngfs_ingest)}')
        
        # Get list of files to download
        csv_str, csv_file = get_csv_lists(days_to_get,ngfs_cfg)
    else:
        # CSV file was passed as a system argument
        today = False
        csv_file = [sys.argv[1]]
        print('Using existing CSV file:')
        print(f'\t{csv_file[0]}')
    
    return today, csv_file

def setup_auto_start(sys_args,ngfs_cfg):
   # Check for auto start conditions
   if 'auto' in str(sys_args) or ngfs_cfg['run_cfg'].get('auto_start_forecasts', False):
        sf.print_question('Autostart detected. Is this what you want? yes/no, default = [no]')
        auto_start = sf.read_boolean('no')
        num_starts = sf.read_integer(25) if auto_start else -1
        force = False
   elif 'force' in str(sys.argv) or ngfs_cfg['run_cfg'].get('force_process', False):
        # Force auto start and avoid questions
        auto_start = True
        num_starts = ngfs_cfg['run_cfg'].get('num_starts', -1)
        force = True
   else:
        auto_start = False
        num_starts = -1
        force = True

   if auto_start:
      print('Simulations will be started automatically. Make sure you have the resources.')
      time.sleep(2)
   else:
      print('Simulations will need to be started manually')
    
   return auto_start, num_starts, force
   
if __name__ == '__main__':


   print()
   print('Starting ngfs script ',pd.Timestamp.now())

   print('Reading the ngfs configuration file')
   with open('etc/ngfs.json','r') as openfile:
         ngfs_cfg = json.load(openfile)
         #should walk user through this if not found  <<<-----------------------------------------------------

   ngfs_directory = ngfs_cfg['ngfs_directory']


   print('Reading the base wrfxpy configuration file')
   with open(ngfs_cfg['wrfxpy_cfg'],'r') as openfile:
         wrfxpy_cfg = json.load(openfile)
   
   #put this in the configuratoin file
   #popultion data available from https://www2.census.gov/programs-surveys/popest/tables/2020-2023/counties/totals/co-est2023-pop.xlsx
   print('Reading in county population data')
   if 'pop_data' in ngfs_cfg.keys():
      pop_file = ngfs_cfg['pop_data']
   else:
      pop_file = 'ingest/NGFS/Population_by_US_County_July_2023.txt'
   try:
      pop_data = pd.read_csv(pop_file,sep='\t',encoding = "ISO-8859-1")
   except:
      print('Error reading population data')
   
   #find any history and old incidents
   started_inc_ids, old_ngfs_incidents  = get_old_incidents(ngfs_directory)
   #print(old_inc_ids)

   #determine if the date is "today" and get list of csv NGFS files to download
   today, csv_file = setup_csv_download(sys.argv,ngfs_cfg)

   #determine the number of forecastrs to possibly run and whether these will start automatically
   auto_start, num_starts, force = setup_auto_start(sys.argv,ngfs_cfg)
   
   ##################################################################
   # configure the job(s), uses questionaire from simple_forecast.p #
   ##################################################################
   #if force = false, confirm base_cfg otherwise load it

   base_cfg = make_base_configuration(force,ngfs_cfg)
   
   ##################################################################
   # read the data from the csv file(s), get date string for naming #
   ##################################################################
   data, csv_date_str = read_NGFS_csv_data(csv_file)
   print('\tFull data shape: ',data.shape)
   
   ##################################################################
   #               initialize an object of the ngfs_day class       #
   ##################################################################
   csv = ngfs_day(ngfs_cfg,csv_date_str,today)
   
   #all data before deleting anything
   full_data = copy.deepcopy(data)

   #removes entries which are not from a known incident
   data = data.dropna(subset=['incident_name'])  #  <------ change to use id_string?

   #get red flag warning data
   try:
      rf_warnings, rf_zones = set_red_flags(csv.date_str)
   except:
      rf_warnings, rf_zones = [],[]
      print('Error in finding the red flag zones')

   ### add detections in selected NWS WFOs that have no incident_id_string, but are possible wildlnad fires ####
   try:
      wfo_list = ngfs_cfg['unknown']['wfo_list']
   except:
      wfo_list = ["KBOU","KSLC","KPDT","KMFR","KBO","KUNR","KBYZ","KGJT","KPUB"]# ['KSLC','KBOU','KPQR']

   #will cause forecasts to be made for unknown, unnamed incidents
   if 'unknown' in str(sys.argv) or ngfs_cfg['unknown']['run']:
      unknown_data = pd.DataFrame()
      #add data for unknown incidents in WFO regions
      wfo_data = subset_wfo_data(full_data,data,wfo_list)
      unknown_data = unknown_data.append(wfo_data,ignore_index=True)
      #add data for unknown fires in red flag areas
      rf_data = subset_red_flag_data(full_data,data,rf_zones)
      unknown_data = unknown_data.append(rf_data,ignore_index=True)
      #add unknown_data to the data for named incidents
      if len(unknown_data) > 0:
         data = data.append(cluster_data(unknown_data),ignore_index=True)


   #print('Incidents in csv file')
   incident_names = data['incident_name'].unique()
   incident_id_strings = data['incident_id_string'].unique()
   num_names = len(incident_names)
   num_id_strings = len(incident_id_strings)
   print('Found ',num_names, 'unique incident names')
   print('Found ',num_id_strings, 'unique id strings')
   time.sleep(5)

   #id strings seem to be less mutable than the names  <<<<---------------------------- probably remove this option, initializing by names doesn't make sense
   try:
      initialize_by_names = ngfs_cfg['run_cfg']['initialize_by_names']
   except:
      initialize_by_names  = False

   #make an array of ngfs_incident objects
   incidents = [ngfs_incident(incident_id_strings[i]) for i in range(num_id_strings)]
   num_incidents = num_id_strings

   #should belong to the ngfs_day class object
   #array for recording which incidents are new
   new_idx = np.zeros(num_incidents, dtype=int)

   print('Acquiring Polar data')
   polar = nh.polar_data(csv.timestamp)
   try:
      firms_days_to_get = ngfs_cfg['firms_cfg']['days_to_get']
      firms_satellites = ngfs_cfg['firms_cfg']['sats'] 
   except:
      firms_days_to_get = 3
      firms_satellites = ['noaa_20', 'noaa_21', 'suomi','landsat','noaa_20_Alaska', 'noaa_21_Alaska','suomi_Alaska']

   def add_firms_data(satellite, csv_timestamp, days_to_get):   #<<------ move into the NGFS_helper module?
      try:
         polar.add_firms_24(sat=satellite, csv_timestamp=csv_timestamp)
      except:
         print(f'Error getting {satellite} 24-hour data')
      try:
         polar.add_firms_urt(sat=satellite, csv_timestamp=csv_timestamp)
      except Exception as e:
         logging.error(f'Error getting {satellite} URT data %s' % repr(e))
         traceback.print_exc()
      ''' Not woking now 9-19-2024
      try:
         polar.add_firms_dates(sat=satellite, csv_timestamp=csv_timestamp, days_to_get=days_to_get)
      except:
         print(f'Error getting {satellite} date data')
      '''
   

   if csv.today:
      print('\tGetting the polar data for the previous 48 hours')
      for satellite in firms_satellites:
         add_firms_data(satellite, csv.timestamp, firms_days_to_get)
   else:
      print(f'Getting the polar data for {csv_date_str}, {csv.timestamp.day_of_year}')
      for satellite in firms_satellites:
         add_firms_data(satellite, csv.timestamp, firms_days_to_get)
      
   
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
      #try:
      incidents[i].process_incident(incident_subset,data,full_data,rf_zones)
      '''
      except:
         print('Error processing incident, saving incident data')
         incident_subset.to_csv(str(incidents[i].incident_id_string)+'_data.csv')
      '''
      #get polar data ignition estimate, if polar data exists, probably should be inside the process_incident function
      if polar:
         new_ign_latlon, new_ign_utc, viirs_ign_data = polar.best_ign_estimate(incidents[i].ignition_pixel,incidents[i].bbox)
         #print(type(new_ign_latlon))
         if not np.any(np.isnan(new_ign_latlon)):
            print('\tChanging ignition location to VIIRS data location')
            # print(type(new_ign_latlon))ingest/FIRMS/J1_VIIRS_C2_USA_contiguous_and_Hawaii_24h_1725456963512_20240904_130000.csv
            # print(type(incidents[i].ign_latlon))
            # print('\tAdding new VIIRS data, ',len(viirs_ign_data))
            incidents[i].set_new_ign_latlon(new_ign_latlon)
            incidents[i].set_new_ign_utc(new_ign_utc)
            incidents[i].add_viirs_data(viirs_ign_data)


      #detect a new incident from today Maybe this should go in the process stage
      if incidents[i].incident_id_string in started_inc_ids:
         print('\tPreviously started incident')
         incidents[i].started = True
      else:
         print('\tUnknown id string or unstarted incident:',incidents[i].incident_name,',',incidents[i].incident_id_string )
      #print information about old incidents
      for old_inc in old_ngfs_incidents:
            if incidents[i].name == old_inc.name:
               #incident in the list from previous, opened pickle file(s), but may be unstarted
               print('\tFound incident in older pickle file data')
               if old_inc.started:
                  incidents[i].started = True
               print('\t\t',old_inc.ign_latlon)
               print('\t\t',old_inc.ign_utc)
               print('\t\t Incident started: ',old_inc.started)
               incidents[i].incident_start_time = min(incidents[i].incident_start_time,old_inc.ign_utc)
            
      #time.sleep(2)
      #start new incidents to run 
      try:
         lookback_time = ngfs_cfg['run_cfg']['lookback_time']
      except:
         lookback_time = 24
      #if the incident is not started already, make new configuration
      if not incidents[i].started and ((csv.timestamp - incidents[i].incident_start_time) < timedelta(hours = lookback_time)): #
         #(incidents[i].incident_start_time.day == csv.timestamp.day and incidents[i].incident_start_time.month == csv.timestamp.month):
         print('\tNew Incident From previous ' + str(lookback_time) + ' hours')
         new_idx[i] = 1
         #change 'new' object attribute
         incidents[i].new = True
         incidents[i].make_incident_configuration(base_cfg,ngfs_cfg)
         json.dump(incidents[i].cfg, open(incidents[i].filename, 'w'), indent=4, separators=(',', ': '))
         print('\tTo start the simulation, execute: ')
         cmd_str = './forecast.sh ' + incidents[i].filename + ' &> logs/' + incidents[i].filename[5:-5] +'.log &'
         incidents[i].set_cmd_str(cmd_str)   
         #print(\n\t\t ./forecast.sh %s' % incidents[i].filename)
         print(cmd_str)

   #prioritize and autostart the new incidents
   #this functionallity should be put in the ngfs_day class
   sort_idx, started = prioritize_incidents(incidents,new_idx,num_starts)

   print('Number of incidents in csv file: ',num_incidents)
   print('Number of new incidents: ',sum(new_idx))
   print('Number of started incidents: ',sum(started))

   #old incidents for tracking
   #add newest started incidents to this lists
   for inc in incidents:
      if inc.started and inc.incident_id_string not in started_inc_ids:
         started_inc_ids.append(inc.incident_id_string)
   csv.started_inc_ids = started_inc_ids

   #update the ngfs_day object
   #maybe some of these could be called inside each other??
   csv.add_data(data) #only csv rows with an incident name or id
   csv.add_full_data(full_data) #all csv data rows
   csv.add_incidents(incidents)

   #red_flag information
   csv.red_flag_warnings = rf_warnings
   csv.red_flag_zones = rf_zones


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
   #don't save pickle file if no new incidents are in the data
   if sum(csv.new):
      csv.save_pickle()
   csv.print_base_map()
