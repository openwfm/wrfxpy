#ngfs_incident class
from __future__ import absolute_import
from __future__ import print_function
import geopandas as gpd
import numpy as np
import pandas as pd
import pickle
import copy
import csv
import os, sys, glob
import subprocess
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
from ngfs import constants as cons
from ingest.downloader import download_url
from fmda.fuel_moisture_model import FuelMoistureModel 
from ngfs import ngfs_api as api


#from ngfs_dictionary import ngfs_dictionary  <-- use to force csv data columns into a type
#import ngfs_helper as nh
import utils, logging, traceback
#dictionary to convert between "CA" and "California", etc
import state_names as sn
import simple_forecast as sf
import ngfs_dictionary as nd
#import shapely
from shapely.geometry import Point, LineString, Polygon
import urban_incident as ui
#depending on the pyproj available
from shapely.geometry import Polygon, Point, box
try:
   from pyproj import Proj, transform, Transformer
except:
   from pyproj import Proj, transform # Transformer
from ingest.NWS_warnings import set_red_flags, red_flag_incident, subset_red_flag_data, subset_wfo_data

def get_fmda_path(cfg):                                                                #<<--------------move elsewhere
   fmda_year = cfg['start_utc'][:4]
   fmda_month = cfg['start_utc'][5:7]
   fmda_day = cfg['start_utc'][8:10]
   fmda_hour = cfg['start_utc'][11:13]

   if 'cawfe' in cfg['fire_namelist_path']:
      base_folder = '/data/WRFXPY/wksp_fmda/CONUS/' + fmda_year + fmda_month + '/'
   elif 'behave' in cfg['fire_namelist_path']:
      base_folder = '/data/jhaley/wrfxpy/wksp_fmda/CONUS/' + fmda_year + fmda_month + '/'
   else:
      return 'ngfs'

   date_folder = 'fmda-CONUS-'+fmda_year+fmda_month+fmda_day+'-'+fmda_hour+'.geo'

   fmda_geo_folder = base_folder + date_folder

   if not os.path.exists(fmda_geo_folder):
      make_geo_folder(fmda_geo_folder)

   return fmda_geo_folder

def make_geo_folder(fmda_geo_folder):                                                  #<<--------------move elsewhere
   #copies netcdf file from older fmda installation into something that can be used by behave model
   #create the folder 
   print(f'\tMaking {fmda_geo_folder} and writing geogrid files')
   os.makedirs(fmda_geo_folder,exist_ok=True)
   #location of netcdf file
   nc_path = fmda_geo_folder.replace('jhaley/wrfxpy','WRFXPY').replace('.geo','.nc')
   if not os.path.exists(nc_path):
      print('FMDA netcdf file not found, returning')
      return
   #load the model 
   fm = FuelMoistureModel.from_netcdf(nc_path)
   # Extend the array to have 6 dimensions
   m_ext = np.zeros(fm.m_ext.shape[:2] + (6, ))
   m_ext[:,:,:-1] = fm.m_ext
   fm.m_ext = m_ext
   #load index
   index = {'projection': 'lambert',
            'dx' : 2539.703,
            'dy' : -2539.703,
            'truelat1' : 25.0,
            'truelat2' : 25.0,
            'stdlon' : 265,
            'radius' : 6371200.0,
            'known_x': 1072.0,
            'known_y': 629.0,
            'known_lat': 39.12699765672619,
            'known_lon': -95.48481787184572
            }
   # Save into Geogrid format
   fm.to_geogrid(fmda_geo_folder, index)

def parse_landcover_string(landcover_str):
   """
   Parse a landcover string into a dictionary with
   landcover types as keys and fractions as values.
   format is like the following:
      'Trees:96,Shrubs:2,Grass/Herbs:1,Water:1'
   """
   #sometimes this may contain a NaN
   if ':' not in str(landcover_str):
      return pd.DataFrame()
   data = {}
   for item in landcover_str.split(","):
      landcover, percent = item.split(":")
      data[landcover.strip()] = float(percent) / 100.0
   return pd.DataFrame(data,index = [0])

def print_landcover(lc):
   for k in lc.keys():
      print(f'\t\t{k}\t{np.round(lc[k][0],3)}')

def incident_landcover(data,verbose = False):
   """
   Determines landcover and fuel types for an incident based on NGFS version 3.8+ csv
   """
   if 'land_cover' not in data.keys():
      print('Older NGFS version does not have landcover types')
      return None,None
   #empty datframes
   land_cover = pd.DataFrame()
   fuel = pd.DataFrame()
   #loop through the detection data
   for i in data.index:
      land_cover = pd.concat([land_cover,parse_landcover_string(data.land_cover.loc[i])]).fillna(0)
      fuel = pd.concat([fuel,parse_landcover_string(data.fuel.loc[i])]).fillna(0)
   #collapse to a one-row datframe, made with the averages which sum to 1.0
   lc = land_cover.mean().to_frame().T
   fl = fuel.mean().to_frame().T
   if verbose:
      print('\tLand cover:')
      print_landcover(lc)
      print('\tFuels:')
      print_landcover(fl)
   return lc,fl

def estimate_viirs_time(goes_data,viirs_pixel):
      #searches goes data for the closest and earliest detection center and assings a time to the viirs pixel
      unique_data = goes_data.drop_duplicates(subset=['latitude','longitude'],keep='first').sort_values(by='acq_date_time').reset_index().copy()
      v_lon = viirs_pixel['longitude']
      v_lat = viirs_pixel['latitude']
      min_dist = 20 #starting distance is 20 degrees off
      t_est = viirs_pixel['acq_date']
      for i in unique_data.index:
         g_lon = unique_data['longitude'].loc[i]
         g_lat = unique_data['latitude'].loc[i]
         d = np.sqrt((v_lon-g_lon)**2 + (v_lat-g_lat)**2)
         if d < min_dist:
            t_est = unique_data['acq_date_time'].loc[i]
            min_dist = d
      return min(t_est,viirs_pixel['acq_date'])

def minimum_perim(data):
   #find a small perim around a set of detections points
      data = data.drop_duplicates(subset=['latitude','longitude'])
      lat_str = 'latitude_cINDEX'
      lon_str = 'longitude_cINDEX'
      #column name of NGFS data
      lat_columns = [lat_str.replace('INDEX',str(i)) for i in range(1, 5)]
      lon_columns = [lon_str.replace('INDEX',str(i)) for i in range(1, 5)]
      lon_center = data.longitude.mean()
      lat_center = data.latitude.mean()
      min_perim = []
      for j in data.index:
         if not 'NGFS' in data.loc[j,'version']:  #take the center of a FIRMS pixel
            min_lat = data.loc[j,'latitude']
            min_lon = data.loc[j,'longitude']
         else:  #take the nearest corner to the center of a NGFS pixel
            min_dist = 20
            for i in range(4):
               lat = data.loc[j,lat_columns[i]]
               lon = data.loc[j,lon_columns[i]]
               d = np.sqrt((lon_center-lon)**2 + (lat_center-lat)**2)
               if d < min_dist:
                  min_dist = d
                  min_lat = lat
                  min_lon = lon
         min_perim.append([min_lat,min_lon])
      return min_perim


class ngfs_incident():
   '''
   The ngfs_incident class provides functionality to organize and configure forecasts for data associated with a named incident
   '''
   def __init__(self,name=None,data=None,base_cfg = None,ngfs_cfg = None):
      self.name = name
      self.incident_name = None
      self.incident_id_string = name
      self.cfg = None #for this incident
      self.base_cfg = base_cfg #default cfg
      self.ngfs_cfg = ngfs_cfg #cfg of ngfs day system
      self.job_filename = None
      self.log_filename = None
      #detection data
      if not data is None:
         self.data = data
         self.incident_name = str(data['known_incident_name'].unique()[0])
         self.incident_id_string = data['known_incident_id'].unique()[0]
      else:
         self.incident_name = None
         self.incident_id_string = name
         self.data = pd.DataFrame()
      self.viirs_data = pd.DataFrame()
      self.viirs_only = False
      self.viirs_ignition_pixel = False
      self.feature_data = pd.DataFrame()
      self.det_count = None
      self.ignition_pixel = None
      self.feature_tracking_id = list()
      self.bbox = None
      self.total_frp = None
      self.unique_latlon = None
      self.unique_latlons = None
      #forecast/incident status
      self.new = False
      self.started = False
      self.continuing = False
      if '.' in self.incident_id_string: #incident id string containing '.' will be for an unknown fire
         self.unknown = True
         print('This is an unknown possible wildland fire')
         #print(f'Loading {ngfs_cfg['unknown']['unknown_cfg_file']} as configuration file')
         #load special configuration
         with open(ngfs_cfg['unknown']['unknown_cfg_file'],'r') as openfile:
            self.base_cfg = json.load(openfile) 
      else:
         self.unknown = False
         self.base_cfg = base_cfg #default cfg
      #incident location
      self.county = None
      self.state = None
      self.loc_str = None
      self.affected_population = 0
      #forcast configuration
      self.force = False
      self.auto_start = False
      #incident ancillary data
      self.red_flag = False
      self.RX = False
      self.urban = None
      self.land_cover = None
      self.fuel = None
      self.nifc_perims = []
      #set current time as estimate ignition time
      self.ign_utc = pd.Timestamp.now(tz='UTC')
      self.start_utc = self.ign_utc - timedelta(hours = 2)
      self.end_utc = self.ign_utc + timedelta(hours = 4)
      self.ign_latlon = [40.0,-120.0]
   '''
   def __setstate__(self,state):
      self.__dict__ = state
      defaults = {
         'cfg':None,
         'base_cfg':None,
         'ngfs_cfg':None,
         'job_filename':None,
         'log_filename':None,
         'viirs_data':None,
         'viirs_only':None,
         'viirs_ignition_pixel':None,
         'feature_data':None,
         'feature_tracking_id':None,
         'continuing':None,
         'unknown':None,
         'red_flag':None,
         'RX':None,
         'urban':None,
         'land_cover':None,
         'fuel':None,
         'nifc_perims':None,
         'json_start_code':None
         }
      for key, value in defaults.items():
         if key not in self.__dict__:
            setattr(self,key,value)

      def merge_incidents(inc1,inc2)

   '''



   def __eq__(self,other):
      #for comparison with other 
      if not isinstance(other,ngfs_incident):
         return NotImplemented
      return self.incident_id_string == other.incident_id_string

   def add_nifc_perims(self):
      if not self.ngfs_cfg is None:
         perim_directory = self.ngfs_cfg['perims_cfg']['perim_dir'] # can be a list
         p = glob.glob(f'{perim_directory}*{self.incident_id_string}*')
         if len(p) > 0:
            print('Found matching perimeters')
            self.nifc_perims.extend(p)  

   def set_RX(self):
      #'known_incident_type':'incident_type' v2,v1 keys
      if 'RX' in self.incident_name.upper(): # or 'RX' in self.data.incident_type.unique():
         self.RX = True

   #these are the columns of the csv file as pandas dataframe
   def add_data(self,df):
      self.data =  pd.concat([self.data,df],ignore_index=True)
      self.data = self.data.drop_duplicates()
      if not hasattr(self.data,'longitude_b5'):
         self.viirs_only = True
      else:
         self.viirs_only = False
   def add_viirs_data(self,viirs_df):
      self.viirs_data =  pd.concat([self.viirs_data,viirs_df],ignore_index=True)
      self.viirs_data = self.viirs_data.drop_duplicates()

   def set_incident_start_time(self):
      if len(self.data) > 0:
         #self.incident_start_time = self.data.iloc[0,3]
         self.incident_start_time = min(self.data['acq_date_time'])
      else:
         self.incident_start_time = 'Unset start time'

   def set_json_start_code(self):
      self.json_start_code = f'./forecast.sh {self.job_filename} &> {self.log_filename} &'

   def start_forecast(self,sleep_time = 1):
      print('Starting forecast')
      cmd = f'sleep {sleep_time}; {self.json_start_code}'
      #subprocess.Popen(cmd,shell=True)
      #os.system(f'(sleep {sleep_time}; {self.json_start_code})')   ### <<<-----------------------------bad? remove ?
      print(f'Starting job after {sleep_time} delay: {self.json_start_code}')
      self.started = True

   
   # These attributes are set when viirs data finds detection pixels within the goes ignition pixel
   def set_new_ign_latlon(self,new_ign_latlon):
      if not any(np.isnan(new_ign_latlon)):
         self.new_ign_latlon = new_ign_latlon
      elif any(np.isnan(new_ign_latlon)) and hasattr(self,'new_ign_latlon'):
         del(self.new_ign_latlon)
         pass
         #print('Not setting new_ign_latlon')
   def set_new_ign_utc(self,new_ign_utc):
      self.new_ign_utc = new_ign_utc
      if new_ign_utc < self.ign_utc:
         print('\tChanging the time parameters of the job for ealier detection information')
         print('\tNew UTC time:',new_ign_utc)
         print('\tOld UTC time:',self.ign_utc)
         forecast_length = self.end_utc - self.start_utc
         self.ign_utc = new_ign_utc
         #the start and finish times of the simulation get adjusted too
         self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=60))
         self.end_utc = self.start_utc + forecast_length
   def set_incident_bounding_box(self):
      #returns bounding box
      try: #ngfs version 2+
         lat_columns = [f'latitude_c{i}' for i in range(1, 5)]
         lon_columns = [f'longitude_c{i}' for i in range(1, 5)]
         min_lat = min(min(self.data[col]) for col in lat_columns)
         max_lat = max(max(self.data[col]) for col in lat_columns)
         min_lon = min(min(self.data[col]) for col in lon_columns)
         max_lon = max(max(self.data[col]) for col in lon_columns)
      except:
         #ngfs_version 1
         lat_columns = [f'lat_c{i}' for i in range(1, 5)]
         lon_columns = [f'lon_c{i}' for i in range(1, 5)]
         min_lat = min(min(self.data[col]) for col in lat_columns)
         max_lat = max(max(self.data[col]) for col in lat_columns)
         min_lon = min(min(self.data[col]) for col in lon_columns)
         max_lon = max(max(self.data[col]) for col in lon_columns)
      self.bbox = min_lon, max_lon, min_lat, max_lat
      print('\tBounding box:', self.bbox)

   def expand_domain_size(self):
      #check to see if domain will cover the the detection footprint
      #make a geoseries  of the ingition point and construct a buffer around it
      #see if buffer contains the bounding box, increas the size of the domain if not
      print('\tChecking for domain size fit')
      def make_detection_buffer(x,y,r):
      #make a buffer region arounf the ignition pixel
         #convert the lon,lat to meter units
         if 'Transformer' not in dir():
            inProj = Proj(init='epsg:4326')
            outProj = Proj(init='epsg:3857')
            xp,yp = transform(inProj,outProj,x,y)
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
      try:
         x = self.data['longitude'].mean()
         y = self.data['latitude'].mean()
      except:
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
         sz = self.cfg['domains']['1']['domain_size'][0] - 1
         new_sz = sz*2**size_increase + 1
         self.cfg['domains']['1']['domain_size'] = [new_sz,new_sz]
         sr = self.cfg['domains']['1']['subgrid_ratio'][0]
         new_sr = int(sr*(1/2)**(size_increase-1))
         self.cfg['domains']['1']['subgrid_ratio'] = [new_sr,new_sr]
         ppn = self.cfg['ppn']
         self.cfg['ppn'] = min(400,ppn*2**(size_increase*2))
         print('\tResizing domain, doublings = ',size_increase)
         #center the domain on the maen of all data
         domain_lat = self.data.latitude.mean()
         domain_lon = self.data.longitude.mean()
         print(f'\tNew domain center after expansion: {domain_lat},{domain_lon}')
         self.cfg['domains']['1']['truelats'] = (domain_lat, domain_lat)
         self.cfg['domains']['1']['center_latlon'] = [domain_lat,domain_lon]     #<<<<<--------------------------------------    Center the domain on the data, not ignition point
         self.cfg['fire_namelist_path'] = 'etc/nlists/default.fire_cawfe_13'
         self.cfg['domains']['1']['time_step'] = 4

      #if size_increase > 2:
      #   print('Domain too large to forecast')
   def random_pix_loc(self,data_pixel):
      #returns a randondom [lat,lon] pair for a detection pixel 
      if ('latitude_b5' in data_pixel.keys()) and (abs(data_pixel['latitude_b5']) < 90): #GOES pixel with SWIR estimate
         lat_str = 'latitude_cINDEX_b5'
         lon_str = 'longitude_cINDEX_b5'
      elif 'longitude_c1'  in data_pixel.keys(): #NGFS pixel with corners
         lat_str = 'latitude_cINDEX'
         lon_str = 'longitude_cINDEX'
      else: #FIRMS ixel without boundary, take center
         return [data_pixel['latitude'],data_pixel['longitude']]
      #colums name of NGFS data
      lat_columns = [lat_str.replace('INDEX',str(i)) for i in range(1, 5)]
      lon_columns = [lon_str.replace('INDEX',str(i)) for i in range(1, 5)]
      #bounds on the pixel area
      min_lat = min([data_pixel[col] for col in lat_columns])
      max_lat = max([data_pixel[col] for col in lat_columns])
      min_lon = min([data_pixel[col] for col in lon_columns])
      max_lon = max([data_pixel[col] for col in lon_columns])

      #choose uniforn random location in boundary
      random_lat = np.random.uniform(min_lat,max_lat)
      random_lon = np.random.uniform(min_lon,max_lon)

      return [random_lat,random_lon]
   
   def forecast_points(self):
      #uses fire_init to make ignitions at all of the points in hotspots points in the domain
      #find unique lat-lon points in data
      #unique_data = self.data.drop_duplicates(subset=['latitude','longitude'],keep='first').copy()
      viirs_sats = ['NOAA-20', 'NOAA-21', 'SNPP']
      if len(self.viirs_data) > 10:
         unique_data = self.viirs_data
      else:   
         unique_data = pd.concat([self.data,self.viirs_data])
      ignitions = []
      for i in unique_data.index:
         data = unique_data.loc[i].copy()
         if 'NGFS' in data['version']:
            points_to_add = 5
         else:
            points_to_add = 1
         for p in range(points_to_add):
            latlon = self.random_pix_loc(data) # [data.latitude,data.longitude]
            if data['satellite'] in viirs_sats:
               data['acq_date_time'] = estimate_viirs_time(self.data,data)
            time_utc =utils.utc_to_esmf(data.acq_date_time)
            d = {
               "latlon" : latlon,
               "time_utc" : time_utc,
               "duration_s" : 1200,
               "radius": 100,
               "ros" : 1
            }
            ignitions.append(d)
      self.cfg['ignitions']['1'] = ignitions
      self.cfg['use_tign_ignition'] = True
      self.cfg['burn_plot_boundary'] = []

   



   def make_incident_configuration(self,base_cfg,ngfs_cfg):
      #change the base configuration file to match the individual incidents parameters
      #print(self.incident_id_string)
      #print(self.data.keys())
      cfg = copy.deepcopy(base_cfg)
      #select the burn model via namelist.fire
      if 'fire_namelist_path' in ngfs_cfg['run_cfg'].keys():
         cfg['fire_namelist_path'] = ngfs_cfg['run_cfg']['fire_namelist_path']
         print('\tUsing ',cfg['fire_namelist_path'], 'for namelist.fire')
         if 'behave_13' in cfg['fire_namelist_path']:
            cfg['domains']['1']['time_step'] = 6

      if 'region_cfg_REMOVE_THIS' in ngfs_cfg.keys():       # <<<< ------------------------------------------------    FIX
         for r in ngfs_cfg['region_cfg']:
            print(r)
            if self.data.state.unique() in ngfs_cfg['region_cfg'][r]['state']:
               cfg['grib_source'] = ngfs_cfg['region_cfg'][r]['grib_source']
               cfg['geo_vars_path'] = ngfs_cfg['region_cfg'][r]['geo_vars_path']
               print(ngfs_cfg['region_cfg'][r]['msg'])
      else:
         update_states = ['CA','AZ','NV','UT', 'NM'] #['OR','WA','ID','MT','WY','CO']
         data_states = self.data.state.unique()
         update_list = [us for us in data_states if us in update_states]
         if len(update_list) > 0:
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
         if any(self.data.state == 'PR'):
            cfg['grib_source'] = 'GFSF'
            print('\tPuerto Rico incident detected, using GFSF and PRVI Landfire data')
            cfg['geo_vars_path'] = 'etc/vtables/geo_vars.json_prvi'
         if any(self.data.state == 'VI'):
            cfg['grib_source'] = 'GFSF'
            print('\tVirgin Islands incident detected, using GFSF and PRVI Landfire data')
            cfg['geo_vars_path'] = 'etc/vtables/geo_vars.json_prvi'
            #maybe use the old adrjrw here because there is so much ocean in the domain?
            #cfg['wrf_namelist_path']  = "etc/nlists/default.input_adjrw"
            #cfg['fire_namelist_path'] = "etc/nlists/default.fire_adjrw"
      try:
         print('\tGrib source: ',cfg['grib_source'])
      except:
         print('\tGrib source is unset')
      #remove the { and } characters at the end of the incident id strings
      gc_string = f'{self.incident_name}_{utils.utc_to_esmf(self.start_utc)}_{self.incident_id_string[1:-1]}'
      #replace weird characters that cause problems with operating system file handling
      replace_chars = ['#','(',')',':',' ']
      for rc in replace_chars:
         gc_string = gc_string.replace(rc,'_')
      cfg['grid_code'] = gc_string

      #time and place of the ignition
      ignitions = cfg['ignitions']
      ign_dur = ignitions['1'][0]['duration_s']
      #change if there is a viirs detection to work with
      #the domain stays the same otherwise so comparision between viirs and goes ignitions may be examined
      if hasattr(self,'new_ign_latlon') and not any(np.isnan(self.new_ign_latlon)):
         print(f'Setting new_ign_latlon {self.new_ign_latlon}')
         t_utc = min(self.ign_utc,self.new_ign_utc)
         cfg['ignitions'] = { '1' : [ { 'time_utc' : utils.utc_to_esmf(t_utc),
                                   'duration_s' : ign_dur,
                                   'latlon' : self.new_ign_latlon } ] }  # <<---- new
         cfg['domains']['1']['center_latlon'] = self.new_ign_latlon
         cfg['domains']['1']['truelats'] = (self.new_ign_latlon[0], self.new_ign_latlon[0])
         cfg['domains']['1']['stand_lon'] = self.new_ign_latlon[1]  
      else:
         cfg['ignitions'] = { '1' : [ { 'time_utc' : utils.utc_to_esmf(self.ign_utc),
                                    'duration_s' : ign_dur,
                                    'latlon' : self.ign_latlon } ] }
         cfg['domains']['1']['center_latlon'] = self.ign_latlon
         cfg['domains']['1']['truelats'] = (self.ign_latlon[0], self.ign_latlon[0])
         cfg['domains']['1']['stand_lon'] = self.ign_latlon[1]
      
      start_utc = utils.esmf_to_utc(base_cfg['start_utc']) #esmf_to_utc
      end_utc = utils.esmf_to_utc(base_cfg['end_utc'])
      forecast_length = end_utc-start_utc

      cfg['start_utc'] = utils.utc_to_esmf(self.start_utc)
      cfg['end_utc'] = utils.utc_to_esmf(utils.round_time_to_hour(self.ign_utc + forecast_length + timedelta(hours =  + 1.5))) #utils.round_time_to_hour(self.ign_utc - timedelta(minutes=30))

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
      if self.ngfs_cfg['fmda_cfg']['use_fmda']:
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
      else:
         print('\tWill use eqilibrium FMC')

      #print(cfg)
      #self.incident_name
      cfg['postproc']['description'] = self.incident_name
      self.job_filename = 'jobs/' + cfg['grid_code'] + '.json'
      self.log_filename = self.job_filename.replace('jobs/','logs/').replace('.json','.log')
      self.cfg = cfg
      #check to see if domain will cover the the detection footprint
      #make a geoseries  of the ingition point and construct a buffer around it
      #see if buffer contains the bounding box, increas the size of the domain if not
      self.expand_domain_size()
      self.set_json_start_code()
      json.dump(self.cfg, open(self.job_filename, 'w'), indent=4, separators=(',', ': '))
      
      del cfg

      #viewprint(cfg)


   def find_old_features(self,goes_data):
      #looks through the GOES detections for earlier occurences of the features tracking ID before the incident was named
      #the feature trakcing ID may not be unique
      time_key = self.data_time_key()
      self.feature_tracking_id = []
      print(self.incident_name)
      for sat in goes_data['satellite'].unique():
         sat_subset = self.data[self.data['satellite'] == sat]
         sat_features = list(sat_subset['feature_tracking_id'].unique())
         print(f'\tIncident has {len(sat_features)} feature tracking ids in {sat}: {sat_features}')
         self.feature_tracking_id.extend(sat_features)
         # Search for older feature tracking IDs in full data, add it to data if found
         for f in sat_features:
            #data with correct feature tracjking id for the specified satellite
            feature_subset = goes_data[(goes_data.feature_tracking_id == f) & (goes_data.satellite == sat)]
            #filter out data that should already part of the incident subet
            feature_subset = feature_subset.drop(feature_subset[feature_subset.known_incident_id == self.incident_id_string].index)
            #only interested in earlier detections
            feature_subset = feature_subset[feature_subset[time_key] < min(self.data[time_key])]
            if not feature_subset.empty:
               print(f'\tFound {len(feature_subset)} earlier hotspots with feature tracking id: {f}')
               print(f'\tMean lat / lon: {np.mean(feature_subset.latitude)} / {np.mean(feature_subset.longitude)}')
               self.feature_data = self.feature_data.append(feature_subset, ignore_index=True)
               #make sure feature_tracking data is in the right place
               feature_subset = feature_subset[
                     (feature_subset.longitude >= min(self.data.longitude)) &
                     (feature_subset.longitude <= max(self.data.longitude)) &
                     (feature_subset.latitude >= min(self.data.latitude)) &
                     (feature_subset.latitude <= max(self.data.latitude))
               ]
               if not feature_subset.empty:
                     feature_subset['known_incident_id'] = self.incident_id_string
                     feature_subset['known_incident_name'] = self.incident_name
                     self.data = pd.concat([self.data,feature_subset],ignore_index=True)
                     self.data = self.data.sort_values(by=time_key)

   def find_old_detections(self,goes_data):
      #finds other detections in the same location(s)   <<< ---------------------------- may be unnecessary, feature tracking id should have this
      unique_latlons = np.asarray(self.data[['latitude', 'longitude']].drop_duplicates())
      added_data = pd.DataFrame()
      for lat,lon in unique_latlons:
         try:
            full_subset = goes_data[
               (goes_data['longitude'] == lat) &
               (goes_data['latitude'] == lon)
            ]
         except KeyError:    #  <<<<------------------------------------------------------- probably can remove this
            full_subset = goes_data[
               (goes_data['lon'] == lat) &
               (goes_data['lat'] == lon)
            ]
         if len(full_subset) > 0:
            added_data = pd.concat([full_subset,added_data],ignore_index = True)
      if len(added_data) > len(self.data):
         time_key = self.data_time_key()
         print('\tFound additional, earlier detections')
         self.data = pd.concat([self.data,added_data],ignore_index=True)
         self.data = self.data.drop_duplicates()
         self.data = self.data.sort_values(by=time_key)

   def data_time_key(self):
      if 'acq_date_time' in self.data.keys():
         return 'acq_date_time'
      for k in self.data.keys():
         if 'time' in k:
            return k

   def set_ignition_point(self):
      print('\tFinding start time and ignition location')
      self.set_incident_start_time()   
      print(f'\tEstimated incident start time: {self.incident_start_time}')
      time_key = self.data_time_key()
      #print(f'Working with time key : {time_key}')
      idx = self.data.index[0]  #should be earliest
      min_time = self.data[time_key].min()
      #print(self.data)
      #print(f'Min time: {min_time}')
      #print(f'Min time type{type(min_time)}')
      try:
         tc_ign_latlon = [self.data.latitude[idx], self.data.longitude[idx]]  #first pixel of the data
         #look at detections within first three hours
         time_msk = (self.data[time_key] < min_time + timedelta(hours=3.0))   # <<<<< ---------------------------- maybe a better way to do this?
         #-999 used as fill values when longitude_b5 not determined
         if hasattr(self.data,'longitude_b5'):
            swir_msk = abs(self.data.longitude_b5) < 180.0 
            self.viirs_only = False
         else:
            #working with viirs data only
            self.viirs_only = True
            swir_msk = abs(self.data.longitude) < 180.0 #should evaluate True for entire array
         #average of nominal pixel locations
         mean_ign_latlon = [np.mean(self.data.latitude[time_msk]), np.mean(self.data.longitude[time_msk])]
         #print(f'mean_ign_latlon {mean_ign_latlon}')
         if sum(time_msk & swir_msk) > 1 and not self.viirs_only:
            swir_ign_latlon = [np.mean(self.data.latitude_b5[time_msk & swir_msk]), np.mean(self.data.longitude_b5[time_msk & swir_msk])]
            print('\tAveraging SWIR pixels from first three hours')
            print(f'\tMean SWIR ign_latlon: {swir_ign_latlon}')
         elif sum(time_msk & swir_msk) > 1 and self.viirs_only:
            swir_ign_latlon = mean_ign_latlon
            print('\tUsing average of VIIRS ignition points from first three hours')
            print(f'\tMean VIIRS ign_latlon: {swir_ign_latlon}')
         else:
            swir_ign_latlon = [6000,3000]  # not on the map
         print(f'\tFirst detection ign_latlon: {tc_ign_latlon}')
         print(f'\tMean nominal pixel ign_latlon: {mean_ign_latlon}')
         #make sure the swir and nominal locations aren't too far apart
         lat_diff = abs(mean_ign_latlon[0] - swir_ign_latlon[0])
         lon_diff = abs(mean_ign_latlon[1] - swir_ign_latlon[1])
         self.ign_latlon = swir_ign_latlon if max(lat_diff, lon_diff) < 0.04 else mean_ign_latlon
      except KeyError:      # #  <<<<------------------------------------------------------- probably can remove this
         print('\tNo terrain corrected lat/lon available')
         self.ign_latlon = [self.data.lat[idx], self.data.lon[idx]]
         mean_ign_latlon = [np.mean(self.data.latitude), np.mean(self.data.longitude)]
         self.ign_latlon = mean_ign_latlon
      print('\tUsing mean ignition point')
      self.unique_latlon = np.unique(self.ign_latlon)
      #print(self.data.acq_date_time[idx],self.incident_start_time,self.data.pixel_date_time[idx])
      self.ign_utc = min(
         self.data.acq_date_time[idx],
         self.incident_start_time,
         self.data.pixel_date_time[idx]
      )

   def NGFS_red_flag(self):
    #looks through a DataFrame with NGFS detection data fro NWS fire weather data
    #df is data associated with an incident
    df = self.data
    keys = ['nws_fire_wx_code', 'event_type']
    for k in keys:
        if k in df.keys():
            codes = df[k].unique()
            for c in codes:
                if (int(c) == 3 or int(c)) == 4:
                    print('\tDetections indicate fire weather')
                    return True
    return False

   def set_incident_demographics(self,pop_data):
      #finds location, population data for the incident
      self.county = self.data['county'].iloc[0]
      self.state = self.data['state'].iloc[0]
      #print(self.state,self.county) 
      if len(self.state) == 2:
         try:
            self.state = sn.abbrev_to_us_state[self.state]
            self.loc_str = f'{self.county}, {self.state}'
         except KeyError:
            print('Error translating the state abbreviation(s)')
            self.state = 'Unknown'
            self.loc_str = 'Unknown'
      #
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
   
   def process_incident(self, goes_data, viirs_data = pd.DataFrame(), rf_zones = [], pop_data = pd.DataFrame()):
      """
      Processes an incident by analyzing the incident data, extracting unique feature tracking IDs,
      updating detection data, estimating start time, and checking for red flag warnings.

      goes_data should be a full set of GOES data going back 24 hours+ before earliest detection in incident subset
      """
      #look though older data for detections associated with incident
      self.find_old_features(goes_data=goes_data)
      #self.find_old_detections(goes_data = goes_data) <<<<<<------------------------- remove this?
      # Sort data by observation time
      time_key = self.data_time_key()
      self.data = self.data.sort_values(by=time_key)   # <<<<< --------------------------------  already sorted?
      
      #some stats and ancillary data about the detections
      self.det_count = len(self.data)
      print(f'\tNumber of detections: {self.det_count}')
      self.unique_latlons = np.asarray(self.data[['latitude', 'longitude']].drop_duplicates())
      print(f'\tNumber of unique detection locations: {len(self.unique_latlons)}')
      self.set_incident_bounding_box()
      try:
         self.total_frp = self.data['frp'].sum()
      except KeyError:
         self.total_frp = self.data['total_frp'].sum()
      #first detection pickel
      self.ignition_pixel = self.data.iloc[0]

      self.set_incident_demographics(pop_data)
      self.set_RX()
      self.land_cover,self.fuel = incident_landcover(self.data,verbose=True)
      #self.red_flag = red_flag_incident(self.ign_latlon[1], self.ign_latlon[0], self.ign_utc, rf_zones, data = self.data)   # <<<<< possible remove since NGFS provides fire weather flags
      self.red_flag = self.NGFS_red_flag()
      if self.red_flag:
         print('\tIncident in a red flag warning zone')

      #find time and place of ignition
      self.set_ignition_point()
   
      #refine  ingition point estimate, if possible
      self.polar_ign_estimate(viirs_data = viirs_data)

      now = pd.Timestamp.now('UTC')
      self.new = now - self.ign_utc < timedelta(hours = 24)
      if not self.new:
         print('\tIncident older than 24 hours')
      
      #configuration of the incident, possible move into make_configuration
      self.start_utc = utils.round_time_to_hour(self.ign_utc - timedelta(minutes=60))
      self.end_utc = self.start_utc + timedelta(hours=24)
      self.time_utc = utils.utc_to_esmf(self.ign_utc)
      print(f'\tLocation of earliest GOES pixel in csv: {self.ign_latlon}')
      print(f'\tUTC ignition time: {self.time_utc}')

      #used for creating string used for system call to start simulation   <<<------------------- unused?
      def set_cmd_str(self,cmd_str):
         self.cmd_str = cmd_str

   def ngfs_viirs_ignition(self,viirs_data = pd.DataFrame()):
      #returns data set with VIIRS data from NGFS detections
      print('\tLooking for VIIRS ignition pixels within NGFS data')
      viirs_ign_data = viirs_data[viirs_data['known_incident_id'] == self.incident_id_string]
      if len(viirs_ign_data) > 0:
         print('\t',len(viirs_ign_data),' viirs ignition polygon detections from NGFS data')
         return viirs_ign_data
      else:
         print('\tNo NGFS VIIRS detections for the incident')
         return pd.DataFrame()

   def firms_viirs_ignition(self,viirs_data = pd.DataFrame()):
      #use shapely to help find the first viirs pixel within the boundaries of the goes pixel
      #self.inc_bb is the incident bounding box (minx, miny, maxx, maxy)
      #setup ignition pixel polygon
      #print(self.ignition_pixel)
      print('\tLooking for VIIRS ignition pixels withing NASA FIRMS data')
      try: #NGFS Version 2.x+
         vertices = [(self.ignition_pixel[f'longitude_c{i}'], self.ignition_pixel[f'latitude_c{i}']) for i in range(1, 5)]
      except: #NGFS Version 1.x
         vertices = [(self.ignition_pixel[f'lon_tc_c{i}'], self.ignition_pixel[f'lat_tc_c{i}']) for i in range(1, 5)]
      ignition_polygon = Polygon(vertices)

      #setup incident polygon, easier construction for square box
      incident_polygon = box(self.bbox[0],self.bbox[2],self.bbox[1],self.bbox[3])
      
      #brute force way to get at it, loop through list of all viirs detections
      #viirs detections within first all GOES pixels boundary
      viirs_inc_data = pd.DataFrame()
      #viirs detections within first GOES ignition pixel
      viirs_ign_data = pd.DataFrame()

      for i in viirs_data.index:
         try:
               #point object
               viirs_point = Point(viirs_data.loc[i,'longitude'],viirs_data.loc[i,'latitude'])
         except:
               print('Error finding viirs point')
               viirs_point = Point(600,3000) # <<------- not in any polygon
         
         if viirs_point.within(incident_polygon):
               if viirs_point.within(ignition_polygon):
                  viirs_ign_data = viirs_ign_data.append(viirs_data.loc[i])
                  #print('\tFound viirs pixel within earliest GOES pixel:',viirs_ign,viirs_data.loc[i,'acq_date'])
               else:
                  #print('\tFound additional incident pixel')
                  viirs_inc_data = viirs_inc_data.append(viirs_data.loc[i])
                  #print('length of viirs INC pixels',len(viirs_inc_data))

      print('\t',len(viirs_inc_data),'viirs incident detections, ',len(viirs_ign_data),' viirs ignition polygon detections from FIRMS data')
      return viirs_inc_data, viirs_ign_data

   def polar_ign_estimate(self,viirs_data = pd.DataFrame()):
      #use shapely to help find the first viirs pixel within the boundaries of the goes pixel
      #self.inc_bb is the incident bounding box (minx, miny, maxx, maxy)
      #GOignition_pixe
      print('\tFinding best viirs pixel for the incident GOES ignition point')
      #print(data)
      try:
         print('\tTime of GOES imaging: ',self.ignition_pixel.loc['observation_time'])
      except:
         print('\tTime of GOES imaging: ',self.ignition_pixel.loc['pixel_date_time'])
      
      #try to find VIIRS pixels within NGFS data
      viirs_ign_data = self.ngfs_viirs_ignition(viirs_data=viirs_data)
      if len(viirs_ign_data) > 0:
         viirs_inc_data = viirs_ign_data.copy()
         ngfs_viirs = True
      else:
         viirs_inc_data, viirs_ign_data = self.firms_viirs_ignition(viirs_data=viirs_data)
         ngfs_viirs = True
      
      '''
      #assign a state of the pixel location
      if len(viirs_ign_data) > 0:
         viirs_ign_data.loc[viirs_ign_data.index,'in_goes_pixel'] = True
      if len(viirs_inc_data) > 0:
         viirs_inc_data.loc[viirs_inc_data.index,'in_goes_pixel'] = False
      '''

      #three possibilities VIIRS data within first GOES pixel, within foorprint of all GOELS pixel, no viirs data at all
      ig_frame = pd.DataFrame()
      if (len(viirs_ign_data) == 0 and len(viirs_inc_data) > 0):
         print('\tUsing earliest VIIRS pixel, but it is not within the earliest GOES Pixel')
         #viirs_ign_data = viirs_inc_data
         print('\t',len(viirs_inc_data),' new incident points found')
         ig_frame = viirs_inc_data
         viirs_inc_data.loc[viirs_inc_data.index,'in_goes_pixel'] = False
      #ignition point from first NGFS data or goes pixel location
      if len(viirs_ign_data) > 0:
         print('\t',len(viirs_ign_data),' new ignition pixel points found')
         ig_frame = viirs_ign_data
         viirs_ign_data.loc[viirs_ign_data.index,'in_goes_pixel'] = True
         if len(viirs_inc_data) > 0:
            if ngfs_viirs:
               viirs_inc_data.loc[viirs_inc_data.index,'in_goes_pixel'] = True
            else:
               viirs_inc_data.loc[viirs_inc_data.index,'in_goes_pixel'] = False

      if len(ig_frame) > 0:
         #sort data by time
         min_ign_time = ig_frame['acq_date'].min()    # <<<--------------------------  not a native field for NGFS data, probably should add the field to firms data instead
         #boolean index
         first_pixels = ig_frame[ig_frame['acq_date'] <= min_ign_time + timedelta(hours=3)]
         print('\tCount of earliest viirs pixels: ',len(first_pixels))
         #replace with better estimate of lat/lon pf ignition
         mean_viirs_lat = first_pixels['latitude'].mean()
         mean_viirs_lon = first_pixels['longitude'].mean()
         # return these to be added to the ngfs_incidents
         new_ign_latlon = [mean_viirs_lat,mean_viirs_lon]
         new_ign_utc = min_ign_time
         print('\tBest viirs ignition estimate location: \n\t\t',new_ign_latlon,new_ign_utc)
         self.viirs_ignition_pixel = True
      else:
         new_ign_latlon = [float('NaN'),float('NaN')]
         new_ign_utc = pd.NaT  # << -------- return something of correct data type+
      
      self.set_new_ign_latlon(new_ign_latlon)
      self.set_new_ign_utc(new_ign_utc)
      #join the incidenet and ingnition data

      #join and sort the two datframes
      viirs_ign_data = pd.concat([viirs_ign_data,viirs_inc_data],ignore_index=True)
      if len(viirs_ign_data) > 0:
         viirs_ign_data = viirs_ign_data.sort_values(by='acq_date',ignore_index=True)
      else:
         print('\tNo VIIRS data found for incident')

      self.add_viirs_data(viirs_ign_data)

   def print_incident(self):
      print(f'{self.incident_id_string}')

      print(f'\t{self.incident_name}')
      print(f'\t{self.loc_str}')
      print(f'\t\t{self.ign_latlon} {self.ign_utc}')
      print(f'\tNumber of GOES detections: {len(self.data)}')
      print(f'\tNumber of VIIRS detections: {len(self.viirs_data)}')
      print(f'\tStarted: {self.started}')
      print()

   def initial_forecast(self):
      from nifc import get_nifc_incident
      
      from sat_webpage import sat_overpass, sat_dict
      try: 
         lat = self.ign_latlon[0]
         lon = self.ign_latlon[1]
      except:
         lat = self.data.latitude.mean()
         lon = self.data.longitude.mean()

      #get viirs times for future overpasses, this function just prins for now
      sat_overpass(lat = lat, lon = lon, alt = 500, horizon = 15, length =8)
      '''
      except:
         print('Error with satellite overpasses')
      '''
      #changes the incident_id_string and runs a forecast when there are no VIIRS detections
      print('Making an initial forecast')
      #change id string and grid_code
      old_id = self.incident_id_string
      new_id = old_id[:-13]+'InitialFcast}'
      self.incident_id_string = new_id
      #rename grid_code
      self.cfg['grid_code'] = self.cfg['grid_code'][:-12] + 'InitialFcast'
      self.job_filename = 'jobs/' + self.cfg['grid_code'] + '.json'
      self.log_filename = self.job_filename.replace('jobs/','logs/').replace('.json','.log')
      self.set_json_start_code()
      print(f'Changing id from {old_id} to {new_id}')
      #change forecast duration
      start_time = utils.esmf_to_utc(self.cfg['start_utc'])
      new_end_time = start_time + timedelta(hours = 8)
      self.cfg['end_utc'] = utils.utc_to_esmf(new_end_time)
      print('Will try to find NIFC location for ignition point')
      try:
         #id = '{'+old_id+'}'
         loc = get_nifc_incident(old_id,feature_type='loc')
         print('Found nifc location', loc)
         if len(loc['features']) > 0:
            #print('nifc api',loc['features'][0]['properties'].keys())
            lon,lat = loc['features'][0]['geometry']['coordinates']
            nifc_ign_latlon = [lat,lon]
         else:
            nifc_ign_latlon = None
      except:
         print('Error getting nifc location data')
         nifc_ign_latlon = None

      if not nifc_ign_latlon:
         #use old ignition point already computed
         print('Will use best estimate from GOES data')
      else:
         print('Setting new ignition point')
         self.cfg['ignitions']['1'][0]['latlon'] = nifc_ign_latlon

      #save the job file
      json.dump(self.cfg, open(self.job_filename, 'w'), indent=4, separators=(',', ': '))

      '''
      cfg['ignitions'] = { '1' : [ { 'time_utc' : utils.utc_to_esmf(t_utc),
                                   'duration_s' : ign_dur,
                                   'latlon' : self.new_ign_latlon } ] } 
      '''
         



if __name__ == '__main__':
    pass
