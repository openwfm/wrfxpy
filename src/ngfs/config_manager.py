#makes configuration for incidents and ngfs setup
from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
import json
import utils
import simple_forecast as sf
import time
from ngfs_start import make_geo_folder



#loading the configuration for NGFS and WRFXPY
def load_cfgs():
   ### ngfs configuration
   if os.path.exists('etc/ngfs.json'):
      print('Reading the ngfs configuration file')
      with open('etc/ngfs.json','r') as openfile:
            ngfs_cfg = json.load(openfile)
            #should walk user through this if not found  <<<-----------------------------------------------------
   else:
      print('The ngfs configuration file was not found')
   ### wrfxpy configuration
   print('Reading the base wrfxpy configuration file')
   with open(ngfs_cfg['wrfxpy_cfg'],'r') as openfile:
         wrfxpy_cfg = json.load(openfile)
   return ngfs_cfg, wrfxpy_cfg
   
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
      print(f'{num_starts} simulations will be started automatically. Make sure you have the resources.')
      time.sleep(2)
   else:
      print('Simulations will need to be started manually')
    
   return auto_start, num_starts, force
   
def make_base_configuration(force,ngfs_cfg):
   #check to see if a  base ngfs configuration file exists
   if utils.file_exists('jobs/base_ngfs_cfg.json'):
      print('Configuration found')
      with open('jobs/base_ngfs_cfg.json','r') as openfile:
         temp_cfg = json.load(openfile)
      #json_print = json.dumps(temp_cfg, indent=4)   
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
   
def get_fmda_path(cfg):
   #
   fmda_year = cfg['start_utc'][:4]
   fmda_month = cfg['start_utc'][5:7]
   fmda_day = cfg['start_utc'][8:10]
   fmda_hour = cfg['start_utc'][11:13]
   #
   if 'cawfe' in cfg['fire_namelist_path']:
      base_folder = f'/data/WRFXPY/wksp_fmda/CONUS/{fmda_year}{fmda_month}/'
   elif 'behave' in cfg['fire_namelist_path']:
      base_folder = f'/data/jhaley/wrfxpy/wksp_fmda/CONUS/{fmda_year}{fmda_month}/'
   else:
      return 'ngfs'
   #
   date_folder = f'fmda-CONUS-{fmda_year}{fmda_month}{fmda_day}-{fmda_hour}.geo'
   #
   fmda_geo_folder = base_folder + date_folder
   #
   if not os.path.exists(fmda_geo_folder):
      make_geo_folder(fmda_geo_folder)
   #
   return fmda_geo_folder

if __name__ == '__main__':
   print('Testsing the configuration manager')
   try:
      ngfs_cfg, wrfxpy_cfg = load_cfgs() #in __main__ fore testimg
      #print_cfg(ngfs_cfg)
      print('NGFS config')
      print(json.dumps(ngfs_cfg, indent=4))
      print('WRFXP config')
      print(json.dumps(wrfxpy_cfg, indent=4))
   except:
      print('Error loading the configurations')
      

