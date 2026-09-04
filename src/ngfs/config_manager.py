"""
Loads the configuration the NGFS system runs from.

load_cfgs reads two layers:

  * etc/ngfs.json -- NGFS-specific settings: data sources and ingest
    directories, run limits, region overrides, unnamed-fire handling. See
    README section 5 for the keys.
  * the wrfxpy configuration it names in 'wrfxpy_cfg', normally etc/conf.json.

Both paths are relative, so the process working directory decides which
installation is configured.

make_base_configuration supplies the base wrfxpy job description that every
incident forecast is derived from, read from jobs/base_ngfs_cfg.json. If that
file is missing it falls through to simple_forecast's interactive
questionnaire, which will hang a cron run -- keep it in place.

get_fmda_path is imported from ngfs_incident at call time rather than at module
scope, so this module stays free of that module's heavy imports and no import
cycle can form. persistence.get_old_incidents uses the same deferred-import
pattern.
"""
from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
import json
import utils
import simple_forecast as sf
import time



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
         #imported here rather than at module scope so this module does not depend
         #on ngfs_incident's heavy imports, and so no import cycle can form; the
         #same deferred-import pattern is used by persistence.get_old_incidents
         from ngfs.ngfs_incident import get_fmda_path
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
      

