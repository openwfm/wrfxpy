import f90nml
import os.path as pth
from wrf_cloner import WRFCloner

def execute(job_id, wksp_path, wps_idir, wrf_idir, grib_source,
            wps_nml_path, wrf_nml_path, fire_nml_path):
  """
  Executes a weather/fire simulation. **args are in a dictionary that is keyed
  by the following names.

  Args:
    job_id - job of ID
    wksp_path - workspace path
  """
  # read in all namelists
  wps_nml = f90nml.read(wps_nml_path)
  wrf_nml = f90nml.read(wrf_nml_path)
  fire_nml = f90nml.read(fire_nml_path)

  # build directories in workspace
  wps_dir = pth.join([wksp_path, job_id, 'wps'])
  wrf_dir = pth.join([wksp_path, job_id, 'wrf'])

  # step 1: clone WPS and WRF directories
  cln = WRFCloner(wrf_idir = wrf_idir, wps_idir = wps_idir)
  cln.clone_wps(wps_dir, grib_source.vtables(), [])
  cln.clone_wrf(wrf_dir, [])

  # step 2: patch namelist for geogrid and execute geogrid
  wps_nml['geogrid']['geog_data_path'] = args['geogrid_path']
  f90nml.write(wps_nml, pth.join(wps_dir, 'namelist.wps'))

  # step 3: download GRIB files from the grib_source that aren't locally available

  # step 4: make GRIBFILE.YYY symlinks, patch namelist for ungrib end execute ungrib

  # step 5: execute metgrid

  # step 6: execute real.exe

  # step 7: execute wrf.exe on parallel backend


  

def test_firejob():
  from grib_source import GribSource

  args = { 'job_id' : 'test-job',
           'workspace_dir' : 'wksp',
           'wps_install_dir' : '/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/WPS',
           'wrf_install_dir' : '/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/WRFV3',
           'grib_source' : GribSource(),
           'wps_namelist' : 'etc/nlists/firejob-colo.wps',
           'wrf_namelist' : 'etc/nlists/firejob-colo.wrf',
           'fire_namelist' : 'etc/nlists/firejob-colo.fire' }
  
  firejob(**args)

