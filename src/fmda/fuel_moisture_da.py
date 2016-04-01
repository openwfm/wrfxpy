# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from wrf.wrf_data import WRFModelData
from trend_surface_model import fit_tsm
from utils import great_circle_distance, find_closest_grid_point
from fuel_moisture_model import FuelMoistureModel
from fm10_observation import FM10Observation
from utils import find_closest_grid_point

import sys
import os
import numpy as np
from datetime import datetime, timedelta
import pytz
import netCDF4
import logging
import json
from MesoPy import Meso

def check_overlap(wrf_path,ts_now):
  """
  Check if the WRF file <wrf_path> timstamps contain <ts_now>.
  """
  wrfout = WRFModelData(wrf_path)
  outts = wrfout['GMT']
  if ts_now in outts:
    return True
  else:
    print("INFO: previous forecast [%s - %s] exists, running DA till %s" % (str(outts[0]),str(outts[-1]),str(ts_now)))
    return False


def retrieve_mesowest_observations(meso_token, tm_start, tm_end, wrf_model):
    """
    Retrieve observation data from Mesowest and repackage them as a time-indexed
    dictionary of lists of observations.

    :param meso_obss: the structure returned from the Mesowest query

    """
    def decode_meso_time(t):
        # example: '2016-03-30T00:30:00Z'
        return datetime.strptime(t, '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=pytz.UTC)

    def meso_time(dt):
        # example: 201603311600
        return '%04d%02d%02d%02d%02d' % (dt.year, dt.month, dt.day, dt.hour, dt.minute)

    # retrieve data from Mesowest API (http://mesowest.org/api/)
    m = Meso(meso_token)
    meso_obss = m.timeseries(meso_time(tm_start - timedelta(minutes=30)),
                          meso_time(tm_end + timedelta(minutes=30)),
                          showemptystations = '0', bbox='%g,%g,%g,%g' % wrf_model.get_domain_extent(),
                          vars='fuel_moisture')


    glat, glon = wrf_model.get_lats(), wrf_model.get_lons()
   
    # repackage all the observations into a time-indexed structure which groups
    # observations at the same time together
    obs_data = {}
    for stinfo in meso_obss['STATION']:
        st_lat, st_lon = float(stinfo['LATITUDE']), float(stinfo['LONGITUDE'])
        elev = float(stinfo['ELEVATION']) / 3.2808
        ngp = find_closest_grid_point(st_lon, st_lat, glon, glat)
        dts = [decode_meso_time(x) for x in stinfo['OBSERVATIONS']['date_time']]
        if 'fuel_moisture_set_1' in stinfo['OBSERVATIONS']:
            fms = stinfo['OBSERVATIONS']['fuel_moisture_set_1']

            for ts,fm_obs in zip(dts,fms):
                o = FM10Observation(ts,st_lat,st_lon,elev,float(fm_obs),ngp)
                obs_t = obs_data.get(ts, [])
                obs_t.append(o)
                obs_data[ts] = obs_t

    return obs_data


def execute_da_step(model, model_time, covariates, fm10):
  """
  Execute a single DA step from the current state/extended parameters and covariance matrix using
  the <covariates> and observations <fm10>.  Assimilation time window is fixed at 60 mins.

  :param model: a FuelMoistureModel
  :param model_time: the current model time
  :param covariates: the covariate fields to take into account to model the spatial structure of the FM field
  :param fm10: the 10-hr fuel moisture observations
  """
  valid_times = [z for z in fm10.keys() if abs((z - model_time).total_seconds()) < 1800]
  logging.info('FMDA there are %d valid times at model time %s' % (len(valid_times), str(model_time)))

  if len(valid_times) > 0:

    # retrieve all observations for current time
    obs_valid_now = []
    for z in valid_times:
        obs_valid_now.extend(fm10[z])

    fmc_gc = model.get_state()
    dom_shape = fmc_gc.shape[:2]

    # construct covariate storage
    Xd3 = min(len(covariates) + 1, len(obs_valid_now))
    X = np.zeros((dom_shape[0], dom_shape[1], Xd3))
    X[:,:,0] = fmc_gc[:,:,1]
    for i,c in zip(range(Xd3-1),covariates):
        X[:,:,i+1] = covariates[i]

    # run the trend surface model (clamp output to [0.0 - 2.5] to be safe)
    Kf_fn, Vf_fn = fit_tsm(obs_valid_now, X)
    Kf_fn[Kf_fn < 0.0] = 0.0
    Kf_fn[Kf_fn > 2.5] = 2.5

    Kg = np.zeros((dom_shape[0], dom_shape[1], 6))

    # run the data assimilation step now
    logging.info("Mean Kf: %g Vf: %g state[0]: %g state[1]: %g state[2]: %g\n" %
      (np.mean(Kf_fn), np.mean(Vf_fn), np.mean(fmc_gc[:,:,0]), np.mean(fmc_gc[:,:,1]), np.mean(fmc_gc[:,:,2])))
    model.kalman_update_single2(Kf_fn[:,:,np.newaxis], Vf_fn[:,:,np.newaxis,np.newaxis], 1, Kg)
    logging.info("Mean Kf: %g Vf: %g state[0]: %g state[1]: %g state[2]: %g\n" %
      (np.mean(Kf_fn), np.mean(Vf_fn), np.mean(fmc_gc[:,:,0]), np.mean(fmc_gc[:,:,1]), np.mean(fmc_gc[:,:,2])))


def run_data_assimilation(wrf_model, fm10, wrf_model_prev = None):
    """
    Run the fuel moisture and DA for all time steps in the model wrf_model.
    If a previous run is available, the fuel moisture values (and covariance if available)
    are transferred.

    :param wrf_model: the current WRF data file to process (wrf input or wrf output)
    :param fm10: a list of the observations of 10-hr fuel moisture available
    :param wrf_model_prev: optional, the previous WRF data file from which fm state may be copied
    :return: the fuel moisture model with the assimilated fields
    """
    tss = wrf_model.get_gmt_times()
    lat, lon = wrf_model.get_lats(), wrf_model.get_lons()
    dom_shape = lat.shape
    T2 = wrf_model['T2']
    Q2 = wrf_model['Q2']
    PSFC = wrf_model['PSFC']
    hgt = wrf_model['HGT']
    rain = wrf_model['RAIN']
    rain = np.log(rain + 1.0)
    constant = np.ones_like(T2)
    Ed,Ew = wrf_model['Ed'], wrf_model['Ew']
    E = 0.5 * (Ed[0,:,:] + Ew[0,:,:])
   
    P0 = np.diag([0.01,0.01,0.01,0.01,0.001,0.001])
    Tk = np.array([1.0, 10.0, 100.0, 1000.0]) * 3600
    Q = np.diag([1e-4,5e-5,1e-5,1e-6,1e-6,1e-6])

    # initialize the grid moisture model with the fuel moisture equilibrium
    model = FuelMoistureModel(E[:,:,np.newaxis][:,:,np.zeros((4,),dtype=np.int)], Tk, P0)

    # if a previous fuel moisture model run is available, copy it's state
    if wrf_model_prev is not None:
        logging.info('FMDA replacing fuel moisture equilibrium with previous calculation from %s.' % wrf_model_prev.path)
        prev_tss = wrf_model_prev.get_gmt_times()
        if tss[0] in prev_tss:
            prev_ndx = prev_tss.index(tss[0])
            model.get_state()[:,:,:3] = wrf_model_prev['FMC_GC'][prev_ndx,:3,:,:].transpose((1,2,0))
            model.get_state()[:,:,3:5] = wrf_model_prev['FMEP'][prev_ndx,:,:,:].transpose((1,2,0))
        
    # precompute static covariates (we assume domains don't move around)
    cov_lon = lon - np.mean(lon)
    cov_lat = lat - np.mean(lat)
    cov_hgt = hgt / 1000.0
    cov_const = np.ones(dom_shape)

    # advance model and run DA for each timestep
    for i, ts in enumerate(tss):
        cov_t2 = T2[i,:,:]
        cov_q2 = Q2[i,:,:]
        cov_psfc = PSFC[i,:,:]
        cov_rain = np.log(rain[i,:,:] + 1.0)
        covariates = [cov_t2, cov_psfc,cov_lon,cov_lat,cov_hgt,cov_t2,cov_q2,cov_const]
        if np.any(rain > 0.0):
            covariates.append(rain)

        if i > 0:
            model.advance_model(Ed[i,:,:], Ew[i,:,:], rain[i,:,:], (ts - tss[i-1]).seconds, Q)

        execute_da_step(model, ts, covariates, fm10)

        # overwrite the WRF model variables for this time step
        d = netCDF4.Dataset(wrf_model.path, 'r+')
        d.variables['FMC_GC'][i,:3,:,:] = model.get_state()[:,:,:3].transpose(2,0,1)
        d.variables['FMEP'][i,:,:,:] = model.get_state()[:,:,4:6].transpose(2,0,1)
        d.close()


def assimilate_fm10_observations(path_wrf, path_wrf0, mesowest_token):
  
  # load the wrfinput file
  wrfin = WRFModelData(path_wrf, ['T2', 'Q2', 'PSFC', 'HGT', 'FMC_GC', 'FMEP'])
  lat, lon = wrfin.get_lats(), wrfin.get_lons()
  tss = wrfin.get_gmt_times()
  tm_start, tm_end = tss[0], tss[-1]
  dom_shape = lat.shape
  logging.info('FMDA domain size is %d x %d grid points with extent (%g to %g) lons (%g to %g)' %
              (dom_shape[0], dom_shape[1],np.amin(lat),np.amax(lat),np.amin(lon),np.amax(lon)))

  # compute the diagonal distance between grid points
  grid_dist_km = great_circle_distance(lon[0,0], lat[0,0], lon[1,1], lat[1,1])
 
  # retrieve fuel moisture observations via the Mesowest API
  fm10 = retrieve_mesowest_observations(mesowest_token, tm_start, tm_end, wrfin)

  logging.info('FMDA retrieved %d observations from Mesowest.' % len(fm10))

  # if a previous cycle is available (i.e. the wrfoutput is a valid file), load the model
  prev_wrf = None
  if path_wrf0 is not None and os.path.exists(path_wrf0) and check_overlap(path_wrf0,tm_start):
      prev_wrf = WRFModelData(path_wrf0)
      outts = prev_wrf['GMT']
      logging.info("FMDA previous forecast [%s - %s] exists" % (str(outts[0]),str(outts[-1])))
  else:
      logging.info("FMDA no previous forecast found, running DA from equilibrium at %s" % str(tm_start))
 
  # run from the start until now (retrieve fuel moisture, extended parameters, covariance matrix)
  run_data_assimilation(wrfin, fm10, prev_wrf)

  return 0


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('usage: %s <wrf-file> <mesowest_token>' % sys.argv[0])
        sys.exit(1)

    #cfg = json.load(open(sys.argv[1]))
    assimilate_fm10_observations(sys.argv[1], None, sys.argv[2])
    sys.exit(0)

