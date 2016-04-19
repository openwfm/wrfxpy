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

from fmda.fuel_moisture_da import execute_da_step, retrieve_mesowest_observations
from fmda.fuel_moisture_model import FuelMoistureModel
from ingest.grib_file import GribFile, GribMessage
from ingest.rtma_source import RTMA
from utils import Dict, ensure_dir, utc_to_esmf
from vis.postprocessor import scalar_field_to_raster
from ssh_shuttle import send_product_to_server


import netCDF4
import numpy as np
import json
import sys
import logging
import os
import os.path as osp

from datetime import datetime, timedelta
import pytz


def postprocess_cycle(cycle, region_cfg, wksp_path):
    """
    Build rasters from the computed fuel moisture.

    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :return: the postprocessing path
    """
    data_path = compute_model_path(cycle, region_cfg.code, wksp_path)
    year_month = '%04d%02d' % (cycle.year, cycle.month)
    cycle_dir = 'fmda-%s-%04d%02d%02d-%02d' %  (region_cfg.code, cycle.year, cycle.month, cycle.day, cycle.hour)
    postproc_path = osp.join(wksp_path, year_month, cycle_dir)

    # open and read in the fuel moisture values
    d = netCDF4.Dataset(data_path)
    fmc_gc = d.variables['FMC_GC'][:,:,:]
    d.close()

    # read in the longitudes and latitudes
    geo_path = osp.join(wksp_path, '%s-geo.nc' % region_cfg.code)
    d = netCDF4.Dataset(geo_path)
    lats = d.variables['XLAT'][:,:]
    lons = d.variables['XLONG'][:,:]
    d.close()

    fm_wisdom = {
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'jet_r',
       'scale' : [0.0, 0.4]
    }

    esmf_cycle = utc_to_esmf(cycle) 
    mf = { "1" : {esmf_cycle : {}}}
    manifest_name = 'fmda-%s-%04d%02d%02d-%02d.json' %  (region_cfg.code, cycle.year, cycle.month, cycle.day, cycle.hour)
    ensure_dir(osp.join(postproc_path, manifest_name))

    for i,name in [(0, '1-hr'), (1, '10-hr'), (2, '100-hr')]:
        fm_wisdom['name'] = '%s fuel moisture' % name
        raster_png, coords, cb_png = scalar_field_to_raster(fmc_gc[:,:,i], lats, lons, fm_wisdom)
        raster_name = 'fmda-%s-raster.png' % name
        cb_name = 'fmda-%s-raster-cb.png' % name
        with open(osp.join(postproc_path, raster_name), 'w') as f:
            f.write(raster_png)
        with open(osp.join(postproc_path, cb_name), 'w') as f:
            f.write(cb_png) 
        mf["1"][esmf_cycle][name] = { 'raster' : raster_name, 'coords' : coords, 'colorbar' : cb_name }
        json.dump(mf, open(osp.join(postproc_path, manifest_name), 'w'))

    return postproc_path


def compute_model_path(cycle, region_code, wksp_path):
    """
    Construct a relative path to the fuel moisture model file
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    year_month = '%04d%02d' % (cycle.year, cycle.month)
    filename = 'fmda-%s-%04d%02d%02d-%02d.nc' %  (region_code, cycle.year, cycle.month, cycle.day, cycle.hour)
    return osp.join(wksp_path,region_code,year_month,filename) 



def find_region_indices(glat,glon,minlat,maxlat,minlon,maxlon):
    """
    Find the indices i1:i2 (lat dimension) and j1:j2 (lon dimension)
    that contain the desired region (minlat-maxlat,minlon-maxlon).

    :param glat: the grid latitudes
    :param glon: the grid longitudes
    :param minlat: the minimum latitude
    :param maxlat: the maximum latitude
    :param minlon: the minimum longitude
    :param maxlon: the maximum longitude
    :return: dim 0 min/max indices and dim1 min/max indices
    """
    i1, i2, j1, j2 = 0, glat.shape[0], 0, glat.shape[1]
    done = False
    while not done:
        done = True
        tmp = np.where(np.amax(glat[:, j1:j2],axis=1) < minlat)[0][-1]
        if i1 != tmp:
            i1 = tmp
            done = False
        tmp = np.where(np.amin(glat[:, j1:j2],axis=1) > maxlat)[0][0]
        if i2 != tmp:
            i2 = tmp
            done = False
        tmp = np.where(np.amax(glon[i1:i2,:],axis=0) < minlon)[0][-1]
        if j1 != tmp:
            j1 = tmp
            done = False
        tmp = np.where(np.amin(glon[i1:i2,:],axis=0) > maxlon)[0][0]
        if j2 != tmp:
            j2 = tmp
            done = False
    return i1,i2,j1,j2



def load_rtma_data(rtma_data, bbox):
    """
    Load relevant RTMA fields and return them
    
    :param rtma_data: a dictionary mapping variable names to local paths
    :param bbox: the bounding box of the data
    :return: a tuple containing t2, rh, lats, lons
    """
    gf = GribFile(rtma_data['temp'])[1]
    lats, lons = gf.latlons()
    
    # bbox format: minlat, minlon, maxlat, maxlon
    i1, i2, j1, j2 = find_region_indices(lats, lons, bbox[0], bbox[2], bbox[1], bbox[3])
    
    t2 = gf.values()[i1:i2,j1:j2] # temperature at 2m in K
    td = GribFile(rtma_data['td'])[1].values()[i1:i2,j1:j2] # dew point in K
    rain = GribFile(rtma_data['precipa'])[1].values()[i1:i2,j1:j2] # precipitation
    hgt = GribFile('static/ds.terrainh.bin')[1].values()[i1:i2,j1:j2]
    
    # compute relative humidity
    rh = 100*np.exp(17.625*243.04*(td - t2) / (243.04 + t2 - 273.15) / (243.0 + td - 273.15))
    
    return t2, rh, rain, hgt, lats[i1:i2,j1:j2], lons[i1:i2,j1:j2]


def compute_equilibria(T, H):
    """
    Compute the drying and wetting equilibrium given temperature and relative humidity.
    
    :param T: the temperature at 2 meters in K
    :param H: the relative humidity in percent
    :return: a tuple containing the drying and wetting equilibrium
    """
    d = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    w = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    d *= 0.01
    w *= 0.01
    return d, w


def fmda_advance_region(cycle, cfg, rtma, wksp_path, lookback_length, meso_token):
    """
    Advance the fuel moisture estimates in the region specified by the configuration.
    The function assumes that the fuel moisture model has not been advanced to this
    cycle yet and will overwrite any previous computations.
    
    Control flow:
    
    1) read in RTMA variables
    2) check if there is a stored FM model for previous cycle
    2a) yes -> load it, advance one time-step, perform DA
    2b) no -> compute equilibrium, use background covariance to do DA
    3) store model
    
    :param cycle: the datetime indicating the processed cycle in UTC
    :param cfg: the configuration dictionary specifying the region
    :param rtma: the RTMA object that can be used to retrieve variables for this cycle
    :param wksp_path: the workspace path for the cycler
    :param lookback_length: number of cycles to search before we find a computed cycle
    :param meso_token: the mesowest API access token
    :return: the model advanced and assimilated at the current cycle
    """
    model = None
    prev_cycle = cycle - timedelta(hours=1)
    prev_model_path = compute_model_path(prev_cycle, cfg.code, wksp_path)
    if not osp.exists(prev_model_path):
        logging.info('CYCLER cannot find model from previous cycle %s' % str(prev_cycle))
        if lookback_length > 0:
            model = fmda_advance_region(cycle - timedelta(hours=1), cfg, rtma, wksp_path, lookback_length - 1, meso_token)
    else:
        logging.info('CYCLER found previous model for cycle %s.' % str(prev_cycle))
        model = FuelMoistureModel.from_netcdf(prev_model_path)
        
    # retrieve the variables and make sure they are available (we should not be here if they are not)
    dont_have_vars, have_vars = rtma.retrieve_rtma(cycle)
    assert not dont_have_vars
    
    logging.info('CYCLER loading RTMA data for cycle %s.' % str(cycle))
    T2, RH, rain, hgt, lats, lons = load_rtma_data(have_vars, cfg.bbox)
    Ed, Ew = compute_equilibria(T2, RH)

    dom_shape = T2.shape

    # store the lons/lats for this domain
    geo_path = osp.join(wksp_path, '%s-geo.nc' % cfg.code)
    if not osp.isfile(geo_path):
        d = netCDF4.Dataset(geo_path, 'w', format='NETCDF4')
        d.createDimension('south_north', dom_shape[0])
        d.createDimension('west_east', dom_shape[1])
        xlat = d.createVariable('XLAT', 'f4', ('south_north', 'west_east'))
        xlat[:,:] = lats
        xlong = d.createVariable('XLONG', 'f4', ('south_north', 'west_east'))
        xlong[:,:] = lons
        d.close()
    
    
    # the process noise matrix
    Q = np.diag([1e-4,5e-5,1e-5,1e-6,1e-6])
    
    # background covariance
    P0 = np.diag([0.01,0.01,0.01,0.001,0.001])

    # check if we must start from equilibrium
    if model is None:
        logging.info('CYCLER initializing from equilibrium for cycle %s.' % (str(cycle)))
        # setup model parameters    
        Nk = 3
        Tk = np.array([1.0, 10.0, 100.0])
        m0 = np.expand_dims(0.5 * (Ed + Ew), axis=2)
        model = FuelMoistureModel(m0[:,:,[0,0,0]], Tk, P0)
    else:
        logging.info('CYCLER advancing model one hour to cycle %s.' % (str(cycle)))
        dt = 3600 # always 1 hr step in RTMA
        model.advance_model(Ed, Ew, rain, dt, Q)
    
    logging.info('CYCLER retrieving fm-10 observations for cycle %s.' % (str(cycle)))
    
    # perform assimilation with mesowest observations
    tm_start = cycle - timedelta(minutes=30)
    tm_end = cycle + timedelta(minutes=30)
    fm10 = retrieve_mesowest_observations(meso_token, tm_start, tm_end, lats, lons)
    fm10v = []
    for fm10_obs in fm10.values():
        for obs in fm10_obs:
            fm10v.append(obs.get_value())
    
    logging.info('CYCLER retrieved %d valid observations, min/mean/max [%g/%g/%g].' %
                 (len(fm10),np.amin(fm10v),np.mean(fm10v),np.amax(fm10v)))
    
    # remove rain that is too small to make any difference
    rain[rain < 0.01] = 0
    
    # run the data assimilation step
    covs = [np.ones(dom_shape), hgt / 2000.0]
    if np.any(rain > 0.01):
        covs.append(rain)
    execute_da_step(model, cycle, covs, fm10)
    
    # store the new model
    model.to_netcdf(ensure_dir(compute_model_path(cycle, cfg.code, wksp_path)))
    
    return model
    
    
def is_cycle_computed(cycle, cfg, wksp_path):
    """
    Check if the fuel model file exists (has been computed) for the
    cycle <cycle> and region configuration <cfg>.
    
    :param cycle: the cycle datetime in UTC
    :param cfg: the region configuration wrapped in a Dict for convenience
    :param wksp_path: the workspace path for the cycler
    :return: True if the model file has been found, False otherwise
    """
    path = compute_model_path(cycle, cfg.code, wksp_path)
    return osp.isfile(path)
    
    
if __name__ == '__main__':
    
    # setup environment
    sys_cfg = Dict(json.load(open('etc/conf.json')))
    cfg = Dict(json.load(open('etc/rtma_cycler.json')))
    meso_token = json.load(open('etc/tokens.json'))['mesowest']
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # current time
    now = datetime.now(pytz.UTC)
    cycle = (now - timedelta(minutes=50)).replace(minute=0,second=0,microsecond=0)
    logging.info('CYCLER activated at %s, will attempt cycle at %s' % (str(now), str(cycle)))
    
    # what is the most recent RTMA data available?
    lookback_length = cfg.lookback_length
    dont_have_vars, have_vars = None, None
    rtma = RTMA('ingest', ['precipa', 'wspd', 'wdir', 'td', 'temp'])
    while lookback_length > 0:
        dont_have_vars, have_vars = rtma.retrieve_rtma(cycle)
        if dont_have_vars:
            logging.info('RTMA variables %s not yet available for cycle %s.' % (str(dont_have_vars), str(cycle)))
            cycle -= timedelta(hours=1)
            lookback_length -= 1
        else:
            break
            
    if dont_have_vars:
        logging.error('CYCLER could not find useable cycle, exiting.')
        sys.exit(1)
        
    logging.info('Have RTMA data for cycle %s.' % str(cycle))
      
    # check for each region, if we are up to date w.r.t. RTMA data available
    for region_id,region_cfg in cfg.regions.iteritems():
        wrapped_cfg = Dict(region_cfg)
        if not is_cycle_computed(cycle, wrapped_cfg, cfg.workspace_path):
            logging.info('CYCLER processing region %s for cycle %s' % (region_id, str(cycle)))
            fmda_advance_region(cycle, wrapped_cfg, rtma, cfg.workspace_path, lookback_length, meso_token)
            pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path)   
            if 'shuttle_remote_host' in sys_cfg:
                sim_code = 'fmda-' + wrapped_cfg.code
                send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, region_id + ' FM')
        else:
            logging.info('CYCLER the cycle %s for region %s is already complete, skipping ...' % (str(cycle), str(region_id)))

    # done
    logging.info('CYCLER cycle %s complete.' % str(cycle))
