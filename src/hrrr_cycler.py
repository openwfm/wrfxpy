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
from ingest.grib_file import GribFile
from ingest.HRRRA import HRRRA
from ingest.HRRR import HRRR
from utils import Dict, ensure_dir, utc_to_esmf, delete, force_copy
from vis.postprocessor import scalar_field_to_raster, scatter_to_raster
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

# setup environment
sys_cfg = Dict(json.load(open('etc/conf.json')))
cfg = Dict(json.load(open('etc/fmda_cycler.json')))
meso_token = json.load(open('etc/tokens.json'))['mesowest']
hrrr_vars = {
    'rain': 635, 'snow': 636, 't2': 616, 'rh': 620, 'psfc': 607, 
    'soil_t_0': 560, 'soil_t_1': 562, 'soil_t_2': 564, 
    'soil_t_3': 566, 'soil_t_4': 568, 'soil_t_5': 570, 
    'soil_t_6': 572, 'soil_t_7': 574, 'soil_t_8': 576, 
    'soil_moist_0': 561, 'soil_moist_1': 563, 'soil_moist_2': 565, 
    'soil_moist_3': 567, 'soil_moist_4': 569, 'soil_moist_5': 571, 
    'soil_moist_6': 573, 'soil_moist_7': 575, 'soil_moist_8': 577,
    'snowh': 615, 'spfh': 618, 'u10': 622, 'v10': 623, 'ws': 624,
    'ust': 643, 'znt': 642, 'swdown': 664 
}

def write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels=None, alpha=None):
    """
    Write postprocessing files.

    :param post: the UTC cycle time
    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :return: the postprocessing path
    """
    raster_name = cycle_dir + '-%s-raster.png' % name
    cb_name = cycle_dir + '-%s-raster-cb.png' % name
    with open(osp.join(postproc_path, raster_name), 'wb') as f:
        f.write(raster_png)
    with open(osp.join(postproc_path, cb_name), 'wb') as f:
        f.write(cb_png) 
    mf["1"][esmf_cycle][name] = { 'raster' : raster_name, 'coords' : coords, 'colorbar': cb_name }
    if levels is not None:
        mf["1"][esmf_cycle][name].update({ 'levels' : levels })
    if alpha is not None:
        mf["1"][esmf_cycle][name].update({ 'alpha' : alpha })


def postprocess_cycle(cycle, region_cfg, wksp_path, fcst_hour, bounds=None):
    """
    Build rasters from the computed fuel moisture.

    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :param bounds: bounding box of the post-processing
    :return: the postprocessing path
    """
    model_path = compute_model_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    cycle_id = compute_fmda_id(cycle, region_cfg.code)
    postproc_path = compute_postproc_path(cycle, region_cfg.code, wksp_path, fcst_hour)
    if fcst_hour > 0:
        cycle_id += f'-f{fcst_hour:02d}'
        prev_cycle = cycle
        prev_fcst_hour = fcst_hour - 1
    else:
        prev_cycle = cycle - timedelta(hours=1)
        prev_fcst_hour = fcst_hour 
    prev_postproc_path = compute_postproc_path(prev_cycle, region_cfg.code, wksp_path, prev_fcst_hour)
    if prev_fcst_hour > 0:
        prev_cycle_id = compute_fmda_id(prev_cycle, region_cfg.code) + f'-f{prev_fcst_hour:02d}'
    else:
        prev_cycle_id = compute_fmda_id(prev_cycle, region_cfg.code)
    
    manifest_name = cycle_id + '.json'
    complete_manifest_name = f'fmda-{region_cfg.code}.json'
    if not is_cycle_computed(cycle, region_cfg, wksp_path, fcst_hour) and not osp.exists(prev_postproc_path):
        logging.warning(f'CYCLER postprocessing misses information for cycle {cycle}')
        return None

    var_wisdom = {
        'dfm' : {
            'native_unit' : '-',
            'colorbar' : '-',
            'colormap' : 'jet_r',
            'scale' : [0.0, 0.4]
        },
        'lfm' : {
            'native_unit' : '-',
            'colorbar' : '-',
            'colormap' : 'jet_r',
            'scale' : [0.0, 3.0],
            'marker' : '^'
        },
        'EQUILd FM' : {
            'name' : 'Drying equilibrium FM',
            'native_unit' : '-',
            'colorbar' : 'i-',
            'colormap' : 'jet_r',
            'scale' : [0.0, 0.4]
        },
        'EQUILw FM' : {
            'name' : 'Wetting equilibrium FM',
            'native_unit' : '-',
            'colorbar' : 'i-',
            'colormap' : 'jet_r',
            'scale' : [0.0, 0.4]
        },
        'RH' : {
            'name' : 'Relative humidity',
            'native_unit' : '%',
            'colorbar' : '%',
            'colormap' : 'jet_r',
            'scale' : [0.0, 100.0]
        },
        'T2' : {
            'name' : 'Temperature at 2m',
            'native_unit' : 'K',
            'colorbar' : 'F',
            'colormap' : 'jet',
            'scale' : [270.0, 320.0]
        },
        'PRECIP' : {
            'name' : 'Precipitation',
            'native_unit' : 'mm/h',
            'colorbar' : 'mm/h',
            'colormap' : 'jet_r',
            'scale' : [0.0, 2.0]
        },
        'HGT' : {
            'name' : 'Terrain height',
            'native_unit' : 'm',
            'colorbar' : 'm',
            'colormap' : 'jet_r',
            'scale' : [-86.0, 4500.0]
        },
    }

    show = ['HGT', 'T2', 'RH', 'PRECIP', 'EQUILd FM', 'EQUILw FM']

    esmf_cycle = utc_to_esmf(cycle) + timedelta(hours=fcst_hour)
    mf = { "1" : {esmf_cycle : {}}}
    ensure_dir(osp.join(postproc_path, manifest_name))
    
    if not is_cycle_computed(cycle, region_cfg, wksp_path, fcst_hour):
        logging.info(f'CYCLER copying postprocessing from cycle {prev_cycle}-f{prev_fcst_hour} to cycle {cycle}-f{fcst_hour}')
        prev_manifest_name = prev_cycle_id + '.json'
        prev_esmf_cycle = utc_to_esmf(prev_cycle) + timedelta(hours=fcst_hour)
        prev_manifest_path = osp.join(prev_postproc_path, prev_manifest_name)
        if not osp.exists(prev_manifest_path):
            logging.error('CYCLER previous post-processing could not be found')
            return None
        prev_mf = json.load(open(prev_manifest_path, 'r')) 
        for name in prev_mf['1'][prev_esmf_cycle].keys():
            prev_raster_name = prev_mf['1'][prev_esmf_cycle][name]['raster']
            prev_cb_name = prev_mf['1'][prev_esmf_cycle][name]['colorbar']
            raster_name = cycle_id + '-%s-raster.png' % name
            cb_name = cycle_id + '-%s-raster-cb.png' % name
            coords = prev_mf['1'][prev_esmf_cycle][name]['coords']
            alpha = prev_mf['1'][prev_esmf_cycle][name].get('alpha',None)
            force_copy(osp.join(prev_postproc_path, prev_raster_name),osp.join(postproc_path, raster_name))
            force_copy(osp.join(prev_postproc_path, prev_cb_name),osp.join(postproc_path, cb_name))
            if alpha:
                mf["1"][esmf_cycle][name] = { 'raster' : raster_name, 'coords' : coords, 'colorbar' : cb_name, 'alpha' : alpha }
            else:
                mf["1"][esmf_cycle][name] = { 'raster' : raster_name, 'coords' : coords, 'colorbar' : cb_name }
    else:
        if bounds is None:
            bounds = (region_cfg.bbox[1],region_cfg.bbox[3],region_cfg.bbox[0],region_cfg.bbox[2])
        # read in the longitudes and latitudes
        geo_path = osp.join(wksp_path, '%s-geo.nc' % region_cfg.code)
        logging.info('CYCLER reading longitudes and latitudes from NetCDF file %s' % geo_path )
        d = netCDF4.Dataset(geo_path)
        lats = d.variables['XLAT'][:,:]
        lons = d.variables['XLONG'][:,:]
        d.close()
        # read and process model variables
        with netCDF4.Dataset(model_path) as d:
            for name in show:
                raster_png, coords, cb_png, levels = scalar_field_to_raster(d.variables[name][:,:], lats, lons, var_wisdom[name])
                write_postprocess(mf, postproc_path, cycle_id, esmf_cycle, name, raster_png, coords, cb_png, levels, .5)
            for i,name in [(0, '1-hr DFM'), (1, '10-hr DFM'), (2, '100-hr DFM'), (3, '1000-hr DFM')]:
                fm_wisdom = var_wisdom['dfm']
                fm_wisdom['name'] = 'Estimated %s' % name
                raster_png, coords, cb_png, levels = scalar_field_to_raster(d.variables['FMC_GC'][:,:,i], lats, lons, fm_wisdom)
                write_postprocess(mf, postproc_path, cycle_id, esmf_cycle, name, raster_png, coords, cb_png, levels, .5)
        if osp.exists('src/ingest/MesoDB'):
            from ingest.MesoDB.mesoDB import mesoDB
            db = mesoDB('ingest/MesoDB')
            db.update['startTime'] = cycle - timedelta(hours=1)
            db.update['endTime'] = cycle + timedelta(hours=1)
            db.params['startTime'] = cycle - timedelta(hours=1)
            db.params['endTime'] = cycle + timedelta(hours=1)
            db.params['longitude1'], db.params['longitude2'], db.params['latitude1'], db.params['latitude2'] = bounds
            if is_cycle_computed(cycle, region_cfg, wksp_path):
                db.params['updateDB'] = False
            df = db.get_DB()
            st = db.sites()
            data = df.groupby('STID').mean().join(st[['LONGITUDE','LATITUDE']])
            meso_wisdom = var_wisdom['dfm']
            meso_wisdom['name'] = 'MesoWest 10-hr DFM'
            meso_wisdom['bbox'] = bounds
            meso_wisdom['text'] = False
            raster_png, coords, cb_png, levels = scatter_to_raster(np.array(data['fm10'])/100., 
                                                   np.array(data['LATITUDE']).astype(float), 
                                                   np.array(data['LONGITUDE']).astype(float), meso_wisdom) 
            name = 'MESO 10-hr DFM'
            write_postprocess(mf, postproc_path, cycle_id, esmf_cycle, name, raster_png, coords, cb_png, levels, 1.)

    logging.info('writing manifest file %s' % osp.join(postproc_path, manifest_name) )
    json.dump(mf, open(osp.join(postproc_path, manifest_name), 'w'), indent=1, separators=(',',':'))
    logging.info(json.dumps(mf))
    if osp.exists(osp.join(prev_postproc_path, complete_manifest_name)):
        complete_mf = json.load(open(osp.join(prev_postproc_path, complete_manifest_name), 'r'))
        complete_mf['1'].update(mf['1'])
        json.dump(complete_mf, open(osp.join(postproc_path, complete_manifest_name), 'w'), indent=1, separators=(',',':'))
    else:
        json.dump(mf, open(osp.join(postproc_path, complete_manifest_name), 'w'), indent=1, separators=(',',':'))

    return postproc_path

def compute_fmda_id(cycle, region_code):
    """
    Construct a fmda id unique for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :return: a unique fmda id
    """
    time_stamp = cycle.strftime('%Y%m%d')
    fmda_id = f'fmda-{region_code}-{time_stamp}-{cycle.hour:02d}'
    return fmda_id

def compute_cycle_path(cycle, region_code, wksp_path):
    """
    Construct a relative path to the cycle path for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :return: a relative path (w.r.t. workspace and region) of the cycle path
    """
    fmda_id = compute_fmda_id(cycle, region_code)
    year_month_folder = cycle.strftime('%Y%m')
    return osp.join(wksp_path, region_code, year_month_folder, fmda_id)

def compute_model_path(cycle, region_code, wksp_path, fcst_hour=0, ext='nc'):
    """
    Construct a relative path to the fuel moisture model file
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :param fcst_hour: forecast hour
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    fmda_id = compute_fmda_id(cycle, region_code)
    cycle_path = compute_cycle_path(cycle, region_code, wksp_path)
    if fcst_hour != 0:
        filename = f'{fmda_id}-f{fcst_hour:02d}.{ext}'
    else:
        filename = f'{fmda_id}.{ext}' 
    return osp.join(cycle_path, filename)

def compute_postproc_path(cycle, region_code, wksp_path, fcst_hour=0, ext='nc'):
    """
    Construct a relative path to the post-processing folder
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :param fcst_hour: forecast hour
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    year_month_folder = cycle.strftime('%Y%m')
    cycle_id = compute_fmda_id(cycle, region_code)
    if fcst_hour != 0:
        fmda_id = f'{cycle_id}-f{fcst_hour:02d}'
    else:
        fmda_id = f'{cycle_id}' 
    return osp.join(wksp_path, year_month_folder, fmda_id)

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
        tmp = np.where(np.amax(glat[:, j1:j2],axis=1) < minlat)[0]
        if len(tmp):
            tmp = tmp[-1]
        else:
            tmp = i1
        if i1 != tmp:
            i1 = tmp
            done = False
        tmp = np.where(np.amin(glat[:, j1:j2],axis=1) > maxlat)[0]
        if len(tmp):
            tmp = tmp[0]
        else:
            tmp = i2
        if i2 != tmp:
            i2 = tmp
            done = False
        tmp = np.where(np.amax(glon[i1:i2,:],axis=0) < minlon)[0]
        if len(tmp):
            tmp = tmp[-1]
        else:
            tmp = j1
        if j1 != tmp:
            j1 = tmp
            done = False
        tmp = np.where(np.amin(glon[i1:i2,:],axis=0) > maxlon)[0]
        if len(tmp):
            tmp = tmp[0]
        else:
            tmp = j2
        if j2 != tmp:
            j2 = tmp
            done = False
    return i1,i2,j1,j2


def compute_hrrr_bounds(bbox):
    """
    Compute bounds from HRRR data even when HRRR data is not available from terrain static data
    
    :param bbox: the bounding box of the data
    :return: a tuple containing bound coordinates (min_lon,max_lon,min_lat,max_lat)
    """
    ds = netCDF4.Dataset('static/hrrr.terrainh.nc')
    lats,lons = ds['XLAT_M'][0], ds['XLONG_M'][0]
    i1, i2, j1, j2 = find_region_indices(lats, lons, bbox[0], bbox[2], bbox[1], bbox[3])
    lats,lons = lats[i1:i2,j1:j2], lons[i1:i2,j1:j2]
    return (lons.min(), lons.max(), lats.min(), lats.max())


def load_hrrr_data(grib_file, bbox):
    """
    Load relevant GRIB fields and return them
    
    :param grib_file: path to HRRR grib file
    :param bbox: the bounding box of the data
    :return: a dictionary with variables
    """
    gf = GribFile(grib_file)
    lats, lons = gf[1].latlons()
    # bbox format: minlat, minlon, maxlat, maxlon
    i1, i2, j1, j2 = find_region_indices(lats, lons, bbox[0], bbox[2], bbox[1], bbox[3])
   
    lats = lats[i1:i2,j1:j2] 
    lons = lons[i1:i2,j1:j2]
    hgt = np.ma.array(netCDF4.Dataset('static/hrrr.terrainh.nc')['HGT_M'][0])[i1:i2,j1:j2]
    data = {'lats': lats, 'lons': lons, 'hgt': hgt}
    
    for v,idx in hrrr_vars.items():
        var = np.ma.array(gf[idx].values())[i1:i2,j1:j2]
        logging.info('{} min {} max {}'.format(v, np.min(var), np.max(var)))
        logging.info('{} min {} max {}'.format(v, np.min(var), np.max(var)))
        data.update({v: var})

    return data


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


def fmda_advance_region(cycle, cfg, grib_files, wksp_path, lookback_length, fcst_hour, meso_token):
    """
    Advance the fuel moisture estimates in the region specified by the configuration.
    The function assumes that the fuel moisture model has not been advanced to this
    cycle yet and will overwrite any previous computations.
    
    Control flow:
    
    1) read in HRRR variables
    2) check if there is a stored FM model for previous cycle
    2a) yes -> load it, advance one time-step, perform DA
    2b) no -> compute equilibrium, use background covariance to do DA
    3) store model
    
    :param cycle: the datetime indicating the processed cycle in UTC
    :param cfg: the configuration dictionary specifying the region
    :param grib_files: path to HRRR grib files to retrieve variables for this cycle (or previous)
    :param wksp_path: the workspace path for the cycler
    :param lookback_length: nubmer of cycles to search before we find a computed cycle
    :param forecast_length: number of cycles to forecast
    :param meso_token: the mesowest API access token or a list of them
    :return: the model advanced and assimilated at the current cycle
    """
    min_num_obs = 10
    max_fm10_value = 50.
    logging.info("hrrr_cycler.fmda_advance_region: %s" % str(cycle))
    model = None
    if fcst_hour == 0:
        prev_cycle = cycle - timedelta(hours=1)
        prev_model_path = compute_model_path(prev_cycle, cfg.code, wksp_path, fcst_hour)
    else:
        prev_cycle = cycle
        prev_model_path = compute_model_path(prev_cycle, cfg.code, wksp_path, fcst_hour-1)
    if not osp.exists(prev_model_path):
        logging.info('CYCLER cannot find model from previous cycle %s' % str(prev_cycle))
        if lookback_length > 0:
            model = fmda_advance_region(cycle - timedelta(hours=1), cfg, grib_files, wksp_path, lookback_length - 1, fcst_hour, meso_token)
    else:
        logging.info('CYCLER found previous model for cycle %s.' % str(prev_cycle))
        model = FuelMoistureModel.from_netcdf(prev_model_path)
    
    grib_file = grib_files[lookback_length]
    # retrieve the variables and make sure they are available (we should not be here if they are not)
    if not osp.exists(grib_file):
        logging.warning('CYCLER could not find useable cycle.')
        logging.error(e)
        logging.warning('CYCLER copying previous post-processing.')
        try:
            bounds = compute_hrrr_bounds(cfg.bbox)
            pp_path = postprocess_cycle(cycle, cfg, wksp_path, fcst_hour, bounds)   
            if pp_path != None:
                if 'shuttle_remote_host' in sys_cfg:
                    sim_code = 'fmda-' + cfg.code
                    send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', cfg.region_id + ' FM')
        except Exception as e:
            logging.warning('CYCLER exception {}'.format(e))
            logging.error('CYCLER skipping region {} for cycle {}'.format(cfg.region_id,str(cycle)))
        sys.exit(1) 
    
    logging.info('CYCLER loading HRRR data for cycle %s.' % str(cycle))
    data = load_hrrr_data(grib_file, cfg.bbox)
    lats = data['lats']
    lons = data['lons']
    hgt = data['hgt']
    Ed, Ew = compute_equilibria(data['t2'], data['rh'])
    
    rain = data['rain'][:,:] + 0
    # remove rain that is too small to make any difference 
    rain[rain < 0.01] = 0
    # remove bogus rain that is too large 
    rain[rain > 1e10] = 0
    # remove masked rain values
    rain[rain.mask] = 0

    dom_shape = data['t2'].shape
    # store the lons/lats for this domain
    geo_path = osp.join(wksp_path, '%s-geo.nc' % cfg.code)
    if not osp.isfile(geo_path):
        logging.info('CYCLER initializing new file %s.' % (geo_path))
        d = netCDF4.Dataset(geo_path, 'w', format='NETCDF4')
        d.createDimension('south_north', dom_shape[0])
        d.createDimension('west_east', dom_shape[1])
        xlat = d.createVariable('XLAT', 'f4', ('south_north', 'west_east'))
        xlat[:,:] = lats
        xlong = d.createVariable('XLONG', 'f4', ('south_north', 'west_east'))
        xlong[:,:] = lons
        d.close()
    else:
        logging.info('CYCLER file already exists:  %s.' % (geo_path))

    # check if we must start from equilibrium
    if model is None:
        logging.info('CYCLER initializing from equilibrium for cycle %s.' % (str(cycle)))
        # setup model parameters    
        Nk = 3
        Tk = np.array([1.0, 10.0, 100.0, 1000.0]) * 3600
        m0 = np.expand_dims(0.5 * (Ed + Ew), axis=2)
        # background covariance
        P0 = np.diag([0.01, 0.01, 0.01, 0.01, 0.001, 0.001])
        model = FuelMoistureModel(m0[:,:,[0, 0, 0, 0]], Tk, P0)
    else:
        logging.info('CYCLER advancing model one hour to cycle %s.' % (str(cycle)))
        # always 1 hr step in HRRR
        dt = 3600 
        # the process noise matrix
        Q = np.diag([1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-5])
        model.advance_model(Ed, Ew, rain, dt, Q)

    logging.info('CYCLER retrieving fm-10 observations for cycle %s.' % (str(cycle)))
    # no assimilation
    if fcst_hour > 0:
        logging.info('CYCLER forecasting mode, skipping data assimilation')
    else:
        # perform assimilation with mesowest observations
        tm_start = cycle - timedelta(minutes=30)
        tm_end = cycle + timedelta(minutes=30)
        fm10 = retrieve_mesowest_observations(meso_token, tm_start, tm_end, lats, lons, hgt)

        # filter fm10 values for statistics
        valid_times = [z for z in fm10.keys() if abs((z - cycle).total_seconds()) < 1800]
        fm10_filter = {}
        obs_valid_now = []
        for z in valid_times:
            vobs = [f for f in fm10[z] if f.obs_val > 0. and f.obs_val < max_fm10_value]
            obs_valid_now.extend(vobs)
            fm10_filter.update({z: vobs})
        fm10v = [obs.get_value() for obs in obs_valid_now]
        fm10 = fm10_filter
        
        logging.info(
            'CYCLER retrieved %d valid observations at %d unique times, min/mean/max [%g/%g/%g].' % (
                len(fm10v), len(valid_times), np.amin(fm10v), np.mean(fm10v), np.amax(fm10v))
        )
        
        if len(obs_valid_now) > min_num_obs:
            # run the data assimilation step
            covs = [np.ones(dom_shape), hgt, lats, lons]
            covs_names = ['const', 'hgt', 'lat', 'lon']
            if np.any(rain > 0.01):
                covs.append(rain)
                covs_names.append('rain')
            if np.any(data['snow'] > 0.01):
                covs.append(data['snow'])
                covs_names.append('snow')
            other_covs = [
                't2', 'rh', 'psfc', 'soil_t_0', 'soil_t_1', 'soil_t_2', 
                'soil_t_3', 'soil_t_4', 'soil_t_5', 'soil_t_6', 'soil_t_7', 
                'soil_t_8', 'soil_moist_0', 'soil_moist_1', 'soil_moist_2', 
                'soil_moist_3', 'soil_moist_4', 'soil_moist_5', 'soil_moist_6', 
                'soil_moist_7', 'soil_moist_8', 'snowh', 'spfh', 'u10', 'v10', 
                'ws', 'ust', 'znt', 'swdown'
            ]
            for cov in other_covs:
                covs_names.append(cov)
                covs.append(data[cov])
                
            execute_da_step(model, cycle, covs, covs_names, fm10, use_lstsq=True)
    
    # make geogrid files for WPS; datasets and lines to add to GEOGRID.TBL
    geo_path = compute_model_path(cycle, cfg.code, wksp_path, fcst_hour, ext="geo")
    index = {
        'projection': 'lambert',
        'dx' : 3000.0,
        'dy' : -3000.0,
        'truelat1' : 38.5,
        'truelat2' : 38.5,
        'stdlon' : 262.5,
        'radius' : 6370000.0
    }
    model.to_geogrid(geo_path, index, lats, lons)

    # make wps format files for WPS
    fmda_path = osp.join(wksp_path, cfg.code, '{:04d}{:02d}'.format(cycle.year,cycle.month))
    time_tag = '{:04d}-{:02d}-{:02d}_{:02d}'.format(cycle.year, cycle.month, cycle.day, cycle.hour)
    model.to_wps_format(fmda_path, index, lats, lons, time_tag)
    
    # store the new model  
    model_path = compute_model_path(cycle, cfg.code, wksp_path)
    logging.info('CYCLER writing model variables to:  %s.' % model_path)
    model.to_netcdf(
        ensure_dir(model_path), {
            'EQUILd FM': Ed, 'EQUILw FM': Ew, 'T2': data['t2'], 
            'RH': data['rh'], 'PRECIP': rain, 'HGT': hgt
        }
    )

    # create visualization and send results
    bounds = (lons.min(), lons.max(), lats.min(), lats.max())
    pp_path = postprocess_cycle(cycle, cfg, wksp_path, fcst_hour, bounds)   
    if pp_path != None:
        if 'shuttle_remote_host' in sys_cfg:
            sim_code = 'fmda-hrrr-' + cfg.code
            send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', cfg.region_id + ' FM')
    
    return model
    
    
def is_cycle_computed(cycle, cfg, wksp_path, fcst_hour=0):
    """
    Check if the fuel model file exists (has been computed) for the
    cycle <cycle> and region configuration <cfg>.
    
    :param cycle: the cycle datetime in UTC
    :param cfg: the region configuration wrapped in a Dict for convenience
    :param wksp_path: the workspace path for the cycler
    :return: True if the model file has been found, False otherwise
    """
    path = compute_model_path(cycle, cfg.code, wksp_path, fcst_hour=fcst_hour)
    return osp.isfile(path)
    
    
if __name__ == '__main__':
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if len(sys.argv) > 1:
        mode = sys.argv[1]
        if mode in ['a', 'A', 'f', 'F']:
            mode = mode.lower()
        else:
            mode = None
    if len(sys.argv) == 3:
        code = sys.argv[2]
        for k,region in cfg.regions.items():
            if region['code'] == code:
                cfg.regions = {
                    k: region
                }
                break
    elif len(sys.argv) == 6:
        mode = sys.argv[1]
        code = 'FIRE'
        cfg.regions = {
             "Fire domain" : {
                  "code" : code,
                  "bbox" : sys.argv[2:6]
             }
        }
        try:
            os.remove(osp.join(cfg.workspace_path,code+'-geo.nc'))
        except Exception as e:
            logging.warning(e)
        try:
            delete(osp.join(cfg.workspace_path,code))
        except Exception as e:
            logging.warning(e)

    if mode is None or len(cfg.regions) < 1:
        print('Usage: to use domains configured in etc/fmda_cycler.json')
        print('{} mode code'.format(sys.argv[0]))
        print('The supported modes are: analysis (a) and forecast (f)')
        print('To use a custom domain named FIRE by giving a bounding box:')
        print('./hrrr_cycler.sh mode lat1 lon1 lat2 lon2')
        print('Example: ./hrrr_cycler.sh a 42 -124.6 49 -116.4')
        exit(1) 

    # get parameters from configuration
    lookback_length = cfg.get('lookback_length', 24)
    forecast_length = cfg.get('forecast_length', 48)
    period_hours = cfg.get('period_hours', 6)
    # get more readable mode
    mode_name = 'analysis' if mode == 'a' else 'forecast'
    # current time
    now = datetime.now(pytz.UTC)
    cycle = (now - timedelta(minutes=59)).replace(minute=0,second=0,microsecond=0)
    # print statements
    logging.info(
        'CYCLER activated at {now}, will attempt cycle at {cycle} with mode {mode_name}'
    )
    logging.info('regions: {}'.format(json.dumps(cfg.regions)))

    # getting necessary data
    try:
        from_utc = cycle - timedelta(hours=lookback_length)
        to_utc = cycle
        hrrra = HRRRA(sys_cfg)
        anl_gribs = hrrra.retrieve_gribs(from_utc, to_utc) 
        grib_files_anl = anl_gribs['grib_files']
        if mode == 'f':
            hrrr = HRRR(sys_cfg)
            shift_hours = cycle.hour % period_hours
            cycle_start = cycle - timedelta(hours=shift_hours)
            from_utc = cycle + timedelta(hours=1)
            to_utc = cycle_start + timedelta(hours=forecast_length)
            fct_gribs = hrrr.retrieve_gribs(from_utc, to_utc, cycle_start=cycle_start)
            grib_files_fct = fct_gribs['grib_files']  
    except:
        logging.warning('CYCLER could not find useable cycle.')
        logging.warning('CYCLER copying previous post-processing.')
        for region_id,region_cfg in cfg.regions.items():
            wrapped_cfg = Dict(region_cfg)
            wrapped_cfg.update({'region_id': region_id})
            try:
                bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, bounds)
                if pp_path != None:
                    if 'shuttle_remote_host' in sys_cfg:
                        sim_code = 'fmda-' + wrapped_cfg.code
                        send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', region_id + ' FM')
            except Exception as e:
                logging.warning('CYCLER exception {}'.format(e))
                logging.error(f'CYCLER skipping region {region_id} for cycle {cycle}')
        sys.exit(1)
    logging.info(f'have necessary HRRR data for cycle {cycle} at mode {mode_name}.')
    
    # check for each region, if we are up to date w.r.t. HRRR data available
    for region_id,region_cfg in cfg.regions.items():
        logging.info(f'CYCLER processing region {region_id} for {cycle}')
        wrapped_cfg = Dict(region_cfg)
        wrapped_cfg.update({'region_id': region_id})
        # real-time part
        if not is_cycle_computed(cycle, wrapped_cfg, cfg.workspace_path):
            logging.info(f'CYCLER real-time processing for region {region_id} at cycle {cycle}')
            try:
                fmda_advance_region(cycle, wrapped_cfg, grib_files_anl, cfg.workspace_path, lookback_length, 0, meso_token)
            except Exception as e:
                logging.warning(f'CYCLER failed real-time processing for region {region_id} at cycle {cycle}')
                logging.warning('CYCLER exception {}'.format(e))
                logging.warning('CYCLER copying previous post-processing or re-trying.')
                try:
                    bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                    pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, bounds)   
                    if pp_path != None:
                        if 'shuttle_remote_host' in sys_cfg:
                            sim_code = 'fmda-hrrr-' + wrapped_cfg.code
                            send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', region_id + ' FM')
                except Exception as e:
                    logging.error(f'CYCLER skipping region {region_id} for cycle {cycle} and mode {mode_name}')
        else:
            logging.info(f'CYCLER already completed real-time processing for region {region_id} at cycle {cycle}, skipping ...')
        # forecasting part
        if mode == 'f':
            logging.info(f'CYCLER forecasting region {region_id} at cycle {cycle}') 
            fcst_hour = 1
            while cycle + timedelta(hours=fcst_hour) < to_utc:
                if not is_cycle_computed(cycle, wrapped_cfg, cfg.workspace_path, fcst_hour=fcst_hour):
                    logging.info(f'CYCLER forecasting region {region_id} at time {cycle}') 
                    try:
                        fmda_advance_region(
                            cycle, wrapped_cfg, 
                            grib_files_fct[fcst_hour-1:fcst_hour],
                            cfg.workspace_path, 0, fcst_hour, meso_token
                        )
                    except Exception as e:
                        logging.warning(f'CYCLER failed forecasting for region {region_id} at time {cycle}')
                        logging.warning('CYCLER exception {}'.format(e))
                        logging.warning('CYCLER copying previous post-processing or re-trying.')
                        try:
                            bounds = compute_hrrr_bounds(wrapped_cfg.bbox)
                            pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, bounds)
                            if pp_path != None:
                                if 'shuttle_remote_host' in sys_cfg:
                                    sim_code = 'fmda-hrrr-' + wrapped_cfg.code
                                    send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', region_id + ' FM')
                        except Exception as e:
                            logging.error(f'CYCLER skipping region {region_id} for cycle {cycle} and mode {mode_name}')
                else:
                    logging.info(f'CYCLER already completed real-time processing for region {region_id} at cycle {cycle}, skipping ...')

                fcst_hour += 1

    # done
    logging.info(f'CYCLER cycle {cycle} complete with mode {mode_name}.')
