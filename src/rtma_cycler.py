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
from ingest.rtma_source import RTMA
from utils import Dict, ensure_dir, utc_to_esmf, delete, force_copy, move
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
cfg = Dict(json.load(open('etc/rtma_cycler.json')))
meso_token = json.load(open('etc/tokens.json'))['mesowest']


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


def postprocess_cycle(cycle, region_cfg, wksp_path, bounds=None):
    """
    Build rasters from the computed fuel moisture.

    :param cycle: the UTC cycle time
    :param region_cfg: the region configuration
    :param wksp_path: the workspace path
    :param bounds: bounding box of the post-processing
    :return: the postprocessing path
    """
    prev_cycle = cycle-timedelta(hours=1)
    post_cycle = cycle+timedelta(hours=1)
    model_path = compute_model_path(cycle, region_cfg.code, wksp_path)
    year_month = '%04d%02d' % (cycle.year, cycle.month)
    prev_year_month = '%04d%02d' % (prev_cycle.year, prev_cycle.month)
    cycle_dir = 'fmda-%s-%04d%02d%02d-%02d' %  (region_cfg.code, cycle.year, cycle.month, cycle.day, cycle.hour)
    prev_cycle_dir = 'fmda-%s-%04d%02d%02d-%02d' %  (region_cfg.code, prev_cycle.year, prev_cycle.month, prev_cycle.day, prev_cycle.hour)
    postproc_path = osp.join(wksp_path, year_month, cycle_dir)
    prev_postproc_path = osp.join(wksp_path, prev_year_month, prev_cycle_dir)
    manifest_name = cycle_dir + '.json'
    complete_manifest_name = 'fmda-%s.json' % region_cfg.code

    if not is_cycle_computed(cycle, region_cfg, wksp_path) and not osp.exists(prev_postproc_path):
        logging.warning('CYCLER postprocessing failed for time {}'.format(str(cycle)))
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
        'TD' : {
            'name' : 'Dew point temperature at 2m',
            'native_unit' : 'K',
            'colorbar' : 'F',
            'colormap' : 'jet',
            'scale' : [270.0, 320.0]
        },
        'T2' : {
            'name' : 'Temperature at 2m',
            'native_unit' : 'K',
            'colorbar' : 'F',
            'colormap' : 'jet',
            'scale' : [270.0, 320.0]
        },
        'PRECIPA' : {
            'name' : 'RTMA precipa',
            'native_unit' : 'kg/m^2/h',
            'colorbar' : 'kg/m^2/h',
            'colormap' : 'jet_r',
            'scale' : [0.0, 2.0]
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

    show = ['TD','PRECIPA','T2','HGT','PRECIP','RH','EQUILd FM','EQUILw FM']
    show = ['T2','HGT','PRECIP','RH']

    esmf_cycle = utc_to_esmf(cycle) 
    mf = { "1" : {esmf_cycle : {}}}
    ensure_dir(osp.join(postproc_path, manifest_name))
    
    if not is_cycle_computed(cycle, region_cfg, wksp_path):
        logging.info('CYCLER copying postprocessing from cycle {} to cycle {}'.format(str(prev_cycle),str(cycle)))
        prev_manifest_name = prev_cycle_dir + '.json'
        prev_esmf_cycle = utc_to_esmf(prev_cycle)
        prev_mf = json.load(open(osp.join(prev_postproc_path, prev_manifest_name), 'r')) 
        for name in prev_mf['1'][prev_esmf_cycle].keys():
            prev_raster_name = prev_mf['1'][prev_esmf_cycle][name]['raster']
            prev_cb_name = prev_mf['1'][prev_esmf_cycle][name]['colorbar']
            raster_name = cycle_dir + '-%s-raster.png' % name
            cb_name = cycle_dir + '-%s-raster-cb.png' % name
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
                write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels, .5)
            for i,name in [(0, '1-hr DFM'), (1, '10-hr DFM'), (2, '100-hr DFM')]:
                fm_wisdom = var_wisdom['dfm']
                fm_wisdom['name'] = 'Estimated %s' % name
                raster_png, coords, cb_png, levels = scalar_field_to_raster(d.variables['FMC_GC'][:,:,i], lats, lons, fm_wisdom)
                write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels, .5)
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
            write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels, 1.)
        # NFMDB observations
        if osp.exists('src/ingest/FMDB'):
            from ingest.FMDB.FMDB import FMDB
            from ingest.FMDB.utils import filter_outliers
            period_length = 7 # period in days
            period_num = np.ceil(cycle.day/period_length)
            db = FMDB('ingest/NFMDB')
            db.params['startYear'] = 2019
            data = db.get_data()
            data = filter_outliers(data) 
            data['fuel_type'] = data['fuel_type'].fillna('None').str.upper()
            data['fuel_variation'] = data['fuel_variation'].fillna('None').str.upper()
            sts = db.sites()
            data = data.join(sts[['lng','lat']],'site_number')
            # mask space
            lats = data['lat']
            lons = data['lng']
            data = data[np.logical_and(lats <= bounds[3],
                            np.logical_and(lats >= bounds[2],
                                np.logical_and(lons <= bounds[1],
                                               lons >= bounds[0])))]
            dates = data['date'].dt.tz_localize(pytz.UTC)
            # calculate top 5 LFM to always plot the same
            top = 5
            hist_data = data[dates.dt.year <= 2020]
            hist_dfm_mask = np.array(['-HOUR' in ft for ft in np.array(hist_data['fuel_type'])]).astype(bool)
            hist_df_lfm = hist_data[~hist_dfm_mask].reset_index(drop=True)
            fts = np.array(hist_df_lfm[['fuel_type','percent']].groupby('fuel_type').count().sort_values(by='percent',ascending=False).index[:top])
            # mask time 
            start = cycle.replace(day=int(period_length*(period_num-1)+1),hour=0,minute=0,second=0,microsecond=0)
            end = cycle
            data = data[np.logical_and(dates >= start, dates <= end)]
            cycle_dir = 'fmda-%s-%04d%02d%02d-%02d' %  (region_cfg.code, start.year, start.month, start.day, start.hour)
            # mask dead and live fuel moisture
            dfm_mask = np.array(['-HOUR' in ft for ft in np.array(data['fuel_type'])]).astype(bool)
            df_dfm = data[dfm_mask].reset_index(drop=True)
            df_lfm = data[~dfm_mask].reset_index(drop=True)
            # plot NFMDB dead fuel moisture
            for i,name in [('1-HOUR','NFMDB 1-hr DFM'),('10-HOUR','NFMDB 10-hr DFM'),('100-HOUR','NFMDB 100-hr DFM'),('1000-HOUR','NFMDB 1000-hr DFM')]:
                fmdb_wisdom = var_wisdom['dfm']
                fmdb_wisdom['name'] = name
                fmdb_wisdom['bbox'] = bounds
                fmdb_wisdom['text'] = True
                fmdb_wisdom['size'] = 40
                fmdb_wisdom['linewidth'] = 1.
                data = df_dfm[df_dfm['fuel_type'] == i]
                raster_png, coords, cb_png, levels = scatter_to_raster(np.array(data['percent'])/100., 
                                                                   np.array(data['lat']), 
                                                                   np.array(data['lng']), fmdb_wisdom) 
                write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels, 1.)
            # plot NFMDB live fuel moisture
            df_lfm = df_lfm.sort_values('date').groupby(['site_number','fuel_type']).last().reset_index()
            for ft in fts:
                name = 'NFMDB {} LFM'.format(ft)
                fmdb_wisdom = var_wisdom['lfm']
                fmdb_wisdom['name'] = name
                fmdb_wisdom['bbox'] = bounds
                fmdb_wisdom['text'] = True
                fmdb_wisdom['size'] = 40
                fmdb_wisdom['linewidth'] = 1.
                data = df_lfm[df_lfm['fuel_type'] == ft]
                raster_png, coords, cb_png, levels = scatter_to_raster(np.array(data['percent'])/100., 
                                                                   np.array(data['lat']), 
                                                                   np.array(data['lng']), fmdb_wisdom)
                write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels, 1.)
            name = 'NFMDB OTHERS LFM'
            fmdb_wisdom = var_wisdom['lfm']
            fmdb_wisdom['name'] = name
            fmdb_wisdom['bbox'] = bounds
            fmdb_wisdom['text'] = True
            fmdb_wisdom['size'] = 40
            fmdb_wisdom['linewidth'] = 1.
            data = df_lfm[~df_lfm['fuel_type'].isin(fts)]
            data = data.groupby('site_number').mean()
            raster_png, coords, cb_png, levels = scatter_to_raster(np.array(data['percent'])/100.,
                                                               np.array(data['lat']),
                                                               np.array(data['lng']), fmdb_wisdom)
            write_postprocess(mf, postproc_path, cycle_dir, esmf_cycle, name, raster_png, coords, cb_png, levels, 1.)

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


def compute_model_path(cycle, region_code, wksp_path, ext='nc'):
    """
    Construct a relative path to the fuel moisture model file
    for the region code and cycle.
    
    :param cycle: the UTC cycle time
    :param region_code: the code of the region
    :param wksp_path: the workspace path
    :return: a relative path (w.r.t. workspace and region) of the fuel model file
    """
    year_month = '%04d%02d' % (cycle.year, cycle.month)
    filename = 'fmda-%s-%04d%02d%02d-%02d.%s' %  (region_code, cycle.year, cycle.month, cycle.day, cycle.hour, ext)
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


def compute_rtma_bounds(bbox):
    """
    Compute bounds from RTMA data even when RTMA data is not available from terrain static data
    
    :param bbox: the bounding box of the data
    :return: a tuple containing bound coordinates (min_lon,max_lon,min_lat,max_lat)
    """
    lats,lons = GribFile('static/ds.terrainh.bin')[1].latlons()
    i1, i2, j1, j2 = find_region_indices(lats, lons, bbox[0], bbox[2], bbox[1], bbox[3])
    lats,lons = lats[i1:i2,j1:j2], lons[i1:i2,j1:j2]
    return (lons.min(), lons.max(), lats.min(), lats.max())


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
    
    t2 = np.ma.array(gf.values())[i1:i2,j1:j2] # temperature at 2m in K
    td = np.ma.array(GribFile(rtma_data['td'])[1].values())[i1:i2,j1:j2] # dew point in K
    precipa = np.ma.array(GribFile(rtma_data['precipa'])[1].values())[i1:i2,j1:j2] # precipitation
    hgt = np.ma.array(GribFile('static/ds.terrainh.bin')[1].values())[i1:i2,j1:j2]
    logging.info('t2 min %s max %s' % (np.min(t2),np.max(t2)))
    logging.info('td min %s max %s' % (np.min(td),np.max(td)))
    logging.info('precipa min %s max %s' % (np.min(precipa),np.max(precipa)))
    logging.info('hgt min %s max %s' % (np.min(hgt),np.max(hgt)))
    
    # compute relative humidity
    rh = 100*np.exp(17.625*243.04*(td - t2) / (243.04 + t2 - 273.15) / (243.0 + td - 273.15))
    
    return td, t2, rh, precipa, hgt, lats[i1:i2,j1:j2], lons[i1:i2,j1:j2]


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
    :param meso_token: the mesowest API access token or a list of them
    :return: the model advanced and assimilated at the current cycle
    """
    logging.info("rtma_cycler.fmda_advance_region: %s" % str(cycle))
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
    try:
        dont_have_vars, have_vars = rtma.retrieve_rtma(cycle)
        assert not dont_have_vars
    except AssertionError or ValueError as e:
        logging.warning('CYCLER could not find useable cycle.')
        logging.error(e)
        logging.warning('CYCLER copying previous post-processing.')
        try:
            bounds = compute_rtma_bounds(cfg.bbox)
            pp_path = postprocess_cycle(cycle, cfg, wksp_path, bounds)   
            if pp_path != None:
                if 'shuttle_remote_host' in sys_cfg:
                    sim_code = 'fmda-' + cfg.code
                    send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', cfg.region_id + ' FM')
        except Exception as e:
            logging.warning('CYCLER exception {}'.format(e))
            logging.error('CYCLER skipping region {} for cycle {}'.format(cfg.region_id,str(cycle)))
        sys.exit(1) 
    
    logging.info('CYCLER loading RTMA data for cycle %s.' % str(cycle))
    TD, T2, RH, precipa, hgt, lats, lons = load_rtma_data(have_vars, cfg.bbox)
    Ed, Ew = compute_equilibria(T2, RH)

    rain = precipa[:,:] + 0
    # remove rain that is too small to make any difference 
    rain[rain < 0.01] = 0
    # remove bogus rain that is too large 
    rain[rain > 1e10] = 0
    # remove masked rain values
    rain[rain.mask] = 0

    dom_shape = T2.shape

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
    fm10 = retrieve_mesowest_observations(meso_token, tm_start, tm_end, lats, lons, hgt)
    fm10v = []
    for fm10_obs in fm10.values():
        for obs in fm10_obs:
            fm10v.append(obs.get_value())
    
    logging.info('CYCLER retrieved %d valid observations, min/mean/max [%g/%g/%g].' %
                 (len(fm10),np.amin(fm10v),np.mean(fm10v),np.amax(fm10v)))
    
    # run the data assimilation step
    covs = [np.ones(dom_shape), hgt / 2000.0]
    covs_names = ['const','hgt/2000']
    if np.any(rain > 0.01):
        covs.append(rain)
        covs_names.append('rain')
    execute_da_step(model, cycle, covs, covs_names, fm10)
    
    # make geogrid files for WPS; datasets and lines to add to GEOGRID.TBL
    geo_path = compute_model_path(cycle, cfg.code, wksp_path,ext="geo")
    index = rtma.geogrid_index()
    print('index',index)
    model.to_geogrid(geo_path,index,lats,lons)

    # make wps format files for WPS
    fmda_path = osp.join(wksp_path,cfg.code,'{:04d}{:02d}'.format(cycle.year,cycle.month))
    time_tag = '{:04d}-{:02d}-{:02d}_{:02d}'.format(cycle.year, cycle.month, cycle.day, cycle.hour)
    model.to_wps_format(fmda_path,index,lats,lons,time_tag)
    
    # store the new model  
    model_path = compute_model_path(cycle, cfg.code, wksp_path)
    logging.info('CYCLER writing model variables to:  %s.' % model_path)
    model.to_netcdf(ensure_dir(model_path),
        {'EQUILd FM':Ed,'EQUILw FM':Ew,'TD':TD,'T2':T2,'RH':RH,'PRECIPA':precipa,'PRECIP':rain,'HGT':hgt})

    # create visualization and send results
    bounds = (lons.min(), lons.max(), lats.min(), lats.max())
    pp_path = postprocess_cycle(cycle, cfg, wksp_path, bounds)   
    if pp_path != None:
        if 'shuttle_remote_host' in sys_cfg:
            sim_code = 'fmda-' + cfg.code
            send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', cfg.region_id + ' FM')
    
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
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if len(sys.argv) == 2:
        pass
    elif len(sys.argv) == 5:
        code = 'FIRE'
        cfg.regions = {
             "Fire domain" : {
                  "code" : code,
                  "bbox" : sys.argv[1:5]
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
    else:
        print('Usage: to use domains configured in etc/rtma_cycler.json:')
        print('./rtma_cycler.sh anything')
        print('To use a custom domain named FIRE by giving a bounding box:')
        print('./rtma_cycler.sh lat1 lon1 lat2 lon2')
        print('Example: ./rtma_cycler.sh 42, -124.6, 49, -116.4')
        exit(1) 

    logging.info('regions: %s' % json.dumps(cfg.regions))
    #logging.info('regions: %s' % json.dumps(cfg.regions, indent=1, separators=(',',':')))


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
        logging.warning('CYCLER could not find useable cycle.')
        logging.warning('CYCLER copying previous post-processing.')
        for region_id,region_cfg in cfg.regions.items():
            wrapped_cfg = Dict(region_cfg)
            wrapped_cfg.update({'region_id': region_id})
            try:
                bounds = compute_rtma_bounds(wrapped_cfg.bbox)
                pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, bounds)
                if pp_path != None:
                    if 'shuttle_remote_host' in sys_cfg:
                        sim_code = 'fmda-' + wrapped_cfg.code
                        send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', region_id + ' FM')
            except Exception as e:
                logging.warning('CYCLER exception {}'.format(e))
                logging.error('CYCLER skipping region {} for cycle {}'.format(region_id,str(cycle)))
        sys.exit(1)
        
    logging.info('Have RTMA data for cycle %s.' % str(cycle))
      
    # check for each region, if we are up to date w.r.t. RTMA data available
    for region_id,region_cfg in cfg.regions.items():
        wrapped_cfg = Dict(region_cfg)
        wrapped_cfg.update({'region_id': region_id})
        #if 1:   # to run every time for debugging
        if not is_cycle_computed(cycle, wrapped_cfg, cfg.workspace_path):
            logging.info('CYCLER processing region %s for cycle %s' % (region_id, str(cycle)))
            try:
                fmda_advance_region(cycle, wrapped_cfg, rtma, cfg.workspace_path, lookback_length, meso_token)
            except Exception as e:
                logging.warning('CYCLER failed processing region {} for cycle {}'.format(region_id,str(cycle)))
                logging.warning('CYCLER exception {}'.format(e))
                logging.warning('CYCLER copying previous post-processing or re-trying.')
                try:
                    bounds = compute_rtma_bounds(wrapped_cfg.bbox)
                    pp_path = postprocess_cycle(cycle, wrapped_cfg, cfg.workspace_path, bounds)   
                    if pp_path != None:
                        if 'shuttle_remote_host' in sys_cfg:
                            sim_code = 'fmda-' + wrapped_cfg.code
                            send_product_to_server(sys_cfg, pp_path, sim_code, sim_code, sim_code + '.json', region_id + ' FM')
                except Exception as e:
                    logging.warning('CYCLER exception {}'.format(e))
                    logging.error('CYCLER skipping region {} for cycle {}'.format(region_id,str(cycle)))
        else:
            logging.info('CYCLER the cycle %s for region %s is already complete, skipping ...' % (str(cycle), str(region_id)))

    # done
    logging.info('CYCLER cycle %s complete.' % str(cycle))
