from __future__ import absolute_import
from __future__ import unicode_literals
try:
    from osgeo import gdal, osr
except:
    import gdal, osr
import pyproj

import os, logging, json, sys, glob, re, traceback
import os.path as osp
import numpy as np
import netCDF4 as nc4

from vis.postprocessor import Postprocessor
from utils import load_sys_cfg, Dict, make_clean_dir
from vis.var_wisdom import get_wisdom, is_windvec, is_fire_var, strip_end
from six.moves import range

def scalar2tiffs(output_path, d, wisdom, projection, geot, times, var, ndv=-9999.0):
    '''
    Creates new GeoTiffs for each time from 2D array
    '''
    # read wisdom variables
    scale = wisdom['scale']
    if 'transparent_values' in wisdom:
        rng = wisdom['transparent_values']

    array = wisdom['retrieve_as'](d,0)
    if is_fire_var(var):
        fm,fn = strip_end(d)
        array = array[:fm,:fn]
    ysize,xsize = array.shape
    zsize = len(times)
    # write each slice of the array along the zsize
    tiffiles = []
    for i in range(zsize):
        # get array
        array = wisdom['retrieve_as'](d,i)
        if is_fire_var(var):
            fm,fn = strip_end(d)
            array = array[:fm,:fn]
        # mask transparent values
        if 'transparent_values' in wisdom:
            array = np.ma.masked_array(array,np.logical_and(array >= rng[0], array <= rng[1]))
        else:
            array = np.ma.masked_array(array)

        # scale data
        if scale != 'original':
            m = array.mask.copy()
            a_min, a_max = scale[0], scale[1]
            array[array < a_min] = a_min
            array[array > a_max] = a_max
            array.mask = m

        # all masked data fill with Nans
        array.data[array.mask] = np.nan
        array.fill_value = np.nan
        # set nans to the original No Data Value
        array[np.isnan(array)] = ndv
        # create a driver
        driver = gdal.GetDriverByName('GTiff')
        tiff_path = osp.join(output_path + times[i] + "-" + var + '.tif')
        logging.info('creating tif file: %s' % tiff_path)
        # set up the dataset with zsize bands
        dataset = driver.Create(tiff_path,xsize,ysize,1,gdal.GDT_Float32)
        dataset.SetGeoTransform(geot)
        dataset.SetProjection(projection.ExportToWkt())
        band = dataset.GetRasterBand(1)
        band.WriteArray(array)
        band.SetUnitType(wisdom['native_unit'])
        band.SetDescription(wisdom['name'])
        band.SetNoDataValue(ndv)
        dataset.FlushCache()
        tiffiles.append(tiff_path)
    return tiffiles

def vector2tiffs(output_path, d, wisdom, projection, geot, times, var, ndv=-9999.0):
    '''
    Creates new GeoTiffs for each time from 2D array
    '''
    w,uw,vw = wisdom
    array = uw['retrieve_as'](d, 0)
    if is_fire_var(var):
        fm,fn = strip_end(d)
        array = array[:fm,:fn]
    ysize,xsize = array.shape
    zsize = len(times)
    # write each slice of the array along the zsize
    tiffiles = []
    for i in range(zsize):
        uarray = uw['retrieve_as'](d, i)
        varray = vw['retrieve_as'](d, i)
        if is_fire_var(var):
            fm,fn = strip_end(d)
            uarray = uarray[:fm,:fn]
            varray = varray[:fm,:fn]
        # set nans to the original No Data Value
        uarray[np.isnan(uarray)] = ndv
        varray[np.isnan(varray)] = ndv
        # create a driver
        driver = gdal.GetDriverByName('GTiff')
        tiff_path = osp.join(output_path + times[i] + "-" + var + '.tif')
        logging.info('creating tif file: %s' % tiff_path)
        # set up the dataset with zsize bands
        dataset = driver.Create(tiff_path,xsize,ysize,2,gdal.GDT_Float32)
        dataset.SetGeoTransform(geot)
        dataset.SetProjection(projection.ExportToWkt())
        band1 = dataset.GetRasterBand(1)
        band1.WriteArray(np.flipud(uarray))
        band1.SetUnitType(w['native_unit'])
        band1.SetDescription(uw['name'])
        band1.SetNoDataValue(ndv)
        band2 = dataset.GetRasterBand(2)
        band2.WriteArray(np.flipud(varray))
        band2.SetUnitType(w['native_unit'])
        band2.SetDescription(vw['name'])
        band2.SetNoDataValue(ndv)
        dataset.FlushCache()
        tiffiles.append(tiff_path)
    return tiffiles

def process_vars_tiff(pp, d, wrfout_path, dom_id, times, vars):
    """
    Postprocess a list of scalar or vector fields for a given wrfout file into TIFF files.

    :param pp: Postprocess class
    :param d: the open netCDF file
    :param dom_id: the domain identifier
    :param times: list of times to process
    :param vars: list of variables to process
    """

    logging.info('process_vars_tiff: looking for file %s' % wrfout_path)
    # netCDF WRF metadata
    projection, gt_atm, gt_fire = ncwrfmeta(d)

    outpath_base = osp.join(pp.output_path, pp.product_name + ("-%02d-" % dom_id))
    # build an output file per variable
    for var in vars:
        logging.info('process_vars_tiff: postprocessing %s' % var)
        try:
            tiff_path, coords, mf_upd = None, None, {}
            wisdom = get_wisdom(var).copy()
            wisdom.update(pp.wisdom_update.get(var, {}))
            if is_fire_var(var):
                gt = gt_fire
            else:
                gt = gt_atm

            if is_windvec(var):
                u_name, v_name = wisdom['components']
                uw, vw = get_wisdom(u_name), get_wisdom(v_name)
                uw.update(pp.wisdom_update.get(u_name, {}))
                vw.update(pp.wisdom_update.get(v_name, {}))
                tiff_path = vector2tiffs(outpath_base, d, (wisdom,uw,vw), projection, gt, times, var)
            else:
                tiff_path = scalar2tiffs(outpath_base, d, wisdom, projection, gt, times, var)

            for idx,time in enumerate(times):
                mf_upd['tiff'] = osp.basename(tiff_path[idx])
                ts_esmf = time.replace('_','T')+'Z'
                pp._update_manifest(dom_id, ts_esmf, var, mf_upd)

        except Exception as e:
            logging.warning("Exception %s while postprocessing %s" % (e, var))
            logging.warning(traceback.print_exc())

def process_outputs_tiff(job_id):
    args = load_sys_cfg()
    jobfile = osp.abspath(osp.join(args.workspace_path, job_id,'job.json'))
    satfile = osp.abspath(osp.join(args.workspace_path, job_id,'sat.json'))
    logging.info('process_tiffs: loading job description from %s' % jobfile)
    try:
        js = Dict(json.load(open(jobfile,'r')))
    except Exception as e:
        logging.error('Cannot load the job description file %s' % jobfile)
        logging.error('%s' % e)
        sys.exit(1)
    logging.info('process_tiffs: loading satellite description from %s' % satfile)
    try:
        jsat = Dict(json.load(open(satfile,'r')))
        available_sats = [sat.upper()+prod for sat in jsat.granules.keys() for prod in _sat_prods]
        not_empty_sats = [sat.upper()+prod for sat in jsat.granules.keys() for prod in _sat_prods if jsat.granules[sat]]
    except:
        logging.warning('Cannot load the satellite data in satellite description file %s' % satfile)
        available_sats = []
        not_empty_sats = []
        pass
    logging.info('process_tiffs: available satellite data %s' % available_sats)
    logging.info('process_tiffs: not empty satellite data %s' % not_empty_sats)

    js.old_pid = js.pid
    js.pid = os.getpid()
    js.state = 'Processing'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))
    js.wrf_dir = osp.abspath(osp.join(args.workspace_path, js.job_id, 'wrf'))

    # set up postprocessing
    pp = None
    js.pp_dir = osp.join(args.workspace_path, js.job_id, "products", "tiffs")
    make_clean_dir(js.pp_dir)
    pp = Postprocessor(js.pp_dir, 'wfc-' + js.grid_code)
    js.manifest_filename= 'wfc-' + js.grid_code + '.json'
    logging.debug('Postprocessor created manifest %s',js.manifest_filename)

    logging.info('Postprocessing all wrfout files.')
    # postprocess all wrfouts
    for wrfout_path in sorted(glob.glob(osp.join(js.wrf_dir,'wrfout_d??_????-??-??_??:??:??'))):
        logging.info("Found %s" % wrfout_path)
        domain_str,wrfout_esmf_time = re.match(r'.*wrfout_d(0[0-9])_([0-9_\-:]{19})',wrfout_path).groups()
        dom_id = int(domain_str)
        d = nc4.Dataset(wrfout_path)
        # extract ESMF string times
        times = [''.join(x) for x in d.variables['Times'][:].astype(str)]
        if js.postproc is not None and str(dom_id) in js.postproc:
            if available_sats:
                sat_list = [sat for sat in available_sats if sat in js.postproc[str(dom_id)]]
                var_list = [str(x) for x in js.postproc[str(dom_id)] if not str(x) in sat_list]
                sat_list = [sat for sat in sat_list if sat in not_empty_sats]
                logging.info("Executing postproc instructions for sats %s for domain %d." % (str(sat_list), dom_id))
            else:
                sat_list = []
                var_list = [str(x) for x in js.postproc[str(dom_id)]]
            logging.info("Executing postproc tiff instructions for vars %s for domain %d." % (str(var_list), dom_id))
            try:
                if sat_list:
                    pass
                    #process_sats_tiff()
                process_vars_tiff(pp, d, wrfout_path, dom_id, times, var_list)
            except Exception as e:
                logging.warning('Failed to postprocess with error %s.' % str(e))
        d.close()

    js.old_pid = js.pid
    js.pid = None
    js.state = 'Completed'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

def ncwrfmeta(ds):
    attrs = ds.ncattrs()
    dims = ds.dimensions
    get_attr = lambda attr,exc=None: ds.getncattr(attr) if attr in attrs else exc
    get_dim = lambda dim,exc=None: dims[dim] if dim in dims else exc
    # getting metadata
    lat1 = get_attr('TRUELAT1')
    lat2 = get_attr('TRUELAT2')
    lat0 = get_attr('MOAD_CEN_LAT')
    lon0 = get_attr('STAND_LON')
    clat = get_attr('CEN_LAT')
    clon = get_attr('CEN_LON')
    # creating CSR onject
    crs = osr.SpatialReference()
    proj4 = '+proj=lcc +lat_1=%.10f +lat_2=%.10f +lat_0=%.10f +lon_0=%.10f +a=6370000.0 +b=6370000.0' % (lat1,lat2,lat0,lon0)
    logging.info('ncwrfmeta - proj4=%s' % proj4)
    crs.ImportFromProj4(proj4)
    wrf_proj = pyproj.Proj(proj4)
    ll_proj = wrf_proj.to_latlong()
    # creating atmospheric geotransform
    e,n = pyproj.transform(ll_proj,wrf_proj,clon,clat)
    dx_atm = get_attr('DX')
    dy_atm = get_attr('DY')
    nx_atm = get_dim('west_east').size
    ny_atm = get_dim('south_north').size
    x0_atm = -nx_atm / 2. * dx_atm + e
    y0_atm = -ny_atm / 2. * dy_atm + n
    gt_atm = (x0_atm,dx_atm,0,y0_atm,0,dy_atm)
    logging.info('ncwrfmeta - GT_atm: (%g,%g,%g,%g,%g,%g)' % gt_atm)

    # creating fire geotransform
    if 'west_east_subgrid' in dims:
        nx = get_dim('west_east_subgrid').size
        ny = get_dim('south_north_subgrid').size
        srx = int(nx/(nx_atm+1))
        sry = int(ny/(ny_atm+1))
        if srx > 0 and sry > 0:
            nx_fire = nx - srx
            ny_fire = ny - sry
            dx_fire = dx_atm/srx
            dy_fire = dy_atm/sry
            x0_fire = -nx_fire / 2. * dx_fire + e
            y0_fire = -ny_fire / 2. * dy_fire + n
            gt_fire = (x0_fire,dx_fire,0,y0_fire,0,dy_fire)
            logging.info('ncwrfmeta - GT_fire: (%g,%g,%g,%g,%g,%g)' % gt_fire)
        else:
            gt_fire = None
    else:
        gt_fire = None

    return crs, gt_atm, gt_fire

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit('usage: ./process_output_tiffs.sh job_id')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    process_outputs_tiff(sys.argv[1])
