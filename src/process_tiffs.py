from osgeo import gdal, osr, gdal_array
import pyproj

import os, logging, json, sys, glob, re, traceback
import os.path as osp
import numpy as np
import netCDF4 as nc4

from utils import load_sys_cfg, Dict, make_clean_dir
from vis.var_wisdom import get_wisdom, is_windvec

def scalar2tiffs(output_path, array, projection, geot, times, var, ndv=-9999.0):
    '''
    Creates new GeoTiffs for each time from 2D array
    '''
    datatype = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)

    if type(datatype)!=np.int:
        if datatype.startswith('gdal.GDT_')==False:
            datatype=eval('gdal.GDT_'+datatype)

    zsize,ysize,xsize = array.shape
    print('dimensions of %s: ' % suffix, xsize, ysize, zsize)
    # set nans to the original No Data Value
    array[np.isnan(array)] = ndv
    # write each slice of the array along the zsize
    tiffiles = []
    for i in range(zsize):
        # create a driver
        driver = gdal.GetDriverByName('GTiff')
        ts_esmf = times[i].replace('_','T')+'Z'
        tiff_path = osp.join(output_path, ts_esmf + "-" + var, '.tif')
        print('> creating tif file: %s' % tiff_path)
        # set up the dataset with zsize bands
        dataset = driver.Create(tiff_path,xsize,ysize,1,datatype)
        dataset.SetGeoTransform(geot)
        dataset.SetProjection(projection.ExportToWkt())
        dataset.GetRasterBand(1).WriteArray(np.flipud(array[i]))
        dataset.GetRasterBand(1).SetNoDataValue(ndv)
        dataset.FlushCache()
        tiffiles.append(tiff_path)
    return tiffiles

def process_vars_tiff(pp, d, dom_id, times, vars):
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
    projection, geotransform = ncwrfmeta(d)

    outpath_base = osp.join(pp.output_path, pp.product_name + ("-%02d-" % dom_id))
    # build an output file per variable
    for var in vars:
        logging.info('process_vars_tiff: postprocessing %s for time %s' % (var, ts_esmf))
        try:
            tiff_path, coords, mf_upd = None, None, {}
            if is_windvec(var):
                tiff_path = vector2tiff(outpath_base, data, projection, geotransform, times, var)
            else:
                wisdom = get_wisdom(var).copy()
                wisdom.update(pp.wisdom_update.get(var, {}))
                data = wisdom['retrieve_as'](d,range(len(times)))
                tiff_path = scalar2tiffs(outpath_base, data, projection, geotransform, times, var)

            for idx,time in enumerate(times):
                mf_upd['tiff'] = osp.basename(tiff_path[idx])
                ts_esmf = time.replace('_','T')+'Z'
                pp._update_manifest(dom_id, ts_esmf, var, mf_upd)
        except Exception as e:
            logging.warning("Exception %s while postprocessing %s for time %s" % (e.message, var, ts_esmf))
            logging.warning(traceback.print_exc())

    d.close()

def process_outputs_tiff(job_id):
    args = load_sys_cfg()
    jobfile = osp.abspath(osp.join(args.workspace_path, job_id,'job.json'))
    logging.info('process_tiffs: loading job description from %s' % jobfile)
    try:
        js = Dict(json.load(open(jobfile,'r')))
    except Exception as e:
        logging.error('Cannot load the job description file %s' % jobfile)
        logging.error('%s' % e)
        sys.exit(1)
    js.old_pid = js.pid
    js.pid = os.getpid()
    js.state = 'Processing'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))
    js.wrf_dir = osp.abspath(osp.join(args.workspace_path, js.job_id, 'wrf'))

    # set up postprocessing
    pp = None
    js.pp_dir = osp.join(args.workspace_path, js.job_id, "products", "tiffs")
    make_clean_dir(js.pp_dir)
    pp = Postprocessor(js.pp_dir, 'wfc-' + js.grid_code + '-tiffs')
    js.manifest_filename= 'wfc-' + js.grid_code + '-tiffs' + '.json'
    logging.debug('Postprocessor created manifest %s',js.manifest_filename)

    logging.info('Postprocessing all wrfout files.')
    # postprocess all wrfouts
    for wrfout_path in sorted(glob.glob(osp.join(js.wrf_dir,'wrfout_d??_????-??-??_??:??:??'))):
        logging.info("Found %s" % wrfout_path)
        domain_str,wrfout_esmf_time = re.match(r'.*wrfout_d(0[0-9])_([0-9_\-:]{19})',wrfout_path).groups()
        dom_id = int(domain_str)
        d = nc4.Dataset(wrfout_path)
        # extract ESMF string times
        times = [''.join(x) for x in d.variables['Times'][:]]
        if js.postproc is not None and str(dom_id) in js.postproc:
            var_list = [str(x) for x in js.postproc[str(dom_id)]]
            logging.info("Executing postproc tiff instructions for vars %s for domain %d." % (str(var_list), dom_id))
            try:
                process_vars_tiff(pp, d, dom_id, times, var_list)
            except Exception as e:
                logging.warning('Failed to postprocess for time %s with error %s.' % (esmf_time, str(e)))
            for esmf_time in sorted(times):
                logging.info("Saving manifest for domain %d for time %s." % (dom_id, esmf_time))
        d.close()

    js.old_pid = js.pid
    js.pid = None
    js.state = 'Completed'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

def ncwrfmeta(d):
    # projection
    lat1 = d.TRUELAT1
    lat2 = d.TRUELAT2
    lat0 = d.MOAD_CEN_LAT # MOAD_CEN_LAT or CEN_LAT? CEN_LAT is from nest
    lon0 = d.STAND_LON # STAND_LON or CEN_LON? CEN_LON is from nest
    csr = osr.SpatialReference()
    proj4 = '+proj=lcc +lat_1=%.10f +lat_2=%.10f +lat_0=%.10f +lon_0=%.10f +a=6370000.0 +b=6370000.0' % (lat1,lat2,lat0,lon0)
    print('proj4: %s' % proj4)
    csr.ImportFromProj4(proj4)
    # geotransform
    dx_atm = nc.DX
    dy_atm = nc.DY
    nx_atm = nc.dimensions['west_east_stag'].size
    ny_atm = nc.dimensions['south_north_stag'].size
    nx = nc.dimensions['west_east_subgrid'].size
    ny = nc.dimensions['south_north_subgrid'].size
    srx = int(nx/nx_atm)
    sry = int(ny/ny_atm)
    dx = dx_atm/srx
    dy = dy_atm/sry
    p = pyproj.Proj(proj4)
    top_left_lon = nc['FXLONG'][0,-1,0]
    top_left_lat = nc['FXLAT'][0,-1,0]
    top_left = p(top_left_lon,top_left_lat)
    geotransform = (top_left[0],dx,0,top_left[1],0,-dy)
    print('geotransform: ',geotransform)

    return csr, geotransform

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit('usage: ./process_output_tiffs.sh job_id')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    process_outputs_tiff(sys.argv[1])
