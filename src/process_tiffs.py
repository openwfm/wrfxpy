from osgeo import gdal, osr, gdal_array
import pyproj

import os, logging, json, sys, glob, re, traceback
import os.path as osp
import numpy as np
import netCDF4 as nc4

from vis.postprocessor import Postprocessor
from utils import load_sys_cfg, Dict, make_clean_dir
from vis.var_wisdom import get_wisdom, is_windvec

def scalar2tiffs(output_path, d, wisdom, projection, geot, times, var, ndv=-9999.0):
    '''
    Creates new GeoTiffs for each time from 2D array
    '''
    array = wisdom['retrieve_as'](d,0)
    datatype = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)
    if type(datatype)!=np.int:
        if datatype.startswith('gdal.GDT_')==False:
            datatype=eval('gdal.GDT_'+datatype)
    ysize,xsize = array.shape
    zsize = len(times)
    # write each slice of the array along the zsize
    tiffiles = []
    for i in range(zsize):
        array = wisdom['retrieve_as'](d,i)
        # set nans to the original No Data Value
        array[np.isnan(array)] = ndv
        # create a driver
        driver = gdal.GetDriverByName('GTiff')
        tiff_path = osp.join(output_path + times[i] + "-" + var + '.tif')
        logging.info('creating tif file: %s' % tiff_path)
        # set up the dataset with zsize bands
        dataset = driver.Create(tiff_path,xsize,ysize,1,datatype)
        dataset.SetGeoTransform(geot)
        dataset.SetProjection(projection.ExportToWkt())
        dataset.GetRasterBand(1).WriteArray(np.flipud(array))
        dataset.GetRasterBand(1).SetNoDataValue(ndv)
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
    projection, geotransform = ncwrfmeta(d)

    outpath_base = osp.join(pp.output_path, pp.product_name + ("-%02d-" % dom_id))
    # build an output file per variable
    for var in vars:
        logging.info('process_vars_tiff: postprocessing %s' % var)
        try:
            tiff_path, coords, mf_upd = None, None, {}
            if is_windvec(var):
                continue
                #tiff_path = vector2shps(outpath_base, data, projection, geotransform, times, var)
            else:
                wisdom = get_wisdom(var).copy()
                wisdom.update(pp.wisdom_update.get(var, {}))
                tiff_path = scalar2tiffs(outpath_base, d, wisdom, projection, geotransform, times, var)

            for idx,time in enumerate(times):
                mf_upd['tiff'] = osp.basename(tiff_path[idx])
                ts_esmf = time.replace('_','T')+'Z'
                pp._update_manifest(dom_id, ts_esmf, var, mf_upd)
        except Exception as e:
            logging.warning("Exception %s while postprocessing %s" % (e.message, var))
            logging.warning(traceback.print_exc())

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
        times = [''.join(x.astype(str)) for x in d.variables['Times'][:]]
        if js.postproc is not None and str(dom_id) in js.postproc:
            var_list = [str(x) for x in js.postproc[str(dom_id)]]
            logging.info("Executing postproc tiff instructions for vars %s for domain %d." % (str(var_list), dom_id))
            try:
                process_vars_tiff(pp, d, wrfout_path, dom_id, times, var_list)
            except Exception as e:
                logging.warning('Failed to postprocess for time %s with error %s.' % (esmf_time, str(e)))
        d.close()

    js.old_pid = js.pid
    js.pid = None
    js.state = 'Completed'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

def ncwrfmeta(d):
    # getting metadata
    lat1 = d.TRUELAT1
    lat2 = d.TRUELAT2
    lat0 = d.MOAD_CEN_LAT
    lon0 = d.STAND_LON
    clat = d.CEN_LAT
    clon = d.CEN_LON
    csr = osr.SpatialReference()
    proj4 = '+proj=lcc +lat_1=%.10f +lat_2=%.10f +lat_0=%.10f +lon_0=%.10f +a=6370000.0 +b=6370000.0' % (lat1,lat2,lat0,lon0)
    logging.info('proj4: %s' % proj4)
    csr.ImportFromProj4(proj4)
    ll_proj = pyproj.Proj('+proj=latlong +datum=WGS84')
    wrf_proj = pyproj.Proj(proj4)
    # geotransform
    dx = d.DX
    dy = d.DY
    nx = d.dimensions['west_east'].size
    ny = d.dimensions['south_north'].size
    if 'west_east_subgrid' in d.dimensions:
        dx_atm,dy_atm,nx_atm,ny_atm = dx,dy,nx,ny
        nx = d.dimensions['west_east_subgrid'].size
        ny = d.dimensions['south_north_subgrid'].size
        srx = int(nx/(nx_atm+1))
        sry = int(ny/(ny_atm+1))
        dx = dx_atm/srx
        dy = dy_atm/sry
    e,n = pyproj.transform(ll_proj,wrf_proj,clon,clat)
    x0 = -nx / 2. * dx + e
    y1 = ny / 2. * dy + n
    geotransform = (x0,dx,0,y1,0,-dy)
    logging.info('geotransform: (%g,%g,%g,%g,%g,%g)' % geotransform)
    return csr, geotransform

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit('usage: ./process_output_tiffs.sh job_id')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    process_outputs_tiff(sys.argv[1])
