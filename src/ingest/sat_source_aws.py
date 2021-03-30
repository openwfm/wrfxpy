#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
import logging
import re
import glob
import json
import subprocess
import os.path as osp
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
from utils import Dict, split_path, utc_to_utcf, ensure_dir
from ingest.sat_source import SatSource

def parse_filename(file_name):
    """
    Parse filename from AWS dataset name

    :param file_name: AWS file name to parse
    :return info: dictionary with information from file name
    """
    info = {}
    pattern = 'FDC([CF]{1})-M([0-9]{1})_G([0-9]{2})_s([0-9]{14})_e([0-9]{14})'
    tfrmt = '%Y%j%H%M%S'
    match = re.search(pattern,file_name)
    if len(match.groups()) == 5:
        info['domain'] = match.group(1)
        info['mode'] = match.group(2)
        info['satellite'] = match.group(3)
        info['start_date'] = datetime.strptime(match.group(4)[0:13],tfrmt)
        info['end_date'] = datetime.strptime(match.group(5)[0:13],tfrmt)
    return info

def parse_projinfo(file_name):
    """
    Parse projection information from NetCDF AWS dataset file

    :param file_name: AWS file name to parse
    :return info: dictionary with information from file name
    """
    d = nc.Dataset(file_name)
    gip = d['goes_imager_projection']
    x = d['x']
    y = d['y']
    return {'shgt': gip.getncattr('perspective_point_height'), 'requ': gip.getncattr('semi_major_axis'),
        'rpol': gip.getncattr('semi_minor_axis'), 'lon0': gip.getncattr('longitude_of_projection_origin'),
        'xscl': x.getncattr('scale_factor'), 'xoff': x.getncattr('add_offset'), 'xdim': x.shape[0], 
        'yscl': y.getncattr('scale_factor'), 'yoff': y.getncattr('add_offset'), 'ydim': y.shape[0],
        'file_name': file_name}

def aws_request(cmd, max_retries=5):
    """
    Request command in AWS storage using AWS CLI retrying if an error occurs

    :param cmd: list with command and flags to run using subprocess
    :param max_retries: max retries to request to AWS CLI
    :return r: output from check_output method of subprocess
    """
    for tries in range(max_retries):
        try:
            logging.debug('aws_request - {}'.format(' '.join(cmd)))
            r = subprocess.check_output(cmd)
            return r
        except:
            if tries == max_retries-1:
                logging.error('error when requesting "{}", no more retries left'.format(' '.join(cmd)))
            else:
                logging.warning('error when requesting "{}", retry {}/{}...'.format(' '.join(cmd),tries+1,max_retries))
    return None

def aws_search(awspaths, time=(datetime(2000,1,1),datetime.now())):
    """
    Search data in AWS storage using AWS CLI

    :param awspaths: list of paths to search using AWS CLI
    :param time: time interval in datetime format to consider the data from
    :return result: metadata of the satellite data found
    """
    ls_cmd = 'aws s3 ls {} --recursive --no-sign-request'
    result = []
    for awspath in awspaths:
        cmd = ls_cmd.format(awspath).split()
        r = aws_request(cmd)
        for line in r.decode().split('\n'):
            if len(line):
                _,_,file_size,file_name = list(map(str.strip,line.split()))
                info = parse_filename(file_name)
                if info['start_date'] <= time[1] and info['start_date'] >= time[0]:
                    base = split_path(awspath)[0]
                    url = osp.join('s3://',base,file_name) 
                    file_basename = osp.basename(file_name)
                    result.append({'file_name': file_basename, 'file_remote_path': file_name,
                        'time_start': utc_to_utcf(info['start_date']), 
                        'time_end': utc_to_utcf(info['end_date']), 
                        'updated': datetime.now(), 'url': url, 'file_size': int(file_size), 
                        'domain': info['domain'], 'satellite': 'G{}'.format(info['satellite']), 
                        'mode': info['mode'], 
                        'granule_id': 'A{:04d}{:03d}{:02d}{:02d}'.format(
                            info['start_date'].year,
                            info['start_date'].timetuple().tm_yday,
                            info['start_date'].hour,
                            info['start_date'].minute
                        )})
    return result

def download_url(url, sat_path):
    """
    Download URL from AWS storage using AWS CLI

    :param url:
    :param sat_path:
    """
    cmd = 'aws s3 cp {} {} --no-sign-request'.format(url, sat_path).split()
    r = aws_request(cmd)

def create_grid(proj_info, out_path):
    """
    Create AWS grid

    :param proj_info: projection information dictionary to create the grid coordinates
    :param out_path: output path to write the grid coordinates
    """
    shgt = proj_info['shgt']
    requ = proj_info['requ']
    rpol = proj_info['rpol']
    lon0 = proj_info['lon0']
    xscl = proj_info['xscl']
    xoff = proj_info['xoff']
    xdim = proj_info['xdim']
    yscl = proj_info['yscl']
    yoff = proj_info['yoff']
    ydim = proj_info['ydim']
    dtor = 0.0174533
    pi = 3.1415926536
    rrat = (requ*requ)/(rpol*rpol)
    bigh = requ + shgt
    i = np.reshape(np.arange(xdim),(1,xdim))
    j = np.reshape(np.arange(ydim),(ydim,1))[::-1]
    x = i*xscl + xoff
    y = j*yscl + yoff
    a = np.sin(x)**2 + np.cos(x)**2*(np.cos(y)**2 + rrat*np.sin(y)**2)
    b = -2*bigh*np.cos(x)*np.cos(y)
    c = bigh**2-requ**2
    with np.errstate(invalid='ignore'):
        rs = (-b - np.sqrt(b*b - 4.*a*c))/(2.*a)
    sx = rs*np.cos(x)*np.cos(y)
    sy = -rs*np.sin(x)
    sz = rs*np.cos(x)*np.sin(y)
    lats = (np.arctan(rrat*sz/np.sqrt((bigh-sx)**2 + sy**2)))/dtor
    lons = (lon0*dtor - np.arctan(sy/(bigh-sx)))/dtor
    nc_file = nc.Dataset(out_path, 'w')
    xx = nc_file.createDimension("x", xdim)
    yy = nc_file.createDimension("y", ydim)
    longitude = nc_file.createVariable("lon", "f4", ("y","x"))
    latitude = nc_file.createVariable("lat", "f4", ("y","x"))
    longitude[:] = lons
    latitude[:] = lats
    nc_file.close()

def process_grid(ingest_dir, meta):
    """
    Process AWS grid

    :param ingest_dir: path to ingest directory for the satellite product
    :param meta: metadata to create the grid for
    """
    current_proj_path = osp.join(ingest_dir,'{}_projinfo.json'.format(meta['file_name'].split('.')[0]))
    current_grid_path = osp.join(ingest_dir,'{}_grid.nc'.format(meta['file_name'].split('.')[0]))
    current_proj = parse_projinfo(meta['local_path'])
    current_proj.update({'grid_path': current_grid_path})
    archived_proj_paths = glob.glob(osp.join(ingest_dir, '*_projinfo.json'))
    if len(archived_proj_paths):
        for archived_proj_path in archived_proj_paths:
            if osp.exists(archived_proj_path):
                archived_proj = eval(json.load(open(archived_proj_path, 'r')))
                if archived_proj == current_proj:
                    archived_grid_path = archived_proj['grid_path']
                    if not osp.exists(archived_grid_path):
                        create_grid(archived_proj, archived_grid_path)
                    return {'proj_path': archived_proj_path, 'grid_path': archived_grid_path}
    json.dump(str(current_proj), open(current_proj_path, 'w'))
    create_grid(current_proj, current_grid_path)
    return {'proj_path': current_proj_path, 'grid_path': current_grid_path}

class SatSourceAWSError(Exception):
    """
    Raised when a SatSourceAWS cannot retrieve satellite data.
    """
    pass

class SatSourceAWS(SatSource):
    """
    The parent class of all satellite sources that implements common functionality, for example

    - local satellite validation (file size check)
    - Satellite retrieval with retries (smart check whether server implements http-range)
    """

    def __init__(self, js):
        super(SatSourceAWS, self).__init__(js)

    def search_api_sat(self, sname, bounds, time, collection=None):
        """
        API search of the different satellite granules and return metadata dictionary

        :param sname: short name satellite product, ex: ''
        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_iso,final_time_iso)
        :param collection: id of the collection to specify

        :return metas: a dictionary with all the metadata for the API search
        """
        # complete_days
        stime = time[0]
        etime = time[1]
        n_complete_days = int((etime-stime).days)+1
        complete_days = [(stime + timedelta(days=x)).strftime('%Y/%j/') for x in range(n_complete_days)]
        awspaths = [osp.join(self.base_url.format(self.platform), self.product.format(self.sector), hour) for hour in complete_days]
        metas = aws_search(awspaths,time)
        logging.info('SatSourceAWS.search_api_sat - {} metas found'.format(len(metas)))
        return metas

    def get_metas_sat(self, bounds, time):
        """
        Get all the meta data for all the necessary products

        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_datetime,final_time_datetime)
        :param maxg: max number of granules to process
        :return granules: dictionary with the metadata of all the products
        """
        metas = Dict(self.search_api_sat(None,None,time))
        return metas

    def group_metas(self, metas):
        """
        Group all the satellite metas before downloading to minimize number of files downloaded

        :param metas: satellite metadatas from API search
        :return metas: groupped and cleanned metadata dictionary
        """
        return metas

    def _download_url(self, url, sat_path, appkey):
        """
        Download a satellite file from a satellite service

        :param url: the URL of the file
        :param sat_path: local path to download the file
        :param appkey: key to use for the download or None if not
        """
        download_url(url, sat_path)

    def retrieve_metas(self, metas):
        """
        Retrieve satellite data from CMR API metadata 

        :return metas: dictonary with all the satellite data to retrieve
        :return manifest: dictonary with all the satellite data retrieved
        """
        logging.info('retrieve_metas - downloading {} products'.format(self.id))
        manifest = Dict({})
        for meta in metas:
            urls = [meta['url']]
            fire_meta = self.download_sat(urls)
            fire_meta.update(meta)
            local_path = fire_meta['local_path']
            remote_size = int(fire_meta['file_size'])
            local_size = osp.getsize(local_path)
            if local_size != remote_size:
                logging.error('retireve_metas - local file with size {} different than remote size {}'.format())
                continue
            info_path = local_path + '.size'
            open(ensure_dir(info_path), 'w').write(str(local_size))
            geo_meta = process_grid(self.ingest_dir,fire_meta)
            manifest.update({fire_meta['granule_id']: {
                'time_start_iso' : fire_meta['time_start'],
                'time_end_iso' : fire_meta['time_end'],
                'geo_url' : fire_meta['url'],
                'geo_local_path' : geo_meta['grid_path'],
                'geo_description' : self.info,
                'fire_url' : fire_meta['url'],
                'fire_local_path' : fire_meta['local_path'],
                'fire_description' : self.info
            }})
        return manifest

    def read_sat(self, files, metas, bounds):
        """
        Read all the satellite files from a specific satellite service

        :param file: dictionary with geolocation (03), fire (14) (, and reflectance (09)) file paths for satellite service
        :param bounds: spatial bounds tuple (lonmin,lonmax,latmin,latmax)
        :return ret: dictionary with Latitude, Longitude and fire mask arrays read
        """
        pass

    def read_files_sat(self, file, bounds):
        """
        Read a specific satellite files from a specific satellite service

        :param file: dictionary with geolocation (03), fire (14) (, and reflectance (09)) file paths for satellite service
        :param bounds: spatial bounds tuple (lonmin,lonmax,latmin,latmax)
        :return ret: dictionary with Latitude, Longitude and fire mask arrays read
        """
        pass

    # instance variables
    base_url='noaa-{}'
    product='ABI-L2-FDC{}'
    sector='F' # full disk

