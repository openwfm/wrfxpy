#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
import logging
import os.path as osp
import re
import requests
import pytz
from datetime import datetime,timedelta
from utils import Dict, utc_to_utcf
from ingest.sat_source import SatSource

def parse_filename(file_name):
    """
    Parse filename from archive dataset name

    :param file_name: archive file name to parse
    :return info: dictionary with information from file name
    """
    info = {}
    tfrmt = '%Y%m%d%H%M%S%f'
    pattern = '_d([0-9]{8})_t([0-9]{7})_e([0-9]{7})_b([0-9]{5})_'
    match = re.search(pattern,file_name)
    if len(match.groups()) == 4:
        start_date = match.group(1)+match.group(2)+'00000' 
        end_date = match.group(1)+match.group(3)+'00000' 
        info['start_date'] = datetime.strptime(start_date,tfrmt).replace(tzinfo=pytz.UTC)
        info['end_date'] = datetime.strptime(end_date,tfrmt).replace(tzinfo=pytz.UTC)
        info['orbit_number'] = match.group(4)
        info['producer_granule_id'] = 'd{:08d}_t{:07d}_e{:07d}_b{:05d}'.format(*map(int,match.groups()))
    return info

def archive_request(url, max_retries=5):
    """
    Listing files in archive URL and retrying if an error occurs

    :param url: url string
    :param max_retries: max retries to request
    :return files: list files in the URL
    """
    pattern = "a href='([^']+)'"
    files = []
    for tries in range(max_retries):
        try:
            r = requests.get(url)
            if r.status_code != 200:
                logging.warning('status code: {}'.format(r.status_code))
            else:
                content = r.content.decode('utf-8')
                for aline in content.split('\n'):
                    amatch = re.search(pattern,aline)
                    if amatch:
                       files.append(osp.basename(amatch.group(1))) 
                return files
        except Exception as e:
            if tries == max_retries-1:
                logging.error('error when requesting {}, no more retries left'.format(url))
                break
            else:
                logging.warning('error when requesting {}, retry {}/{}'.format(url,tries+1,max_retries))
    return files

def archive_search(archive_paths, time):
    """
    Search files in the path

    :param archive_paths: list of archive paths to construct the manifest from
    :return metas: metadata of the satellite data found
    """
    metas = {}
    for url in archive_paths:
       files = archive_request(url) 
       for f in files:
           info = parse_filename(f)
           if info['start_date'] <= time[1].replace(tzinfo=pytz.UTC) and info['start_date'] >= time[0].replace(tzinfo=pytz.UTC):
               product_id = 'A{:04d}{:03d}_{:02d}{:02d}{:02d}'.format(info['start_date'].year,
                                            info['start_date'].timetuple().tm_yday,
                                            info['start_date'].hour,
                                            info['start_date'].minute,
                                            info['start_date'].second)
               remote_path = osp.join(url,f)
               metas.update({product_id: {'file_name': f, 'file_remote_path': remote_path, 
                   'time_start': utc_to_utcf(info['start_date']), 
                   'time_end':  utc_to_utcf(info['end_date']),
                   'updated': datetime.now(), 'orbit_number': info['orbit_number'],
                   'url': remote_path, 'product_id': product_id}})

    return metas

class SatSourceArchiveError(Exception):
    """
    Raised when a SatSourceArchive cannot retrieve satellite data.
    """
    pass

class SatSourceArchive(SatSource):
    """
    The parent class of all satellite sources from archive that implements common functionality, for example

    - local satellite validation (file size check)
    - Satellite retrieval with retries (smart check whether server implements http-range)
    """

    def __init__(self, js):
        super(SatSourceArchive, self).__init__(js)

    def search_api_sat(self, sname, bounds, time, collection=None):
        """
        Archive URL construction

        :param sname: short name satellite product, ex: ''
        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_iso,final_time_iso)
        :param collection: id of the collection to specify
        """
        nrt = False if collection is None else collection
        base_url = self.base_url if not nrt else self.base_url_nrt
        stime,etime = time
        n_complete_days = int((etime.replace(hour=0,minute=0,second=0,microsecond=0)-stime.replace(hour=0,minute=0,second=0,microsecond=0)).days)+1
        complete_days = [(stime + timedelta(days=x)).strftime('%Y/%j/') for x in range(n_complete_days)]
        archive_paths = [osp.join(self.base_url, sname, day) for day in complete_days]
        metas = archive_search(archive_paths,time)
        return metas

    def get_metas_sat(self, bounds, time):
        """
        Get all the meta data for all the necessary products

        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_datetime,final_time_datetime)
        :return granules: dictionary with the metadata of all the products
        """
        metas=Dict({})
        if hasattr(self,'geo_prefix') and self.geo_prefix is not None:
            metas.geo=self.search_api_sat(self.geo_prefix,bounds,time,collection=False)
            metas.fire=self.search_api_sat(self.fire_prefix,bounds,time,collection=False)
        else:
            metas.geo={}
            metas.fire={}
        if hasattr(self,'geo_nrt_prefix') and self.geo_nrt_prefix is not None:
            metas.geo_nrt=self.search_api_sat(self.geo_nrt_prefix,bounds,time,collection=True)
            metas.fire_nrt=self.search_api_sat(self.fire_nrt_prefix,bounds,time,collection=True)
        else:
            metas.geo_nrt={}
            metas.fire_nrt={}
        return metas

    def group_metas(self, metas):
        """
        Group all the satellite metas before downloading to minimize number of files downloaded

        :param metas: satellite metadatas from API search
        :return result: groupped and cleanned metadata dictionary
        """
        gmetas = Dict({})
        gmetas.geo = Dict({})
        gmetas.fire = Dict({})
        for k,m in metas['geo'].items():
            gmetas.geo.update({k: m})
        # add fire metas, if geo available
        for k,m in metas['fire'].items():
            if k in gmetas.geo.keys():
                gmetas.fire.update({k: m})
            else:
                logging.warning('group_metas - geolocation meta not found for id {}, eliminating fire meta'.format(k))
        # add geolocation NRT metas, if necessary
        for k,m in metas['geo_nrt'].items():
            if k not in gmetas.geo.keys():
                gmetas.geo.update({k: m})
        # add fire NRT metas, if necessary
        for k,m in metas['fire_nrt'].items():
            if k not in gmetas.fire.keys():
                if k in gmetas.geo.keys():
                    gmetas.fire.update({k: m})
                else:
                    logging.warning('group_metas - geolocation not found for id {}'.format(k))
        # delete geolocation if not fire on it
        exc = []
        for k in gmetas.geo.keys():
            if k not in gmetas.fire.keys():
                logging.warning('group_metas - fire meta not found for id {}, eliminating geolocation meta'.format(k))
                exc.append(k)
        for k in exc:
            gmetas.geo.pop(k)

        return gmetas

    def retrieve_metas(self, metas):
        """
        Retrieve satellite data from CMR API metadata 

        :return metas: dictonary with all the satellite data to retrieve
        :return manifest: dictonary with all the satellite data retrieved
        """
        logging.info('retrieve_metas - downloading {} products'.format(self.id))
        manifest = Dict({})
        for g_id, geo_meta in metas['geo'].items():
            if g_id in metas['fire'].keys():
                fire_meta = metas['fire'][g_id]
                logging.info('retrieve_metas - downloading product id {}'.format(g_id))
                urls = [geo_meta['url']]
                m_geo = self.download_sat(urls)
                if m_geo:
                    geo_meta.update(m_geo)
                    urls = [fire_meta['url']]
                    m_fire = self.download_sat(urls)
                    if m_fire:
                        fire_meta.update(m_fire)
                        manifest.update({g_id: {
                            'time_start_iso' : geo_meta['time_start'],
                            'time_end_iso' : geo_meta['time_end'],
                            'geo_url' : geo_meta['url'],
                            'geo_local_path' : geo_meta['local_path'],
                            'fire_url' : fire_meta['url'],
                            'fire_local_path' : fire_meta['local_path']
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
    base_url=None
    base_url_nrt=None

