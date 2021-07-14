#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
import logging
import json
import os.path as osp
from datetime import datetime 
from utils import Dict
from .downloader import download_url, DownloadError

class SatError(Exception):
    """
    Raised when a SatSource cannot retrieve satellite data.
    """
    pass

class SatSource(object):
    """
    The parent class of all satellite sources that implements common functionality, for example

    - local satellite validation (file size check)
    - Satellite retrieval with retries (smart check whether server implements http-range)
    """

    def __init__(self, js):
        """
        Initialize satellite source with ingest directory (where satellite files are stored).

        :param js: job structure with at least ingest_path root of satellite storage and sys_install_path
        """
        self.ingest_dir=osp.abspath(osp.join(js.get('ingest_path','ingest'),self.prefix))
        self.cache_dir=osp.abspath(js.get('cache_path','cache'))
        self.sys_dir=osp.abspath(js.get('sys_install_path',None))
        self.tokens=js.get('tokens',{})
        if not self.tokens:
            try:
                tokens = json.load(open('etc/tokens.json'))
                self.tokens = tokens.get('tokens',{})
            except:
                logging.warning('Any etc/tokens.json specified, any token is going to be used.')


    def available_locally_sat(self, path):
        """
        Check if a satellite file is available locally and if it's file size checks out.

        :param path: the satellite file path
        """
        info_path = path + '.size'
        if osp.exists(path) and osp.exists(info_path):
            content_size = int(open(info_path).read())
            if content_size > 0:
                return osp.getsize(path) == content_size
        return False

    def search_api_sat(self, sname, bounds, time, collection=None):
        """
        API search of the different satellite granules and return metadata dictionary

        :param sname: short name satellite product, ex: 'MOD03'
        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_iso,final_time_iso)
        :param collection: id of the collection to specify

        :return metas: a dictionary with all the metadata for the API search
        """
        return {}

    def archive_url(self, meta, path_col, nrt=False):
        """
        Archive URL construction

        :param meta: metadata from CMR API search
        :param path_col: collection folder for the path
        :param nrt: near real time (nrt) flag
        :return url: url of the archive reconstruction
        """
        return ''

    def get_metas_sat(self, bounds, time):
        """
        Get all the meta data for all the necessary products

        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_datetime,final_time_datetime)
        :param maxg: max number of granules to process
        :return granules: dictionary with the metadata of all the products
        """
        return Dict({})

    def group_metas(self, metas):
        """
        Group all the satellite metas before downloading to minimize number of files downloaded

        :param metas: satellite metadatas from API search
        :return result: groupped and cleanned metadata dictionary
        """
        return Dict({})

    def _download_url(self, url, sat_path, token=None, min_size=1):
        """
        Download a satellite file from a satellite service

        :param url: the URL of the file
        :param sat_path: local path to download the file
        :param token: key to use for the download or None if not
        """
        download_url(url, sat_path, token=token, min_size=min_size)

    def download_sat(self, urls, token=None, min_size=10000):
        """
        Download all satellite file from a satellite service

        :param urls: the URL of the file
        :param token: key to use for the download or None if not
        """
        for url in urls:
            logging.info('downloading %s satellite data from %s' % (self.prefix, url))
            sat_name = osp.basename(url)
            sat_path = osp.join(self.ingest_dir,sat_name)
            if self.available_locally_sat(sat_path):
                logging.info('%s is available locally' % sat_path)
                return {'url': urls[0],'local_path': sat_path}
            else:
                try:
                    self._download_url(url, sat_path, token=token, min_size=min_size)
                    return {'url': url,'local_path': sat_path,'downloaded': datetime.now()}
                except Exception as e:
                    logging.warning('download_sat - {0} cannot download satellite file {1}'.format(self.prefix, url))
                    logging.warning(e)

        logging.error('%s cannot download satellite file using %s' % (self.prefix, urls))
        logging.warning('Please check %s for %s' % (self.info_url, self.info))
        return {}

    def retrieve_metas(self, metas):
        """
        Retrieve satellite data from CMR API metadata 

        :return metas: dictonary with all the satellite data to retrieve
        :return manifest: dictonary with all the satellite data retrieved
        """
        return Dict({})

    def retrieve_data_sat(self, bounds, from_utc, to_utc):
        """
        Retrieve satellite data in a bounding box coordinates and time interval

        :param bounds: polygon with the search bounding box
        :param time: time interval (init_time_iso,final_time_iso)
        :return data: dictonary with all the data
        """
        if not osp.exists(osp.join(osp.expanduser('~'),'.netrc')):
            logging.warning('satellite acquisition can fail because some data centers require to have $HOME/.netrc specified from an existent Earthdata account')
        
        time = (from_utc, to_utc)
        metas = self.group_metas(self.get_metas_sat(bounds,time))
        logging.info('retrieve_data_sat - found {0} metas for {1} satellite service'.format(sum([len(m) for m in metas.values()]),self.prefix))

        manifest = self.retrieve_metas(metas)
        logging.info('retrieve_data_sat - adquiered {0} granules for {1} satellite service'.format(len(manifest.keys()),self.prefix))

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

    def datacenter_to_token(self,data_center):
        """
        From data center to token to use for that data center

        :param data_center: string with the data center information
        """
        return {'LAADS': self.tokens.get('laads',None),
                'LPDAAC_ECS': None,
                'LANCEMODIS': self.tokens.get('nrt',None)
            }.get(data_center,None)

    # instance variables
    id=None
    platform=None
    info_url=None
    info=None
    prefix=None
    geo_prefix=None
    fire_prefix=None
    geo_nrt_prefix=None
    fire_nrt_prefix=None
    geo_collection_id=None
    fire_collection_id=None
    geo_nrt_collection_id=None
    fire_nrt_collection_id=None
    geo_col=None
    fire_col=None
    geo_nrt_col=None
    fire_nrt_col=None
    base_url=None
    base_url_nrt=None
    product=None
    sector=None

