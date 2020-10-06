#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
import glob, re, datetime, logging, requests, json
import os.path as osp
import numpy as np
from datetime import datetime 
from six.moves.urllib.request import urlopen
from cmr import GranuleQuery
from utils import Dict, esmf_to_utc, utc_to_utcf, duplicates
from .downloader import download_url, DownloadError
from six.moves import filter
from six.moves import range

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
        self.appkey=js.get('appkey',None)
        if not self.appkey:
            try:
                tokens = json.load(open('etc/tokens.json'))
                self.appkey = tokens.get('appkey',None)
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
            if content_size != 0:
                return osp.getsize(path) == content_size
        return False

    def search_api_sat(self, sname, bbox, time, collection=None):
        """
        API search of the different satellite granules and return metadata dictionary

        :param sname: short name satellite product, ex: 'MOD03'
        :param bbox: polygon with the search bounding box
        :param time: time interval (init_time_iso,final_time_iso)

        :return metas: a dictionary with all the metadata for the API search
        """
        maxg=1000
        time_utcf=(utc_to_utcf(time[0]),utc_to_utcf(time[1]))
        api = GranuleQuery()
        search = api.parameters(
                short_name=sname,
                downloadable=True,
                polygon=bbox,
                temporal=time_utcf)
        sh=search.hits()
        if sh>maxg:
            logging.warning("The number of hits %s is larger than the limit %s." % (sh,maxg))
            logging.warning("Any satellite data with prefix %s used." % self.prefix)
            logging.warning("Use a reduced bounding box or a reduced time interval.")
            metas = []
        else:
            metas = api.get(sh)
        if collection:
            metas = [m for m in metas if m['collection_concept_id'] == collection]
        logging.info('search_api_sat - CMR API gives {0} hits for {1} of collection {2}'.format(len(metas),sname,collection))
        return metas

    def archive_url(self, meta, path_col, nrt=False):
        """
        Archive URL construction

        :param meta: metadata from CMR API search
        :param path_col: collection folder for the path
        :param nrt: near real time (nrt) flag
        :return url: url of the archive reconstruction
        """
        base_url = self.base_url if not nrt else self.base_url_nrt
        g_id = osp.basename(meta['links'][0]['href'])
        split = g_id.split('.')
        time_str = '{0}-{1}_{2}:{3}'.format(split[1][1:5],split[1][5:8],split[2][:2],split[2][2:4])
        tt = datetime.strptime(time_str,'%Y-%j_%H:%M').timetuple()
        folder = '{0:04d}/{1:03d}'.format(tt.tm_year,tt.tm_yday)
        url = osp.join(base_url,path_col,folder,g_id) 
        return url

    def get_metas_sat(self, bbox, time):
        """
        Get all the meta data for all the necessary products

        :param bbox: polygon with the search bounding box
        :param time: time interval (init_time_datetime,final_time_datetime)
        :param maxg: max number of granules to process
        :return granules: dictionary with the metadata of all the products
        """
        metas=Dict({})
        metas.geo=self.search_api_sat(self.geo_prefix,bbox,time,collection=self.geo_collection_id)
        metas.fire=self.search_api_sat(self.fire_prefix,bbox,time,collection=self.fire_collection_id)
        metas.geo_nrt=self.search_api_sat(self.geo_nrt_prefix,bbox,time,collection=self.geo_nrt_collection_id)
        metas.fire_nrt=self.search_api_sat(self.fire_nrt_prefix,bbox,time,collection=self.fire_nrt_collection_id)
        return metas

    def group_metas(self, metas):
        """
        Group all the satellite metas before downloading to minimize number of files downloaded

        :param metas: satellite metadatas from API search
        :return result: groupped and cleanned metadata dictionary
        """
        g_id = lambda m: '_'.join(m['producer_granule_id'].split('.')[1:3])
        gmetas = Dict({})
        gmetas.geo = Dict({})
        gmetas.fire = Dict({})
        for m in metas['geo']:
            m.update({'archive_url': self.archive_url(m,osp.join(self.geo_col,self.geo_prefix))})
            gmetas.geo.update({g_id(m): m})
        # add fire metas, if geo available
        for m in metas['fire']:
            k = g_id(m)
            if k in gmetas.geo.keys():
                m.update({'archive_url': self.archive_url(m,osp.join(self.fire_col,self.fire_prefix))})
                gmetas.fire.update({k: m})
            else:
                logging.warning('group_metas - geolocation meta not found for id {}, eliminating fire meta'.format(k))
        # add geolocation NRT metas, if necessary
        for m in metas['geo_nrt']:
            k = g_id(m)
            if k not in gmetas.geo.keys():
                m.update({'archive_url': self.archive_url(m,osp.join(self.geo_nrt_col,self.geo_nrt_prefix),True)})
                gmetas.geo.update({k: m})
        # add fire NRT metas, if necessary
        for m in metas['fire_nrt']:
            k = g_id(m)
            if k not in gmetas.fire.keys():
                if k in gmetas.geo.keys():
                    m.update({'archive_url': self.archive_url(m,osp.join(self.fire_nrt_col,self.fire_nrt_prefix),True)})
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

    def download_sat(self, urls, appkey):
        """
        Download a satellite file from a satellite service

        :param urls: the URL of the file
        :param appkey: key to use for the download or None if not
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
                    download_url(url, sat_path, appkey=appkey)
                    return {'url': url,'local_path': sat_path,'downloaded': datetime.now}
                except DownloadError as e:
                    logging.warning('download_sat - {0} cannot download satellite file {1}'.format(self.prefix, url))

        logging.error('%s cannot download satellite file using %s' % (self.prefix, urls))
        logging.warning('Please check %s for %s' % (self.info_url, self.info))
        return {}
        raise SatError('SatSource: failed to download file %s' % url)

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
                urls = [geo_meta['links'][0]['href'],geo_meta.get('archive_url')]
                m_geo = self.download_sat(urls,self.datacenter_to_appkey(geo_meta['data_center']))
                if m_geo:
                    geo_meta.update(m_geo)
                    urls = [fire_meta['links'][0]['href'],fire_meta.get('archive_url')]
                    m_fire = self.download_sat(urls,self.datacenter_to_appkey(fire_meta['data_center']))
                    if m_fire:
                        fire_meta.update(m_fire)
                        manifest.update({g_id: {
                            'time_start_iso' : geo_meta['time_start'],
                            'time_end_iso' : geo_meta['time_end'],
                            'geo_url' : geo_meta['url'],
                            'geo_local_path' : geo_meta['local_path'],
                            'geo_description' : geo_meta['dataset_id'],
                            'fire_url' : fire_meta['url'],
                            'fire_local_path' : fire_meta['local_path'],
                            'fire_description' : fire_meta['dataset_id']
                        }})
                

        return manifest

    def retrieve_data_sat(self, bounds, from_utc, to_utc):
        """
        Retrieve satellite data in a bounding box coordinates and time interval

        :param bounds: polygon with the search bounding box
        :param time: time interval (init_time_iso,final_time_iso)
        :return data: dictonary with all the data
        """
        if not osp.exists(osp.join(osp.expanduser('~'),'.netrc')):
            logging.warning('satellite acquisition can fail because some data centers require to have $HOME/.netrc specified from an existent Earthdata account')
        
        lonmin,lonmax,latmin,latmax = bounds
        bbox = [(lonmin,latmax),(lonmin,latmin),(lonmax,latmin),(lonmax,latmax),(lonmin,latmax)]
        time = (from_utc, to_utc)

        metas = self.group_metas(self.get_metas_sat(bbox,time))
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

    def datacenter_to_appkey(self,data_center):
        """
        From data center to appkey to use for that data center

        :param data_center: string with the data center information
        """
        return {'LAADS': self.appkey,
                'LPDAAC_ECS': None,
                'LANCEMODIS': self.appkey
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
    base_url='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData'
    base_url_nrt='https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/allData'

