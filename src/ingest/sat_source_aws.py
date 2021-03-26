#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
import logging
import re
import os.path as osp
from datetime import datetime, timedelta
from utils import Dict, split_path, utc_to_utcf
from ingest.sat_source import SatSource

def parse_filename(file_name):
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

def aws_search(awspaths):
    result = []
    for awspath in awspaths:
        cmd = 'aws s3 ls {} --recursive --no-sign-request'.format(awspath).split()
        logging.debug('aws_search - {}'.format(' '.join(cmd)))
        max_retries = 5
        while max_retries:
            try:
                r = subprocess.check_output(cmd)
                for line in r.decode().split('\n'):
                    if len(line):
                        file_date,file_time,file_size,file_name = list(map(str.strip,line.split()))
                        file_info = parse_filename(file_name)
                        base = split_path(file_name)[0]
                        url = 's3://{}'.format(osp.join(base,file_name)) 
                        result.append({'producer_granule_id': file_name, 
                            'time_start': utc_to_utcf(file_info['start_date']), 
                            'updated': datetime.now(), 'url': url, 
                            'time_end': utc_to_utcf(file_info['end_date']), 
                            'domain': info['domain'], 'satellite': 'G{}'.format(info['satellite']),
                            'mode': info['mode']})
                        break
            except:
                if not max_retries:
                    logging.error('error when requesting "{}", no more retries left'.format(' '.join(cmd)))
                else:
                    logging.warning('error when requesting "{}", retrying...'.format(' '.join(cmd)))
                    max_retries -= 1
    return result

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
        metas = []
        # hours first day
        stime = time[0]
        etime = min(time[0].replace(hour=23),time[1])
        n_first_hours = int((etime-stime).total_seconds()/3600)+1
        first_hours = [(stime + timedelta(hours=x)).strftime('%Y/%j/%H/') for x in range(n_first_hours)]
        awspaths = [osp.join(self.base_url.format(self.platform), self.product.format(self.sector), hour) for hour in first_hours]
        metas += aws_search(awspaths)
        # complete_days
        stime = (time[0]+timedelta(days=1)).replace(hour=0)
        etime = (time[1]-timedelta(days=1)).replace(hour=23)
        if stime < etime:
            n_complete_days = int((etime-stime).days)+1
            complete_days = [(stime + timedelta(days=x)).strftime('%Y/%j/') for x in range(n_complete_days)]
            awspaths = [osp.join(self.base_url.format(self.platform), self.product.format(self.sector), hour) for hour in complete_days]
            metas += aws_search(awspaths)
        # hours last day
        stime = max(time[0],time[1].replace(hour=0))
        etime = time[1]
        n_last_hours = int((etime-stime).total_seconds()/3600)+1
        last_hours = [h for h in ((stime + timedelta(hours=x)).strftime('%Y/%j/%H/') for x in range(n_last_hours)) if h not in first_hours]
        awspaths = [osp.join(self.base_url.format(self.platform), self.product.format(self.sector), hour) for hour in last_hours]
        metas += aws_search(awspaths)
        return metas

    def get_metas_sat(self, bounds, time):
        """
        Get all the meta data for all the necessary products

        :param bounds: bounding box as (lon_min,lon_max,lat_min,lat_max)
        :param time: time interval (init_time_datetime,final_time_datetime)
        :param maxg: max number of granules to process
        :return granules: dictionary with the metadata of all the products
        """
        metas=Dict({})
        metas.geo=self.search_api_sat(self.geo_prefix,bounds,time,collection=self.geo_collection_id)
        metas.fire=self.search_api_sat(self.fire_prefix,bounds,time,collection=self.fire_collection_id)
        metas.geo_nrt=self.search_api_sat(self.geo_nrt_prefix,bounds,time,collection=self.geo_nrt_collection_id)
        metas.fire_nrt=self.search_api_sat(self.fire_nrt_prefix,bounds,time,collection=self.fire_nrt_collection_id)
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

