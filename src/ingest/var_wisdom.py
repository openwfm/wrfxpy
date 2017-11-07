#
# Dalton Burke, CU Denver
#

from utils import ensure_dir, symlink_unless_exists
from downloader import download_url, DownloadError, get_dList

from datetime import datetime, timedelta
import sys
import logging
import os.path as osp
from bisect import bisect



_gran_wisdom = {
    # MODIS_AQUA
      'MYD03': {
            'name' : 'MYD03',
            'type' : 'geolocation',
            'satellite' : 'MODIS_AQUA',
            'source' : 'laadsweb',
            'rel_path' : 'allData/6/',
            'lat_path' : 'Latitude',
            'lon_path' : 'Longitude',
      },
      'MYD021KM' : {
            'name' : 'MYD021KM',
            'type' : 'radiance',
            'satellite' : 'MODIS_AQUA',
            'source' : 'laadsweb',
            'rel_path' : 'allData/6/',
            'lat_path' : 'Latitude',
            'lon_path' : 'Longitude',
      },
      #'MYD02HKM' : {
      #      'name' : 'MYD02HKM',
      #      'type' : 'radiance'
      #      'satellite' : 'MODIS_AQUA',
      #      'source' : 'laadsweb',
      #      'lat_path' : 'Latitude',
      #      'lon_path' : 'Longitude',
      #},
      #'MYD02QKM' : {
      #      'name' : 'MYD02QKM',
      #      'type' : 'radiance'
      #      'satellite' : 'MODIS_AQUA',
      #      'source' : 'laadsweb',
      #      'rel_path' : 'allData/6/',
      #      'lat_path' : 'Latitude',
      #      'lon_path' : 'Longitude',
      #},
      'MYD14': {
            'name' : 'MYD14',
            'type' : 'active fire',
            'satellite' : 'MODIS_AQUA',
            'source' : 'laadsweb',
            'rel_path' : 'allData/6/'

      },

    # MODIS_TERRA
      'MOD03' : {
            'name' : 'MOD03',
            'type' : 'geolocation',
            'satellite' : 'MODIS_TERRA',
            'source' : 'laadsweb',
            'rel_path' : 'allData/6/',
            'lat_path' : 'Latitude',
            'lon_path' : 'Longitude',
      },
      'MOD021KM' : {
            'name' : 'MOD021KM',
            'type' : 'radiance',
            'satellite' : 'MODIS_TERRA',
            'source' : 'laadsweb',
            'rel_path' : 'allData/6/',
            'lat_path' : 'Latitude',
            'lon_path' : 'Longitude',
      },
      #'MOD02HKM' : {
      #      'name' : 'MOD02HKM',
      #      'type' : 'radiance',
      #      'source' : 'laadsweb'
      #      'satellite' : 'MODIS_TERRA',
      #      'lat_path' : 'Latitude',
      #      'lon_path' : 'Longitude',
      #},
      #'MOD02QKM' : {
      #      'name' : 'MOD02QKM',
      #      'type' : 'radiance',
      #      'source' : 'laadsweb',
      #      'satellite' : 'MODIS_TERRA',
      #      'rel_path' : 'allData/6/',
      #      'lat_path' : 'Latitude',
      #      'lon_path' : 'Longitude',
      #},
      'MOD14' : {
            'name' : 'MOD14',
            'type' : 'active fire',
            'satellite' : 'MODIS_TERRA',
            'source' : 'laadsweb',
            'rel_path' : 'allData/6/',

      },

    # VIIRS_NPP
      'NPP_IMFTS_L1' : {
            'name' : 'NPP_IMFTS_L1',
            'type' : 'geolocation',
            'satellite' : 'VIIRS_NPP',
            'source' : 'laadsweb',
            'rel_path' : 'allData/5000/',
      },

}
def get_gran_wisdom(var_name):
    """Return wisdom for the variable <var_name>."""
    return _gran_wisdom[var_name]

def get_gran_wisdom_variables():
    """Return the variables for which wisdom is available."""
    return _gran_wisdom.keys()


def laads_retr_manifest(ingest_dir, manifest):
      ingest_dir = osp.abspath(osp.expanduser(ingest_dir))

      nonlocals = filter(lambda x: not available_locally(osp.join(ingest_dir, x)), manifest)

      logging.info('Retrieving data from laadsweb...')

      map(lambda x: download_file(ingest_dir, laads_file_to_url(x), x), nonlocals)

      return manifest


def laads_geo_from_geoMeta(ingest_dir, satellite, from_utc, to_utc, lonlat):

      geoMeta_manifest = retrieve_geoMeta(ingest_dir, satellite, from_utc, to_utc)

      geoMeta = []
      for fname in geoMeta_manifest:
            with open(fname) as f:
                  # first line contains description of file structure, we don't want it
                  next(f)
                  # split each line into a list, items seperated by commas
                  for line in f:
                        geoMeta.append(line.split(','))

      # each line in the file has southbounding coord in pos 5, north bound in pos 8, west bound in 7, east bound in 6.
      # using these to determine if the particular file intersects the lonlat range
      # we select the file name (in position 0) using list comprehension method
      geo_manifest = [item[0] for item in filter(\
                                                 lambda x:\
                                                 lonlat_intersect(lonlat, [float(x[5]), float(x[8]), float(x[7]), float(x[6])])\
                                                 and (x[1] > str(from_utc) and x[1] < str(to_utc)), geoMeta)]

      return geo_manifest

def laads_range_manifest(gran, from_utc, to_utc):
      """
      creates manifest for files of a granule from given time range
      """
      start_year = from_utc.year
      if start_year != to_utc.year:
          return laads_range_manifest(gran, from_utc, datetime(year=start_year, month=12,day=31,hour=23,minute=59)) + \
                 laads_range_manifest(gran, datetime(year=start_year+1, month=1, day=1, hour=0, minute=0), to_utc)

      # The source has data for different days in different folders, we'll need to get their paths for each day
      start_day = (from_utc - datetime(start_year, 1,1)).days + 1
      end_day = (to_utc - datetime(start_year, 1, 1)).days + 1

      file_list = []
      url = 'ftp://ladsweb.nascom.nasa.gov/' + gran['rel_path'] + gran['name'] + '/%s/%s'

      for day in range(start_day, end_day + 1):
          file_list.extend(get_dList(url % (str(start_year), str(day))))

      # we now have a list with all of the filenames during the days that the query requested, so now we'll trim the stuff at the front and back we don't need
      # invent a sample filename for the start time, they look like this:
      # MOD03.AYYYYDDDD.HHMM.006.#############.hdf
      start_filename = '%s.A%04d%03d.%02d%02d.006.9999999999999.hdf' % (gran['name'], start_year, start_day, from_utc.hour, from_utc.minute)

      # bisect searches for that sample name and returns the index of where that file should go
      # to make sure we get that data we start at the file before it (-1)
      start_index = bisect(file_list, start_filename) - 1

      # we'll do the same for the last one
      end_filename =  '%s.A%04d%03d.%02d%02d.006.9999999999999.hdf' % (gran['name'], start_year, end_day, to_utc.hour, to_utc.minute)
      end_index = bisect(file_list, end_filename)

      return file_list[start_index:end_index]


def laads_list_manifest(gran, gran_list):
      """
      creates manifest for a given granule matching time signatures of given granules.
      """

      url = 'ftp://ladsweb.nascom.nasa.gov/' + gran['rel_path'] + gran['name'] + '/%s/%s'
      prefix = ''
      file_list = []

      i = gran_list[0].find('.')

      for g in gran_list:
            if g[:i+10] != prefix:
                prefix = g[:i+10]
                file_list.extend(get_dList(url % (prefix[i+2:i+6], prefix[i+6:i+9])))

      search_string = gran['name'] + '%s.9999999999999.hdf'
      manifest = map(lambda x: file_list[bisect(file_list, search_string % x[i:(i+18)])-1], gran_list)
      return manifest

def laads_file_to_url(filename):
      """
      converts a filename into it's respective url
      :filename: given filename
      :return: url for the file
      """
      attr = filename.split('.')
      g_name = attr[0]
      rel_path = get_gran_wisdom(g_name)['rel_path']
      year = attr[1][1:5]
      day = attr[1][5:9]
      return 'ftp://ladsweb.nascom.nasa.gov/%s%s/%s/%s/%s' % (rel_path, g_name, year, day, filename)


_source_wisdom = {
      'laadsweb' : {
            'url' : 'ftp://ladsweb.nascom.nasa.gov',
            'retr_manifest' : lambda x,y: laads_retr_manifest(x,y),
            'list_manifest' : lambda x,y: laads_list_manifest(x,y),
            'range_manifest' : lambda x,y,z: laads_range_manifest(x,y,z),
            'geo_from_geoMeta' : lambda v,w,x,y,z: laads_geo_from_geoMeta(v,w,x,y,z),
            'file_to_url' : lambda x: laads_file_to_url(x),
      },
}

def get_source_wisdom(source_name):
      return _source_wisdom[source_name]

def get_source_wisdom_variables():
      return _source_wisdom.keys()

_satellite_wisdom = {
      # Geolocation granule should come first here
      'MODIS_AQUA'  : ['MYD03','MYD021KM','MYD14'],
      'MODIS_TERRA' : ['MOD03','MYD021KM','MYD14'],
      # 'VIIRS_NPP'   : ['NPP_IMFTS_L1'],
}

def get_sat_wisdom(satellite_name):
      return _satellite_wisdom[satellite_name]

def get_sat_wisdom_variables():
      return _satellite_wisdom.keys()


def retr_sat(ingest_dir, satellite, from_utc, to_utc, lonlat = []):
      """
      Downloads all implemented granules of a satellite in a given time and lonlat (optional) range
      """
      ingest_dir = osp.abspath(osp.expanduser(ingest_dir))
      sat_gran_names = get_sat_wisdom(satellite)
      sat_gran_wis = map(lambda x: get_gran_wisdom(x), sat_gran_names)

      if lonlat != []:
            for gran in sat_gran_wis:
                if gran['type'] == 'geolocation':
                    geo_gran = gran

            geo_source_wis = get_source_wisdom(geo_gran['source'])

            if geo_gran['name'] == 'MOD03':
                  geo_manifest = geo_source_wis['geo_from_geoMeta'](ingest_dir, 'TERRA', from_utc, to_utc, lonlat)
            elif geo_gran['name'] == 'MYD03':
                  geo_manifest = geo_source_wis['geo_from_geoMeta'](ingest_dir, 'AQUA', from_utc, to_utc, lonlat)
            else:
                  geo_manifest = geo_source_wis['range_manifest'](geo_gran, from_utc, to_utc)
      else:
            geo_manifest = geo_source_wis['range_manifest'](geo_gran, from_utc, to_utc)

      geo_source_wis['retr_manifest'](ingest_dir, geo_manifest)

      # need copy here (shallow is fine), use list()
      # this caused a lot of pain
      manifest = list(geo_manifest)

      for gran in sat_gran_wis:
          if gran['name'] != geo_gran['name']:
                  source_wis = get_source_wisdom(gran['source'])
                  logging.info('Downloading %s from source %s' % (gran['name'], gran['source']))
                  gran_manifest = source_wis['list_manifest'](gran, geo_manifest)
                  source_wis['retr_manifest'](ingest_dir, gran_manifest)
                  manifest.extend(gran_manifest)

      return manifest


def available_locally(path):
        """
        Check in a level0 file is available locally and if it's file size checks out.
        :param path: the level0 file path
        """
        info_path = path + '.size'
        if osp.exists(path) and osp.exists(info_path):
            content_size = int(open(info_path).read())
            return osp.getsize(path) == content_size
        else:
            return False


def download_file(ingest_dir, url, rel_path, max_retries=3):
    """
    Download a file and stream to <rel_path> in ingest_dir.

    :param url_base: the base URL where the file is hosted
    :param rel_path: the relative path of the file
    :param max_retries: how many times we may retry to download the file
    """
    # logging.info("Downloading %s from %s" % (rel_path, url))
    path = osp.join(ingest_dir, rel_path)
    try:
        download_url(url, path, max_retries)
    except DownloadError as e:
        raise data_sourceError('data_source: failed to download file %s' % url)


def geoMeta_manifest(from_utc, to_utc, satellite):
      """
      retrieves geoMeta files for a given granule in a given time frame
      geoMeta files give the lonlat range of all geolocation files for the given day
      """
      logging.info('Computing manifest of geoMeta files needed...')

      if satellite == 'TERRA':
            gran_name = 'MOD03'
      elif satellite == 'AQUA':
            gran_name = 'MYD03'
      else:
            raise Exception(ValueError)

      i = from_utc.replace(hour=0,minute=0,second=0,microsecond=0)
      end_date = to_utc.replace(hour=0,minute=0,second=0,microsecond=0)

      geoMeta_manifest = []
      url = 'ftp://ladsweb.nascom.nasa.gov/geoMeta/6/' + satellite + '/%04d/' + gran_name + '_%04d-%02d-%02d.txt'

      while i <= end_date:
          geoMeta_manifest.append(url % (i.year, i.year, i.month, i.day))
          i = i + timedelta(days=1)

      return geoMeta_manifest


def retrieve_geoMeta(ingest_dir, satellite, from_utc, to_utc):
      """
      attempts to retrieve geoMeta files in the time range
      """

      logging.info("Retrieving geoMeta files...")

      if satellite != 'TERRA' and satellite != 'AQUA':
            raise Exception(ValueError)

      manifest = geoMeta_manifest(from_utc, to_utc, satellite)

      nonlocals = filter(lambda x: not available_locally(osp.join(ingest_dir, 'geoMeta/' + x.split('/')[-1])), manifest)

      logging.info('Retrieving geolocation data from %s' % ('ftp://ladsweb.nascom.nasa.gov/geoMeta/6/' + satellite))

      map(lambda x: download_file(ingest_dir, x, osp.join(ingest_dir, 'geoMeta/' + x.split('/')[-1])), nonlocals)

      return [(osp.join(ingest_dir, 'geoMeta/' + x.split('/')[-1])) for x in manifest]



def lonlat_intersect(lonlat1, lonlat2):
      """
      Check if two lonlat boxes intersect
      :lonlat1: list, [leftlon, rightlon, botlat, toplat]
      :lonlat2: list, [leftlon, rightlon, botlat, toplat]
      :return: true if lonlat boxes overlap, false otherwise
      """
      # Before continuing, note that a common way to check if two intervals [a,b] and [x,y] overlap
      # is to check b > x AND a < y. In this function we're checking if two BOXes overlap
      # which only happens if there's a vertical overlap AND a horizontal overlap.

      # if the latitude intervals don't overlap we have nothing further to discuss
      if not (lonlat1[2] < lonlat2[3] and lonlat1[3] > lonlat2[2]):
            return False

      # Now all we have to do is determine if the longitude intervals intersect

      # if both cross dateline, then the longitude intervals intersect at the dateline
      if lonlat1[0] > lonlat1[1] and lonlat2[0] > lonlat2[1]:
            return True

      # if only lonlat1 crosses dateline, check intervals from lonlat1[0] to 180 and -180 to lonlat1[1].
      # the computation simplifies because 180 will always be greater than lonlat2[1] and -180 always smaller than lonlat2[0]
      elif lonlat1[0] > lonlat1[1]:
            return lonlat1[0] < lonlat2[1] or lonlat1[1] > lonlat2[0]

      elif lonlat2[0] > lonlat2[1]:
            return lonlat2[0] < lonlat1[1] or lonlat2[1] > lonlat1[0]

      # neither cross dateline, check if longitudes overlap like normal
      else:
            return lonlat1[0] < lonlat2[1] and lonlat1[1] > lonlat2[0]
