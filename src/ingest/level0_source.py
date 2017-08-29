#
# Dalton Burke, CU Denver
#
# CONUS = [-124.7844079,-66.9513812,24.7433195,49.3457868]

from utils import ensure_dir, symlink_unless_exists
from downloader import download_url, DownloadError, get_dList

# fast searching of dList
from bisect import bisect

from datetime import datetime, timedelta
from pyhdf import SD

import pytz
import requests
import os
import os.path as osp
import sys
import logging

class data_sourceError(Exception):
    """
    Raised when a level0 source cannot retreive files
    """
    pass

class data_source(object):
    """
    Parent class of all data sources that implement common functionality, for example
    - local validation (file size check)
    - HDF retrieval with retries (smart check whether server implements http-range)
    """

    def __init__(self, ingest_dir):
        """
        Initialize level0 source with ingest directory (where level0 files are stored).

        :param ingest_dir: root of level0 storage
        """
        self.ingest_dir = osp.abspath(osp.expanduser(ingest_dir))


    def retrieve_data(self, from_utc, to_utc, lonlat):
        """
        Retrieves all data (geo and active fire) in the given time range and longitude/latitude box. This function is what end users will use to get data

        :param from_utc: start time
        :param to_utc: end time
        :param lonlat: list of form [lowlon, highlon, lowlat, highlat] describing longitude/latitude box
        :return: list of paths to local files that were retrieved
        """
        # I think all data should be ingested into one directory, then whichever files
        # are needed for a given job can be copied to a new folder with a job name

        two_weeks_ago = datetime.utcnow() - timedelta(days=14)
        manifest = []

        if from_utc > two_weeks_ago:
            manifest.extend(self.retrieve_l0(from_utc, to_utc))

        elif to_utc < two_weeks_ago:
            # filter geo_list on intersection with lonlat, the hdf library i'd want to use here isn't ready yet
            geo_list = filter(lambda x: geo_intersects(self.ingest_dir + '/' + x, lonlat), self.retrieve_geo(from_utc, to_utc))
            # geo_list = retrieve_geo(from_utc, to_utc)
            manifest.extend(geo_list)
            manifest.extend(self.retrieve_af(geo_list))

        else:
            manifest.extend(self.retrieve_l0(two_weeks_ago + timedelta(minutes=10), to_utc))

            # filter geo_list on intersection with lonlat
            geo_list = filter(lambda x: geo_intersect(self.ingest_dir + '/' + x, lonlat), self.retrieve_geo(from_utc, two_weeks_ago))
            # geo_list = retrieve_geo(from_utc, two_weeks_ago)
            manifest.extend(geo_list)
            manifest.extend(self.retrieve_af(geo_list))

        return manifest


    def retrieve_geo(self, from_utc, to_utc, ref_utc = None):
        """
        Attempts to retrieve geolocation files in the time range
        First, check if they're available locally, if unavailable proceed to download

        :param from_utc: start time
        :param to_utc: end time
        :return: a list of paths to local geolocation files
        """
        pass


    def compute_geo_manifest(from_utc, to_utc):
        """
        Get list of geolocation file names for the given time frame

        :param from_utc: start time UTC
        :param to_utc: end time UTC
        :return: list of file names as strings
        """
        pass


    def retrieve_af(self, geo_list):
        """
        Attempts to retrieve active fire files in the time range and latitude/longitude box

        :param geo_list: list containing the relevant geolocation file names
        :return: a list of paths to the local active fire files
        """
        pass


    def compute_af_manifest(geo_list):
        """
        get list of active fire file names from a set of geolocation files

        :param geo_list: list containing geolocation file names
        """
        pass


    def retrieve_l0(self, from_utc, to_utc, ref_utc = None):
        """
        Attempts to retrieve the firedata files for the time range.
        It should be first verified whether the firedata files are available locally.
        For any unavailable files, downloads should be initiated.

        :param from_utc: start time
        :param to_utc: end time
        :return: a list of paths to local level0 files
        """
        pass


    def compute_l0_manifest(self, from_utc, to_utc):
        """
        Compute list of files in the source for the given time frame

        :param from_utc: time UTC format
        :param to_utc: time UTC format
        :return: list of file names as strings
        """
        pass


    def download_file(self, url_base, rel_path, max_retries=3):
        """
        Download a file and stream to <rel_path> in ingest_dir.

        :param url_base: the base URL where the file is hosted
        :param rel_path: the relative path of the file
        :param max_retries: how many times we may retry to download the file
        """
        url = url_base + '/' + rel_path
        path = osp.join(self.ingest_dir, rel_path)
        try:
            print 'downloading', url
            download_url(url, path, max_retries)
            print 'done'
        except DownloadError as e:
            raise data_sourceError('data_source: failed to download file %s' % url)


    def available_locally(self, path):
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


class MODIS_TERRA(data_source):
    """
    750m data from the MODIS instrument on the Terra satellite
    """

    def __init__(self, ingest_dir):
        super(MODIS_TERRA, self).__init__(ingest_dir)


    def retrieve_geo(self, from_utc, to_utc):
        """
        Attempts to retrieve geolocation files in the time range
        First, check if they're available locally, if unavailable proceed to download

        :param from_utc: start time
        :param to_utc: end time
        :return: a list of paths to local geolocation files
        """

        manifest = self.compute_geo_manifest(from_utc, to_utc)

        nonlocals = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest)

        logging.info('Retrieving geolocation data from %s' % (self.url_base_hdf + '/' + self.filepath_geo))

        map(lambda x: self.download_file(self.url_base_hdf + '/' + self.filepath_geo + '/' + x[7:11] + '/' + x[11:14], x), nonlocals)

        return manifest


    def compute_geo_manifest(self, from_utc, to_utc):
        """
        Get list of geolocation file names for the given time frame

        :param from_utc: start time UTC
        :param to_utc: end time UTC
        :return: list of file names as strings
        """

        # I don't really want to deal with splitting it on years, so we'll recurse on that
        # from now on we can assume that to and from occur in the same year
        start_year = from_utc.year
        if start_year != to_utc.year:
            return compute_geo_manifest(from_utc, datetime(year=start_year, month=12,day=31,hour=23,minute=59)) + \
                   compute_geo_manifest(datetime(year=start_year+1, month=1, day=1, hour=0, minute=0), to_utc)

        # The source has data for different days in different folders, we'll need to get their paths for each day
        start_day = (from_utc - datetime(start_year, 1,1)).days + 1
        end_day = (to_utc - datetime(start_year, 1, 1)).days + 1

        file_list = []

        for day in range(start_day, end_day + 1):
            file_list.extend(get_dList(self.url_base_hdf + '/' + self.filepath_geo + '/' + str(start_year) + '/' + str(day)))

        # we now have a list with all of the filenames during the days that the query requested, so now we'll trim the stuff at the front and back we don't need
        # invent a sample filename for the start time, they look like this:
        # MOD03.AYYYYDDDD.HHMM.006.#############.hdf
        start_filename = 'MOD03.A%04d%03d.%02d%02d.006.9999999999999.hdf' % (start_year, start_day, from_utc.hour, from_utc.minute)

        # bisect searches for that sample name and returns the index of where that file should go
        # to make sure we get that data we start at the file before it (-1)
        start_index = bisect(file_list, start_filename) - 1

        # we'll do the same for the last one
        end_filename =  'MOD03.A%04d%03d.%02d%02d.006.9999999999999.hdf' % (start_year, end_day, to_utc.hour, to_utc.minute)
        end_index = bisect(file_list, end_filename)

        return file_list[start_index:end_index]


    def retrieve_af(self, geo_list):
        """
        Attempts to retrieve active fire files in the time range and latitude/longitude box

        :param geo_list: list containing the relevant geolocation file names
        :return: a list of paths to the local active fire files
        """
        manifest = self.compute_af_manifest(geo_list)

        nonlocals = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest)

        logging.info('Retrieving active fire data from %s' % (self.url_base_hdf + '/' + self.filepath_af))

        map(lambda x: self.download_file(self.url_base_hdf + '/' + self.filepath_af + '/' + x[7:11] + '/' + x[11:14], x), nonlocals)

        return manifest


    def compute_af_manifest(self, geo_list):
        """
        get list of active fire file names from a set of geolocation files

        :param geo_list: list containing geolocation file names
        """

        prefix = ''
        file_list = []

        for g in geo_list:
            if g[:19] != prefix:
                prefix = g[:19]
                file_list.extend(get_dList(self.url_base_hdf + '/' + self.filepath_af + '/' + str(prefix[7:11]) + '/' + str(prefix[11:14])))

        manifest = []

        # Search for what the name should look like and use that index to add that name to the manifest
        # this takes n*log(n) time, which I think is pretty good
        for g in geo_list:
            manifest.append(file_list[bisect(file_list, 'MOD14' + g[5:24] + '99999999999999.hdf') - 1])

        return manifest


    def retrieve_l0(self, from_utc, to_utc):
        """
        Attempts to retrieve the files to satisfy the simulation request from from_utc to to_utc.

        :param from_utc: start time
        :param to_utc: end time
        :return: list of paths to local level0 files
        """

        # This only works for requests going back about the last two weeks
        # can add a source for older data later, but I don't think it's needed,
        # given the purpose of the project.
        if from_utc < datetime.utcnow() - timedelta(days=14):
            raise data_sourceError('Requested data older than two weeks')

        manifest = self.compute_l0_manifest(from_utc, to_utc)

        nonlocals = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest)


        logging.info('Retrieving level0s from %s' % (self.url_base_l0 + '/' + self.filepath_l0))

        map(lambda x:self.download_file(self.url_base_l0 + '/' + self.filepath_l0, x), nonlocals)

        return manifest

    def compute_l0_manifest(self, from_utc, to_utc):
        """
        Compute list of files in the source for the given time frame

        :param from_utc: time UTC format
        :param to_utc: time UTC format
        :return: list of file names as strings
        """

        # Retrieve the directory listing
        dList = get_dList(self.url_base_l0 + '/' + self.filepath_l0)

        # We want a list of all of the filenames which land between from_utc and to_utc

        # Gameplan:
        # What would a file that starts exactly at from_utc look like?
        # Filenames have this pattern: P0420064AAAAAAAAAAAAAAyyDDDhhmmss000.PDS
        current_time = from_utc

        days = (current_time - datetime(current_time.year, 1, 1)).days + 1
        year = current_time.year % 100

        filename = 'P0420064AAAAAAAAAAAAAA%02d%03d%02d%02d%02d000.PDS' % (year, days,
                                                                          current_time.hour,
                                                                          current_time.minute,
                                                                          current_time.second)

        # Then, we find out where that filename would go in the dList
        # This call binary searches dList for filename, and returns it's index (pretty efficient)
        # If the filename is not found, it returns the index of the first file larger than it
        index = bisect(dList, filename)

        # If the filename we made up is not in the list (very likely), we actually want the first file
        # smaller than the filename, so we still get the data for that time period
        # (-2 since the files come in pairs, one that ends in 000.PDS and one that ends in 001.PDS)

        if index == len(dList):
            index = index - 2

        elif dList[index] != filename:
            index = index - 2

        level0manifest = []

        # Now that we know where to start, we'll begin filling the manifest with relevant files
        while current_time < to_utc:
            # Add 000.PDS file to manifest
            level0manifest.append(dList[index])
            # Add 001.PDS file to manifest
            level0manifest.append(dList[index+1])

            # Move the index to the next pair, if we run out of files just break
            index = index + 2
            if index >= len(dList):
                break

            current_file = dList[index]

            # Change time to match the next file, use that time to compare to to_utc
            # If the time that we get from this exceeds to_utc, we have all the data we want
            current_time = current_time.replace(year = 2000 + int(current_file[22:24]))
            current_time = current_time.replace(day=1, month=1)
            current_time = current_time + timedelta(days=int(current_file[24:27]) - 1)
            current_time = current_time.replace(hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest

    url_base_l0 = 'ftp://is.sci.gsfc.nasa.gov'
    filepath_l0 = 'gsfcdata/terra/modis/level0'

    url_base_hdf = 'ftp://ladsweb.nascom.nasa.gov'
    filepath_geo = 'allData/6/MOD03'
    filepath_af = 'allData/6/MOD14'

# Near clone of MODIS_TERRA, only changes to url and filename
class MODIS_AQUA(data_source):
    """
    750m data from the MODIS instrument on the Aqua satellite
    Uniqueness- Requires data from two directories on the source server,
    modis data denoted with _m, and gbad data denoted with _g
    """

    def __init__(self, ingest_dir):
        super(MODIS_AQUA, self).__init__(ingest_dir)

    def retrieve_geo(self, from_utc, to_utc):
        """
        Attempts to retrieve geolocation files in the time range
        First, check if they're available locally, if unavailable proceed to download

        :param from_utc: start time
        :param to_utc: end time
        :return: a list of paths to local geolocation files
        """
        manifest = self.compute_geo_manifest(from_utc, to_utc)

        nonlocals = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest)

        logging.info('Retrieving geolocation data from %s' % (self.url_base_hdf + '/' + self.filepath_geo))

        map(lambda x: self.download_file(self.url_base_hdf + '/' + self.filepath_geo + '/' + x[7:11] + '/' + x[11:14], x), nonlocals)

        return manifest


    def compute_geo_manifest(self, from_utc, to_utc):
        """
        Get list of geolocation file names for the given time frame

        :param from_utc: start time UTC
        :param to_utc: end time UTC
        :return: list of file names as strings
        """


        # I don't really want to deal with splitting it on years, so we'll recurse on that
        # from now on we can assume that to and from occur in the same year
        start_year = from_utc.year
        if start_year != to_utc.year:
            return compute_geo_manifest(from_utc, datetime(year=start_year, month=12,day=31,hour=23,minute=59)) + \
                   compute_geo_manifest(datetime(year=start_year+1, month=1, day=1, hour=0, minute=0), to_utc)

        # The source has data for different days in different folders, we'll need to get their paths for each day
        start_day = (from_utc - datetime(start_year, 1,1)).days + 1
        end_day = (to_utc - datetime(start_year, 1, 1)).days + 1

        file_list = []

        for day in range(start_day, end_day + 1):
            file_list.extend(get_dList(self.url_base_hdf + '/' + self.filepath_geo + '/' + str(start_year) + '/' + str(day)))


        # we now have a list with all of the filenames during the days that the query requested, so now we'll trim the stuff at the front and back we don't need
        # invent a sample filename for the start time, they look like this:
        # MYD03.AYYYYDDDD.HHMM.006.#############.hdf
        start_filename = 'MYD03.A%04d%03d.%02d%02d.006.9999999999999.hdf' % (start_year, start_day, from_utc.hour, from_utc.minute)

        # bisect searches for that sample name and returns the index of where that file should go
        # to make sure we get that data we start at the file before it (-1)
        start_index = bisect(file_list, start_filename) - 1

        # we'll do the same for the last one
        end_filename =  'MYD03.A%04d%03d.%02d%02d.006.9999999999999.hdf' % (start_year, end_day, to_utc.hour, to_utc.minute)
        end_index = bisect(file_list, end_filename)


        return file_list[start_index:end_index]


    def retrieve_af(self, geo_list):
        """
        Attempts to retrieve active fire files in the time range and latitude/longitude box

        :param geo_list: list containing the relevant geolocation file names
        :return: a list of paths to the local active fire files
        """
        manifest = self.compute_af_manifest(geo_list)

        nonlocals = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest)

        logging.info('Retrieving active fire data from %s' % (self.url_base_hdf + '/' + self.filepath_af))

        map(lambda x: self.download_file(self.url_base_hdf + '/' + self.filepath_af + '/' + x[7:11] + '/' + x[11:14], x), nonlocals)

        return manifest




    def compute_af_manifest(self, geo_list):
        """
        get list of active fire file names from a set of geolocation files

        :param geo_list: list containing geolocation file names
        """

        prefix = ''
        file_list = []

        for g in geo_list:
            if g[:19] != prefix:
                prefix = g[:19]
                file_list.extend(get_dList(self.url_base_hdf + '/' + self.filepath_af + '/' + str(prefix[7:11]) + '/' + str(prefix[11:14])))

        manifest = []

        # Search for what the name should look like and use that index to add that name to the manifest
        # this takes n*log(n) time, which I think is pretty good
        for g in geo_list:
            manifest.append(file_list[bisect(file_list, 'MYD14' + g[5:24] + '99999999999999.hdf') - 1])

        return manifest



    def retrieve_l0(self, from_utc, to_utc):
        """
        Attempts to retrive the files to satisfy the simulation request from from_utc to to_utc.

        :param from_utc: start time
        :param to_utc: end time
        :return: list of paths to local level0 files
        """

        # This only works for requests going back about the last two weeks
        # can add a source for older data later, but I don't think it's needed,
        # given the purpose of the project.
        if from_utc < datetime.utcnow() - timedelta(days=14):
            raise data_sourceError('Requested data older than two weeks')

        manifest_m = self.compute_l0_manifest_m(from_utc, to_utc)
        manifest_g = self.compute_l0_manifest_g(from_utc, to_utc)

        nonlocals_m = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest_m)
        nonlocals_g = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest_g)

        logging.info('Retrieving level0s from %s' % self.url_base_l0 + '/' + self.filepath_l0_m)
        map(lambda x:self.download_file(self.url_base_l0 + '/' + self.filepath_l0_m, x), nonlocals_m)

        logging.info('Retrieving level0s from %s' % self.url_base_l0 + '/' + self.filepath_l0_g)
        map(lambda x:self.download_file(self.url_base_l0 + '/' + self.filepath_l0_g, x), nonlocals_g)

        return manifest_m + manifest_g

    def compute_l0_manifest_m(self, from_utc, to_utc):
        """
        Compute list of MODIS files in the source for the given time frame

        :param from_utc: time UTC format
        :param to_utc: time UTC format
        :return: list of file names as strings
        """

        # We want a list of all of the filenames which land between from_utc and to_utc

        # Retrieve the directory listing
        dList = get_dList(self.url_base_l0 + '/' + self.filepath_l0_m)

        # Gameplan:
        # What would a file that starts exactly at from_utc look like?
        # Filenames have this pattern: P1540064AAAAAAAAAAAAAAyyDDDhhmmss000.PDS
        current_time = from_utc

        days = (current_time - datetime(current_time.year, 1, 1)).days + 1
        year = current_time.year % 100

        filename = 'P1540064AAAAAAAAAAAAAA%02d%03d%02d%02d%02d000.PDS' % (year, days,
                                                                          current_time.hour,
                                                                          current_time.minute,
                                                                          current_time.second)

        # Then, we find out where that filename would go in the dList
        # This call binary searches dList for filename, and returns it's index (pretty efficient)
        # If the filename is not found, it returns the index of the first file larger than it
        index = bisect(dList, filename)

        # If the filename we made up is not in the list (very likely), we actually want the first file
        # smaller than the filename, so we still get the data for that time period
        # (-2 since the files come in pairs, one that ends in 000.PDS and one that ends in 001.PDS)


        if index == len(dList):
            index = index - 2
        elif dList[index] != filename:
            index = index - 2

        level0manifest = []

        while current_time < to_utc:
            # Add 000.PDS File
            level0manifest.append(dList[index])
            # Add 001.PDS file
            level0manifest.append(dList[index+1])

            # Move index to start of next pair,

            index = index + 2
            if index >= len(dList):
                break

            current_file = dList[index]

            # Change time to match the next file, use that time to compare to to_utc
            # If the time on the next file is bigger than to_utc, then we have all the files we care about
            current_time = current_time.replace(year = 2000 + int(current_file[22:24]))
            current_time = current_time.replace(day=1, month=1)
            current_time = current_time + timedelta(days=int(current_file[24:27]) - 1)
            current_time = current_time.replace(hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest

    def compute_l0_manifest_g(self, from_utc, to_utc):
        """
        Compute list of GBAD files (AQUA specific) in the source for the given time frame

        :param from_utc: time UTC format
        :param to_utc: time UTC format
        :return: list of file names as strings
        """

        # We want a list of all of the filenames which land between from_utc and to_utc

        # Retrieve the directory listing
        dList = get_dList(self.url_base_l0 + '/' + self.filepath_l0_g)

        # Gameplan:
        # What would a file that starts exactly at from_utc look like?
        # Filenames have this pattern: P1540064AAAAAAAAAAAAAAyyDDDhhmmss000.PDS
        current_time = from_utc

        days = (current_time - datetime(current_time.year, 1, 1)).days + 1
        year = current_time.year % 100

        filename = 'P1540957AAAAAAAAAAAAAA%02d%03d%02d%02d%02d000.PDS' % (year, days,
                                                                          current_time.hour,
                                                                          current_time.minute,
                                                                          current_time.second)

        # Then, we find out where that filename would go in the dList
        # This call binary searches dList for filename, and returns it's index (pretty efficient)
        # If the filename is not found, it returns the index of the first file larger than it
        index = bisect(dList, filename)

        # If the filename we made up is not in the list (very likely), we actually want the first file
        # smaller than the filename, so we still get the data for that time period
        # (-4 because for each time there are 4 GBAD files, however there are only 2 we care for)
        if index == len(dList):
            index = index - 4
        elif dList[index] != filename:
            index = index - 4

        level0manifest = []

        while current_time < to_utc:
            # Add 000.PDS file
            level0manifest.append(dList[index])
            # Add 001.PDS file
            level0manifest.append(dList[index+1])

            # Move index to next pair, (remember, there are 4 GBAD files, we only care about 2 of them)
            # If we run out of filenames before reaching to_utc, that's fine, just break
            index = index + 4
            if index >= len(dList):
                break

            current_file = dList[index]

            # Change time to match the next file, use that time to compare to to_utc
            # If the new time is bigger than to_utc, we have all of the files we care about
            current_time = current_time.replace(year = 2000 + int(current_file[22:24]))
            current_time = current_time.replace(day=1, month=1)
            current_time = current_time + timedelta(days=int(current_file[24:27]) - 1)
            current_time = current_time.replace(hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest


    url_base_l0 = 'ftp://is.sci.gsfc.nasa.gov'
    filepath_l0_m = 'gsfcdata/aqua/modis/level0'
    filepath_l0_g = 'gsfcdata/aqua/gbad'

    url_base_hdf = 'ftp://ladsweb.nascom.nasa.gov'
    filepath_geo = 'allData/6/MYD03'
    filepath_af = 'allData/6/MYD14'



class VIIRS_NPP(data_source):
    """
    375m data from VIIRS instrument on the NPP satellite
    """
    def __init__(self, ingest_dir):
        super(VIIRS_NPP, self).__init__(ingest_dir)

    def retrieve_l0(self, from_utc, to_utc):
        """
        Attempts to retrive the files to satisfy the simulation request from from_utc to to_utc.

        :param from_utc: start time
        :param to_utc: end time
        :return: list of paths to local level0 files
        """

        # This only works for requests going back about the last two weeks
        # can add a source for older data later, but I don't think it's needed,
        # given the purpose of the project.
        if from_utc < datetime.utcnow() - timedelta(days=14):
            raise data_sourceError('Requested data older than two weeks')

        manifest = self.compute_l0_manifest(from_utc, to_utc)

        nonlocals = filter(lambda x: not self.available_locally(osp.join(self.ingest_dir, x)), manifest)


        logging.info('Retrieving level0s from %s' % self.url_base_l0 + '/' + self.filepath_l0)

        map(lambda x:self.download_file(self.url_base_l0 + '/' + self.filepath_l0, x), nonlocals)

        return manifest

    def compute_l0_manifest(self, from_utc, to_utc):
        """
        Compute list of files in the source for the given time frame

        :param from_utc: time UTC format
        :param to_utc: time UTC format
        :return: list of file names as strings
        """
        # We want a list of all of the filenames which land between from_utc and to_utc

        # Retrieve the directory listing
        dList = get_dList(self.url_base_l0 + '/' + self.filepath_l0)

        # Gameplan:
        # What would a file that starts exactly at from_utc look like?
        # format:   RNSCA-RVIRS_npp_dYYYYMMdd_thhmmssS_ehhmmssS_bnnnnn_cnnnnnnnnnnnnnnnnnnnn_aaaa_aaa.h5
        filename = 'RNSCA-RVIRS_npp_d%04d%02d%02d_t%02d00000_e000000_b00000_c00000000000000000000_aaaa_aaa.h5' % (from_utc.year,
                                                                                                            from_utc.month,
                                                                                                            from_utc.day,
                                                                                                            from_utc.hour)

        # Then, we find out where that filename would go in the dList
        # This call binary searches dList for filename, and returns it's index (pretty efficient)
        # If the filename is not found, it returns the index of the first file larger than it
        index = bisect(dList, filename)

        # If the filename we made up is not in the list (very likely), we actually want the first file
        # smaller than the filename, so we still get the data for that time period
        if index == len(dList):
            index = index - 1
        elif dList[index] != filename:
            index = index - 1

        current_time = from_utc

        level0manifest = []

        # there are strange gaps in times between files that I can't reconcile
        # so I just take the start of the next file as current_time
        while current_time < to_utc:
            # Get the file
            level0manifest.append(dList[index])

            index = index + 1
            if index >= len(dList):
                break

            current_file = dList[index]
            # Change time to match the next file, use that time to compare to to_utc
            # If the time of the next file is bigger than to_utc, then we have all of the files we care about
            current_time = current_time.replace(year=int(current_file[17:21]),
                                                month=int(current_file[21:23]),
                                                day=int(current_file[23:25]),
                                                hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest


    url_base_l0 = 'ftp://is.sci.gsfc.nasa.gov'
    filepath_l0 = 'gsfcdata/npp/viirs/level0'




def geo_intersects(filename, lonlat):
    """
    Checks a geolocation file for overlap with a latitude longitude box
    :filename: name of file to check
    :lonlat: list, [leftlon, rightlon, botlat, toplat]
    :return: boolean, true if there was overlap
    """
    print "Checking",filename, "..."

    if filename[-4:] != '.hdf':
        print "ERROR: File", filename, "is not the correct filetype (require hdf)"
        return False

    if lonlat[0] > lonlat[1]:
        print "ERROR: Requested box crosses dateline, no support for this yet"
        return False

    try:
        hdf = SD.SD(filename)
    except:
        print "ERROR: Could not load file " + filename + '\n'
        return False

    lon = hdf.select('Longitude')
    lat = hdf.select('Latitude')

    dim1 = len(lon[:])
    dim2 = len(lon[0])

    minlon = float(lon[0][0])
    maxlon = float(lon[dim1 - 1][dim2 - 1])

    minlat = float(lat[dim1 - 1][dim2 - 1])
    maxlat = float(lat[0][0])

    if minlon > maxlon:
        #print "ERROR: File " + filename + " crosses dateline, no support for this yet\n"
        return False


    lonoverlap = minlon < lonlat[1] and maxlon > lonlat[0]
    latoverlap = minlat < lonlat[3] and maxlat > lonlat[2]

    return lonoverlap and latoverlap
