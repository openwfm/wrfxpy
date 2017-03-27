#
# Dalton Burke, UCDenver
#

from utils import ensure_dir, symlink_unless_exists
from downloader import download_url, DownloadError, get_dList

# fast searching of dList
from bisect import bisect

from datetime import datetime, timedelta

import pytz
import requests
import os
import os.path as osp
import sys
import logging

class level0Error(Exception):
    """
    Raised when a level0 source cannot retreive files
    """
    pass

class level0Source(object):
    """
    Parent class of all level0 sources that implement common functionality, for example
    - local validation (file size check)
    - HDF retrieval with retries (smart check whether server implements http-range)
    - symlinking files
    """

    def __init__(self, ingest_dir):
        """
        Initialize level0 source with ingest directory (where level0 files are stored).

        :param ingest_dir: root of level0 storage
        """
        self.ingest_dir = osp.abspath(ingest_dir)

    def retrieve_level0s(self, from_utc, to_utc, ret_utc = None):
        """
        Attempts to retrieve the level0 files for the time range.
        It should be first verified whether the level0 files are available locally.
        For any unavailable files, downloads should be initiated.

        :param from_utc: start time
        :param to_utc: end time
        :param ref_utc: reference time which defines 'now,' None means datetime.utcnow().
        :return: a list of paths to local level0 files
        """
        pass

    def download_level0(self, url_base, rel_path, max_retries=3):
        """
        Download a level0 file from a level0 service and stream to <rel_path> in ingest_dir.

        :param url_base: the base URL part of the level0 service
        :param rel_path: the relative path of the file
        :param max_retries: how many times we may retry to download the file
        """
        url = url_base + '/' + rel_path
        level0_path = osp.join(self.ingest_dir, rel_path)
        try:
            download_url(url, level0_path, max_retries)
        except DownloadError as e:
            raise level0Error('level0Source: failed to download file %s' % url)

    def level0_available_locally(self, path):
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

    def symlink_level0s(self, manifest, wps_dir):
        """
        Make symlinks in the from LEVEL0FILE.XYZ to all manifest files into wps_dir.

        :param manifest: relative paths (wrt ingest_dir) to level0 files we want linked
        :param wps_dir: the WPS directory where we want the symlinks to appear
        :return:
        """
        for rel_path, level0_name in zip(manifest, generate_level0_names()):
            logging.info('Linking %s -> %s' % (osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, level0_name)))
            symlink_unless_exists(osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, level0_name))

class MODIS_TERRA(level0Source):
    """
    750m data from the MODIS instrument on the Terra satellite
    """

    def __init__(self, ingest_dir):
        super(MODIS_TERRA, self).__init__(ingest_dir)


    def retrieve_level0s(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrive the files to satisfy the simulation request from from_utc to to_utc.

        :param from_utc: start time
        :param to_utc: end time
        :return: list of paths to local level0 files
        """

        # This only works for requests going back about the last two weeks
        # can add a source for older data later, but I don't think it's needed,
        # given the purpose of the project.
        if from_utc < datetime.now(pytz.UTC) - timedelta(days=14):
            raise level0Error('Requested data older than two weeks')

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        manifest = self.compute_manifest(from_utc, to_utc)

        nonlocals = filter(lambda x: not self.level0_available_locally(osp.join(self.ingest_dir, x)), manifest)


        logging.info('Retrieving level0s from %s' % (self.url_base + '/' + self.filepath))

        map(lambda x:self.download_level0(self.url_base + '/' + self.filepath, x), nonlocals)

        return manifest

    def compute_manifest(self, from_utc, to_utc):
        current_time = from_utc

        days = (current_time - datetime(current_time.year, 1, 1, tzinfo=pytz.UTC)).days + 1
        year = current_time.year % 100

        # Give an estimate for the first file that we want, we'll search for it later
        filename = 'P0420064AAAAAAAAAAAAAA%02d%03d%02d%02d%02d001.PDS' % (year, days,
                                                                          current_time.hour,
                                                                          current_time.minute,
                                                                          current_time.second)

        # Filenames have this pattern: P0420064AAAAAAAAAAAAAAyyDDDhhmmss000.PDS
        # Retrieve the directory listing
        dList = get_dList(self.url_base + '/' + self.filepath)

        # Search for where our made up file would go, and start at the pair before it
        index = bisect(dList, filename) - 2

        level0manifest = []

        while current_time < to_utc and index < len(dList):
            level0manifest.append(dList[index])
            level0manifest.append(dList[index+1])

            index = index + 2
            current_file = dList[index]
            # Change time to match the next file, use that time to compare to to_utc
            current_time = current_time.replace(year = 2000 + int(current_file[22:24]))
            current_time = current_time.replace(day=1, month=1)
            current_time = current_time + timedelta(days=int(current_file[24:27]) - 1)
            current_time = current_time.replace(hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest

    url_base = 'ftp://is.sci.gsfc.nasa.gov'
    filepath = 'gsfcdata/terra/modis/level0'

# Near clone of MODIS_TERRA, only changes to url and filename
class MODIS_AQUA(level0Source):
    """
    750m data from the MODIS instrument on the Aqua satellite
    Uniqueness- Requires data from two directories on the source server,
    modis data denoted with _m, and gbad data denoted with _g
    """

    def __init__(self, ingest_dir):
        super(MODIS_AQUA, self).__init__(ingest_dir)
    def retrieve_level0s(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrive the files to satisfy the simulation request from from_utc to to_utc.

        :param from_utc: start time
        :param to_utc: end time
        :return: list of paths to local level0 files
        """

        # This only works for requests going back about the last two weeks
        # can add a source for older data later, but I don't think it's needed,
        # given the purpose of the project.
        if from_utc < datetime.now(pytz.UTC) - timedelta(days=14):
            raise level0Error('Requested data older than two weeks')

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        manifest_m = self.compute_manifest_m(from_utc, to_utc)
        manifest_g = self.compute_manifest_g(from_utc, to_utc)

        nonlocals_m = filter(lambda x: not self.level0_available_locally(osp.join(self.ingest_dir, x)), manifest_m)
        nonlocals_g = filter(lambda x: not self.level0_available_locally(osp.join(self.ingest_dir, x)), manifest_g)

        logging.info('Retrieving level0s from %s' % self.url_base + '/' + self.filepath_m)
        map(lambda x:self.download_level0(self.url_base + '/' + self.filepath_m, x), nonlocals_m)

        logging.info('Retrieving level0s from %s' % self.url_base + '/' + self.filepath_g)
        map(lambda x:self.download_level0(self.url_base + '/' + self.filepath_g, x), nonlocals_g)

        return manifest_m + manifest_g

    def compute_manifest_m(self, from_utc, to_utc):
        current_time = from_utc

        days = (current_time - datetime(current_time.year, 1, 1, tzinfo=pytz.UTC)).days + 1
        year = current_time.year % 100

        # Give an estimate for the first file that we want, we'll search for it later
        filename = 'P1540064AAAAAAAAAAAAAA%02d%03d%02d%02d%02d001.PDS' % (year, days,
                                                                          current_time.hour,
                                                                          current_time.minute,
                                                                          current_time.second)

        # Filenames have this pattern: P1540064AAAAAAAAAAAAAAyyDDDhhmmss000.PDS
        # Retrieve the directory listing
        dList = get_dList(self.url_base + '/' + self.filepath_m)

        # Search for where our made up file would go, and start at the pair before it
        index = bisect(dList, filename) - 2
        current_file = dList[index]

        level0manifest = []

        while current_time < to_utc and index < len(dList):
            level0manifest.append(dList[index])
            level0manifest.append(dList[index+1])

            index = index + 2
            current_file = dList[index]
            # Change time to match the next file, use that time to compare to to_utc
            current_time = current_time.replace(year = 2000 + int(current_file[22:24]))
            current_time = current_time.replace(day=1, month=1)
            current_time = current_time + timedelta(days=int(current_file[24:27]) - 1)
            current_time = current_time.replace(hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest

    def compute_manifest_g(self, from_utc, to_utc):
        current_time = from_utc

        days = (current_time - datetime(current_time.year, 1, 1, tzinfo=pytz.UTC)).days + 1
        year = current_time.year % 100

        # Give an estimate for the first file that we want, we'll search for it later
        filename = 'P1540957AAAAAAAAAAAAAA%02d%03d%02d%02d%02d001.PDS' % (year, days,
                                                                          current_time.hour,
                                                                          current_time.minute,
                                                                          current_time.second)

        # Filenames have this pattern: P1540957AAAAAAAAAAAAAAyyDDDhhmmss000.PDS
        # Retrieve the directory listing
        dList = get_dList(self.url_base + '/' + self.filepath_g)

        # Search for where our made up file would go, and start at the pair before it
        index = bisect(dList, filename) - 4
        current_file = dList[index]

        level0manifest = []

        while current_time < to_utc and index < len(dList):
            level0manifest.append(dList[index])
            level0manifest.append(dList[index+1])

            index = index + 4
            current_file = dList[index]
            # Change time to match the next file, use that time to compare to to_utc
            current_time = current_time.replace(year = 2000 + int(current_file[22:24]))
            current_time = current_time.replace(day=1, month=1)
            current_time = current_time + timedelta(days=int(current_file[24:27]) - 1)
            current_time = current_time.replace(hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest
    url_base = 'ftp://is.sci.gsfc.nasa.gov'
    filepath_m = 'gsfcdata/aqua/modis/level0'
    filepath_g = 'gsfcdata/aqua/gbad'


class VIIRS_NPP(level0Source):
    """
    375m data from VIIRS instrument on the NPP satellite
    """
    def __init__(self, ingest_dir):
        super(VIIRS_NPP, self).__init__(ingest_dir)

    def retrieve_level0s(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrive the files to satisfy the simulation request from from_utc to to_utc.

        :param from_utc: start time
        :param to_utc: end time
        :return: list of paths to local level0 files
        """

        # This only works for requests going back about the last two weeks
        # can add a source for older data later, but I don't think it's needed,
        # given the purpose of the project.
        if from_utc < datetime.now(pytz.UTC) - timedelta(days=14):
            raise level0Error('Requested data older than two weeks')

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        manifest = self.compute_manifest(from_utc, to_utc)

        nonlocals = filter(lambda x: not self.level0_available_locally(osp.join(self.ingest_dir, x)), manifest)


        logging.info('Retrieving level0s from %s' % self.url_base + '/' + self.filepath)

        map(lambda x:self.download_level0(self.url_base + '/' + self.filepath, x), nonlocals)

        return manifest

    def compute_manifest(self, from_utc, to_utc):

        # Give an estimate for the first file that we want, we'll search for it later
        filename = 'RNSCA-RVIRS_npp_d%04d%02d%02d_t%02d00000_e000000_b00000_c00000000000000000000_aaaa_aaa.h5' % (from_utc.year,
                                                                                                            from_utc.month,
                                                                                                            from_utc.day,
                                                                                                            from_utc.hour)
        # format:   RNSCA-RVIRS_npp_dYYYYMMdd_thhmmssS_ehhmmssS_bnnnnn_cnnnnnnnnnnnnnnnnnnnn_aaaa_aaa.h5
        # ex.       RNSCA-RVIRS_npp_d20170220_t0536181_e0547409_b00001_c20170220054740972000_all-_dev.h5

        # Retrieve the directory listing
        dList = get_dList(self.url_base + '/' + self.filepath)

        # Search for where our made up file would go, and start at the file before it, so we can include that data
        index = bisect(dList, filename) - 1
        current_file = dList[index]
        current_time = from_utc

        level0manifest = []

        # there are strange gaps in times between files that I can't reconcile
        # so I just take the start of the next file as current_time
        while current_time < to_utc and index < len(dList):
            level0manifest.append(dList[index])

            index = index + 1
            current_file = dList[index]
            # Change time to match the next file, use that time to compare to to_utc
            current_time = current_time.replace(year=int(current_file[17:21]),
                                                month=int(current_file[21:23]),
                                                day=int(current_file[23:25]),
                                                hour=int(current_file[27:29]),
                                                minute=int(current_file[29:31]),
                                                second=int(current_file[31:33]))

        return level0manifest


    url_base = 'ftp://is.sci.gsfc.nasa.gov'
    filepath = 'gsfcdata/npp/viirs/level0'
