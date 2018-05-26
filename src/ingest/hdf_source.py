#
# Dalton Burke, CU Denver
#


from utils import ensure_dir, symlink_unless_exists
from downloader import download_url, DownloadError

from datetime import datetime, timedelta
import pytz
import requests
import os
import os.path as osp
import sys
import logging

class HDFError(Exception):
    """
    Raised when an HDFSource cannot retrieve HDFs.
    """
    pass

class HDFSource(object):
    """
    the parent class of all HDF sources that implement common functionality, for example

    - local HDF validation (file size check)
    - HDF retrieval with retries (smart check whether server implements http-range)
    - symlinking HDF files for unhdf
    """

    def __init(self, ingest_dir):
        """
        Initialize HDF source with ingest directory (where HDF files are stored).

        :param ingest_dir: root of HDF storage
        """
        self.ingest_dir = osp.abspath(ingest_dir)

    def vtables(self):
        """
        Returns the vtables that must be used with this source as a table with keys:
        __________________________________

        :return: a dictionary mapping vtable keys to specific table files
        """

        return {}

    def namelist_keys(self):
        """
        Some grib2 source filse require that in namelist.input, certain parameters have particular values, such keys should be returned here

        :return: a dictionary mapping section names that must be modified
        """
        return {}

    def retrieve_hdfs(self, from_utc, to_utc, ret_utc = None):
        """
        Attempts to retrieve the HDF files for the time range.
        It should be first verified whether the HDF files are available locally.
        For any unavailable files, downloads should be initiated.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :param ref_utc: reference time which defines 'now,' None means datetime.utcnow().
        :return: a list of paths to local HDF files
        """
        pass

    def download_hdf(self, url_base, rel_path, max_retries=3):
        """
        Download an HDF file from an HDF service and stream to <rel_path> in ingest_dir.

        :param url_base: the base URL part of the HDF service
        :param rel_path: the relative path of the file
        :param max_retries: how many times we may retry to download the file
        """
        url = url_base + '/' + rel_path
        hdf_path = osp.join(self.ingest_dir, rel_path)
        try:
            download_url(url, hdf_path, max_retries)
        except DownloadError as e:
            raise HdfError('HDFSource: failed to download file %s' % url)

    def hdf_available_locally(self, path):
        """
        Check if an HDF file is available locally and if it's file size checks out.

        :param path: the HDF file path
        """
        info_path = path + '.size'
        if osp.exists(path) and osp.exists(info_path):
            content_size = int(open(info_path).read())
            return osp.getsize(path) = content_size
        else:
            return False

        def symlink_gribs(self, manifest, wps_dir):
            """
            Make symlinks in the form HDFFILE.XYZ to all manifest files into wps_dir.

            :param manifest: relative paths (wrt ingest_dir) to HDF files we want linked
            :param wps_dir: the WPS directory where we want the symlinks to appear
            :return:
            """
            for rel_path, hdf_name in zip(manifest, generate_hdf_names()):
                logging.info('Linking %s -> %s' % (osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, hdf_name)))
                symlink_unless_exists(osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, hdf_name))


class MODIS_TERRA(HDFSource):
    """
    750m data from the MODIS instrument on the Terra satellite
    """

    def __init__(self, ingest_dir):
        super(MODIS, self).__init__(ingest_dir)


    # do we need vtable stuff? seems specific to gribs.
    def vtables(self):
        """
        returns the variable tables that must be linked in for use with the MODIS data source.
        :return:
        """
        return {}

    def namelist_keys(self):
        """
        returns the namelist keys that must be modified in namelist.input with MODIS.
        (seems specific to gribs)
        """
        return{}

    def retrieve_hdfs(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrieve the files to satisfy the the simulation request from_utc = to_utc.

        Starts with the most recent cycle avaiable an hour ago, then moves further
        into the past. For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there. The presence of the remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local HDF files
        """
        # ensure seconds are zero, makes comparisons more straightforward later
        from_utc = from_utc.replace(second=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(second=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        # 'allData/6/MOD[DATA]/[YEAR]/[DAY]/MOD[DATA].A[YEAR][DAY].[HOUR][MINUTE].006.*.hdf'

        # compute_hdf_manifest function
        # requires from_utc and to_utc
        filepath = 'allData/6/MOD%02d/%04d/%03d/MOD%02d.A%04d%03d.%02d%02d.006.*.hdf'
        # Make sure that the minutes in the file name is a multiple of 5
        current = from_utc - datetime.timedelta(minutes=from_utc.minute % 5)

        hdfManifest = []
        while current <= to_utc:
            # Gives Julian day, used in the URL
            days = (current - datetime.datetime(current.year, 1,1)).days + 1
            year, hour, minute = current.year, current.hour, current.minute

            # Add url for both firedata and geolocation data
            # 'allData/6/MOD[DATA]/[YEAR]/[DAY]/MOD[DATA].A[YEAR][DAY].[HOUR][MINUTE].006.*.hdf'
            hdfManifest.append(filepath % (03, year, days, 03, year, days, hour, minute))
            hdfManifest.append(filepath % (14, year, days, 14, year, days, hour, minute))

            current = current + datetime.timedelta(minutes=5)


        # compute_meta_manifest function
        # requires from_utc and to_utc

        # Turns out that the metadata is stored in a strange path on the FTP server
        # Located at geoMeta/6/{TERRA, AQUA}/[YEAR]/{MOD, MYD}03_[YEAR]-[MONTH]-[DAY].txt
        # each file contains data for all of the MOD/MYD 03 granules collected on that day

        # Generate file paths for relevant metadata files
        filepath = 'geoMeta/6/TERRA/%04d/MOD03_%04d-%02d-%02d.txt'
        current = from_utc
        # Start at the beginning of the day so that the easy comparison is correct
        current = current.replace(second=0, minute=0, hour=0)
        metaManifest = []

        while current <= to_utc:
            metaManifest.append(filepath % (current.year, current.year, current.month, current.day))
            current = current + datetime.timedelta(days = 1)

        # Next step is to remove all locally available, and all that aren't in the lonlat box

        url_base = 'ftp://ladsweb.nascom.nasa.gov/'

        # This will all need to be done again for MODIS_AQUA and VIIRS
