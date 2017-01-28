#
# Dalton Burke, UC Denver
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

    def download_grib(self, url_base, rel_path, max_retries=3):
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


class MODIS(HDFSource):
    """
    750m data from the MODIS instrument on the Terra/Aqua satellites
    """

    def __init__(self, ingest_dir):
        super(MODIS, self).__init__(ingest_dir)


    # do we need vtable stuff? seems specific to gribs.
