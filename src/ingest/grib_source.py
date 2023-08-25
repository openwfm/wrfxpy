from utils import ensure_dir, symlink_unless_exists, readhead
from .downloader import download_url, DownloadError
import os.path as osp
import logging

class GribError(Exception):
    """
    Raised when a GribSource cannot retrieve GRIBs.
    """
    pass


class GribSource(object):
    """
    The parent class of all GRIB2 sources that implements common functionality, for example

    - local GRIB2 validation (file size check)
    - GRIB2 retrieval with retries (smart check whether server implements http-range)
    - symlinking GRIB2 files for ungrib
    """

    def __init__(self, js):
        """
        Initialize grib source with ingest directory (where GRIB files are stored).

        :param js: job structure with at least ingest_path root of GRIB storage and sys_install_path 
        """
        self.ingest_dir = osp.abspath(osp.join(js.get('ingest_path','ingest'),self.id))
        self.cache_dir = osp.abspath(js.get('cache_path','cache'))
        self.sys_dir = osp.abspath(js.get('sys_install_path',None))

        self.interval_seconds = 3600 * self.period_hours


    def colmet_files(self, colmet_files_utc):
        """
        Compute the names of list of colment files from their datetime_utc.
        :param colmet_files_utc: list of datetime uct times
        :return: file names
        """
        return ['%s:%04d-%02d-%02d_%02d' % (self.prefix, x.year, x.month, x.day, x.hour) for x in colmet_files_utc]

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified for this source
        return: a dictionary of namelist entries
        """
        return {}

    def vtables(self):
        """
        Returns the vtables that must be used with this source as a table with keys:
        geogrid_vtable, ungrib_table, metgrid_table.

        :return: a dictionary mapping vtable keys to specific table files
        """
        return {}

    def namelist_keys(self):
        """
        Some GRIB2 source files require that in namelist.input, certain parameters have
        particular values.  Such keys should be returned here.

        :return: a dictionary mapping section names to keys that must be modified.
        """
        return {}

    def clone_vtables(self, tgt):
        """
        Clone all vtables (build symlink name, ensure directories exist, create the symlink)
        :param tgt: target directory into which WPS is cloned
        """

        # where are the symlink locations for vtable files (name of symlink)
        vtable_locs = {'geogrid_vtable': 'geogrid/GEOGRID.TBL',
                        'ungrib_vtable': 'Vtable',
                       'metgrid_vtable': 'metgrid/METGRID.TBL'}
        vtables = self.vtables()
        # vtables: a dictionary with keys from list ['geogrid_vtable', 'ungrib_vtable', 'metgrid_vtable'],
        #               which contain paths of the variable tables relative to 'etc/vtables'

        for vtable_id, vtable_path in vtables.items():
            # build path to link location
            symlink_path = osp.join(tgt, vtable_locs[vtable_id])

            if not osp.exists(symlink_path):
                symlink_tgt = osp.join(self.sys_dir, "etc/vtables", vtable_path)
                symlink_unless_exists(symlink_tgt, ensure_dir(symlink_path))


    def retrieve_gribs(self, from_utc, to_utc, ref_utc = None, cycle_start_utc = None, download_all_gribs = False):
        """
        Attempts to retrieve the GRIB files for the forecast time range.
        It should be first verified whether the GRIB2 files are available locally.
        For any unavailable files, downloads should be initiated.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :param ref_utc: a reference time which defines 'now' for the purpose of
                        retrieval, None means datetime.utcnow().
        :return: a list of paths to local GRIB files
        """
        pass

    def download_grib(self, url_base, rel_path):
        """
        Download a GRIB file from a GRIB service and stream to <rel_path> in ingest_dir.

        :param url_base: the base URL part of the GRIB service
        :param rel_path: the relative path of the file (w.r.t GRIB base url and w.r.t self.ingest_dir)
        :param max_retries: how many times we may retry to download the file
        """
        url = osp.join(url_base,rel_path)
        logging.info('downloading %s grib from %s' % (self.id, url))
        grib_path = osp.join(self.ingest_dir, rel_path)
        try:
            download_url(url, grib_path)
        except DownloadError as e:
            logging.error('%s cannot download grib file %s' % (self.id, url))
            logging.warning('Please check %s for %s' % (self.info_url, self.info))
            raise GribError('GribSource: failed to download file %s' % url)

    
    def grib_available_locally(self, path):
        """
        Check if a GRIB2 file is available locally and if it's file size checks out.

        :param path: the GRIB2 file path
        """
        info_path = path + '.size' 
        if osp.exists(path) and osp.exists(info_path):
            content_size = int(open(info_path).read())
            return osp.getsize(path) == content_size
        else:
            return False


    def symlink_gribs(self, manifest, wps_dir):
        """
        Make symlinks in the form GRIBFILE.XYZ to all manifest files into wps_dir.

        :param manifest: relative paths (w.r.t. ingest_dir) to GRIB files we want linked
        :param wps_dir: the WPS directory where we want the symlinks to appear
        :return:
        """
        for rel_path, grib_name in zip(manifest, generate_grib_names()):
            logging.info('Linking %s -> %s' % ( osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, grib_name)) )
            symlink_unless_exists(osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, grib_name))

    def colmet_missing(self,colmet_prefix,colmet_files):
        """
        Make list of files missing in the cache
        :param colmet prefix: the cache subdirectory the files should be in
        :param colmet_files: List of all files needed
        :return: list of all files not in cache
        """

        logging.info('%s: %d COLMET intermediate files needed' % (self.id,len(colmet_files)) )
        for f in  colmet_files:
            logging.info('Will need file   ' +f)
        # check what colmet files are available locally
        colmet_missing = [f for f in colmet_files if not osp.isfile(osp.join(self.cache_dir, colmet_prefix, f))]
        logging.info('%s: %d COLMET intermediate files not in cache' % (self.id,len(colmet_missing)) )
        for f in  colmet_missing:
            logging.info('Missing in cache ' +f)
        return colmet_missing

    def available_online(self,links):
        """
        Make link available online from a list of links
        :param links: list of links
        :return: link available online, '' otherwise
        """

        url_bases = self.remote_url
        if isinstance(url_bases,str):
            url_bases = [url_bases]
        for url_base in url_bases:
            if url_base[:5] == 's3://':
                available = [x for x in links if readhead(osp.join(self.browse_aws, x), msg_level=0).status_code == 200]
            else:
                available = [x for x in links if readhead(osp.join(url_base, x), msg_level=0).status_code == 200]
            if len(available) > 0:
                self.remote_url = url_base
                return available[0]
        logging.error('grib_source.available_online - online file not existent, urls tried:\n {}'.format([osp.join(url_base,link) for url_base in url_bases for link in links]))
        raise GribError('GribSource: failed to find an available online file')


    # instance variables  
    prefix = 'COLMET'
    id = None
    period_hours = None    
    
    

## Utility functions

def generate_grib_names():
    """
    Keeps generating gribfile names from GRIBFILE.AAA to ZZZ.
    """
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for c1 in alphabet:
        for c2 in alphabet:
            for c3 in alphabet:
                yield "GRIBFILE." + c1 + c2 + c3
