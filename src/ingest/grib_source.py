# Copyright (C) 2013-2016 Martin Vejmelka, CU Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from utils import ensure_dir, symlink_unless_exists, timedelta_hours, readhead, Dict
from downloader import download_url, DownloadError
from datetime import datetime, timedelta
import pytz
import requests
import os
import os.path as osp
import sys
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
        self.ingest_dir = osp.abspath(js.get('ingest_path','ingest'))
        self.cache_dir = osp.abspath(js.get('cache_path','cache'))
        self.sys_dir = osp.abspath(js.get('sys_install_path',None))

        self.interval_seconds = 3600 * self.period_hours

    def prefix(self):
        """
        Return the prefix of ungribbed files 
        return: string  
        """
        return "COLMET"

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

        for vtable_id, vtable_path in vtables.iteritems():
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
        url = url_base + '/' + rel_path
        logging.info('downloading %s grib from %s' % (self.id, url))
        grib_path = osp.join(self.ingest_dir, rel_path)
        try:
            download_url(url, grib_path)
        except DownloadError as e:
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

    # instance variables  
    # id = "not specified"
    id = None
    period_hours = None
    forecast_cycle_hours = 6
    
    
class GribForecast(GribSource):
    """
    Common part for all grib forecasts.
    """

    def __init__(self, arg):
        super(GribForecast, self).__init__(arg)
        
    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None, cycle_start = None, download_whole_cycle=False):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned, or if retrieval fails, an error is raised.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files and hopefully cached COLMET files
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)
 
        logging.info('retrieve_gribs %s from_utc=%s to_utc=%s ref_utc=%s cycle_start=%s download_whole_cycle=%s' %
            (self.id, from_utc, to_utc, ref_utc, cycle_start, download_whole_cycle ))

        # it is possible that a cycle output is delayed and unavailable when we expect it (3 hours after cycle time)
        # in this case, the NAM grib source supports using previous cycles (up to 2)
        cycle_shift = 0
        while cycle_shift < 3:
    
            if cycle_start is not None:
                logging.info('forecast cycle start given as %s' % cycle_start)
            else:
                # select cycle (at least two hours behind real-time), out of the four NAM 218 cycles per day
                # which occurr at [0, 6, 12, 18] hours
                ref_utc_2 = ref_utc - timedelta(hours=3)
                ref_utc_2 = ref_utc_2.replace(minute=0,second=0,microsecond=0)
                cycle_start = min(from_utc, ref_utc_2)
                cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % 6)
                cycle_start -= timedelta(hours=6*cycle_shift)
                logging.info('forecast cycle start selected as %s' % cycle_start)

            if download_whole_cycle:
                logging.info('%s downloading whole cycle' % self.id)
                fc_start, fc_hours = 0, 84
            else:
                logging.info('%s downloading from %s to %s' % (self.id, from_utc, to_utc))
                fc_start, fc_hours = self.forecast_times(cycle_start, from_utc, to_utc)

            logging.info('%s downloading cycle %s forecast hours %d to %d' % (self.id, cycle_start, fc_start, fc_hours))

            # computes the relative paths of the desired files (the manifest)
            fc_list, colmet_files_utc = self.file_times(cycle_start,fc_start, fc_hours)
            grib_files, colmet_prefix, colmet_files = self.file_names(cycle_start, fc_list, colmet_files_utc)

            for f in grib_files:
               logging.info('NAM218 will retrive ' + f) 
            for f in colmet_files:
               logging.info('UNGRIB will create  ' + colmet_prefix + '/' +f) 

            # check what colmet files are available locally
            colmet_missing = [f for f in colmet_files if not osp.isfile(osp.join(self.cache_dir, colmet_prefix, f))] 
            logging.info('%d COLMET intermediate files not in cache' % len(colmet_missing) )
            for f in  colmet_missing:
                logging.info('Missing in cache    ' +f) 
            if len(colmet_missing) > 0:

                # check what's available locally
                nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), grib_files)
    
                # check if GRIBs we don't are available remotely
                url_base = self.remote_url
                logging.info('Retrieving NAM218 GRIBs from %s' % url_base)
                unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
                if len(unavailables) > 0:
                    logging.warning('NAM218: failed retrieving cycle data for cycle %s, unavailables %s' % (cycle_start, repr(unavailables)))
                    cycle_shift += 1
                    continue
    
                # download all gribs we need
                map(lambda x: self.download_grib(url_base, x), nonlocals)

            # return manifest

            return Dict({'grib_files': grib_files, 
                'colmet_files_utc': colmet_files_utc, 
                'colmet_prefix': colmet_prefix, 
                'colmet_files': colmet_files, 
                'colmet_missing': colmet_missing})

        raise GribError('Unsatisfiable: failed to retrieve GRIB2 files in eligible cycles %s' % repr(unavailables))



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



