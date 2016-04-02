# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
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

from utils import ensure_dir, symlink_unless_exists
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
    Parent class of all grib sources.
    """

    def __init__(self, ingest_dir):
        """
        Initialize grib source with ingest directory (where GRIB files are stored).

        :param ingest_dir: root of GRIB storage
        """
        self.ingest_dir = osp.abspath(ingest_dir)

    def vtables(self):
        """
        Returns the vtables that must be set for use with this source.
        """
        return {}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with this GRIB source.
        """
        return {}

    def retrieve_gribs(self, from_utc, to_utc, ref_utc = None):
        """
        Attempts to retrieve the GRIB files for the forecast range.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :param ref_utc: a reference time which defines 'now' for the purpose of retrieval, None results in ref_utc being replaced with the UTC time when requesting, so the retrieval functions simply 'do the right thing'
        :return: a list of paths to local GRIB files
        """
        pass

    def download_grib(self, url_base, rel_path, max_retries=3):
        """
        Download a GRIB file from a GRIB service and stream to rel_path in ingest_dir.

        :param url_base: the base URL part of the GRIB service
        :param rel_path: the relative path of the file (w.r.t GRIB base url and w.r.t self.ingest_dir)
        """
        url = url_base + '/' + rel_path
        r = requests.get(url, stream=True)
        content_size = int(r.headers['Content-Length'])

        # dump the correct file size to an info file next to the grib file
        # when re-using the GRIB2 file, we check its file size against this record
        # to avoid using partial files
        info_path = osp.join(self.ingest_dir, rel_path + '.size')
        open(ensure_dir(info_path), 'w').write(str(content_size))

        # stream the download to file
        grib_path = osp.join(self.ingest_dir, rel_path)
        with open(ensure_dir(grib_path), 'wb') as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)
        
        # does the server accept byte range queries? e.g. the NOMADs server does
        accepts_ranges = 'bytes' in r.headers.get('Accept-Ranges', '')
        
        file_size = osp.getsize(grib_path)
        retries_available = max_retries
        while file_size < content_size:
            if retries_available > 0:
                if accepts_ranges:
                    # if range queries are supported, try to download only the missing portion of the file
                    headers = { 'Range' : 'bytes=%d-%d' % (file_size, content_size) }
                    r = requests.get(url, headers=headers, stream=True)
                    with open(grib_path, 'ab') as f:
                        for chunk in r.iter_content(1024 * 1024):
                            f.write(chunk)
                    retries_available -= 1
                else:
                    # call the entire function recursively, this will attempt to redownload the entire file
                    # and overwrite previously downloaded data
                    self.download_grib(url_base, rel_path, max_retries-1)
            else:
                os.remove(grib_path)
                os.remove(info_path)
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
            symlink_unless_exists(osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, grib_name))


class HRRR(GribSource):
    """
    The HRRR (High Resolution Rapid Refresh) grib source as provided by NOMADS.
    """

    def __init__(self, ingest_dir):
        super(HRRR, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the HRRR data source.
        :return:
        """
        return {'geogrid_vtable': 'GEOGRID.TBL.HRRR',
                'ungrib_vtable': 'Vtable.HRRR',
                'metgrid_vtable': 'METGRID.TBL.HRRR'}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with HRRR.

        HRRR requires that ''num_metgrid_soil_levels'' is set to 8.
        """
        #return { 'domains' : { 'num_metgrid_soil_levels': 8 }}
        return {}

    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        # select cycle (at least one hour behind)
        cycle_start = min(from_utc, ref_utc - timedelta(hours=1))

        # check if the request is even satisfiable
        delta = to_utc - cycle_start
        fc_hours = delta.days * 24 + delta.seconds / 3600

        if fc_hours > 15:
            raise GribError('Unsatisfiable: HRRR only forecasts 15 hours ahead.')

        # computes the relative paths of the desired files (the manifest)
        manifest = self.compute_manifest(cycle_start, fc_hours)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't are available remotely
        url_base = self.remote_url
        unavailables = filter(lambda x: requests.head(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs we need
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def compute_manifest(self, cycle_start, fc_hours):
        """
        Computes the relative paths of required GRIB files.

        Note, the system is built so that relative paths are the same in local cache
        and in remote system w.r.t. URL base.

        :param cycle_start: UTC time of cycle start
        :param fc_hours: final forecast hour 
        """
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: self.path_tmpl % (year, mon, day, hour, x), range(fc_hours+1))

    # instance variables
    remote_url = 'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod'
    path_tmpl = 'hrrr.%04d%02d%02d/hrrr.t%02dz.wrfprsf%02d.grib2'


class NAM218(GribSource):
    """
    The NAM (North American Mesoscale) 218 grib source as provided by NOMADS.
    """

    def __init__(self, ingest_dir):
        super(NAM218, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NAM data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL.NAM',
                'ungrib_vtable':'Vtable.NAM',
                'metgrid_vtable':'METGRID.TBL.NAM'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NAM.

        NAM 218 requires that ''num_metgrid_soil_levels'' is set to 8.
        """
        return { 'domains' : { 'num_metgrid_levels': 40, 'num_metgrid_soil_levels' : 4 }}

    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned, or if retrieval fails, an error is raised.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        # select cycle (at least two hours behind real-time), out of the four NAM 218 cycles per day
        # which occurr at [0, 6, 12, 18] hours
        ref_utc_2 = ref_utc - timedelta(hours=2)
        ref_utc_2 = ref_utc_2.replace(minute=0,second=0,microsecond=0)
        cycle_start = min(from_utc, ref_utc_2)
        cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % 6)

        # check if the request is even satisfiable
        delta_end = to_utc - cycle_start
        fc_hours = delta_end.days * 24 + int(delta_end.seconds / 3600)

        if fc_hours > 84:
            raise GribError('Unsatisfiable: NAM 218 initial and boundary conditions are only available for 84 hours.')

        delta_start = from_utc - cycle_start
        fc_start = delta_start.days * 24 + int(delta_start.seconds / 3600)

        logging.info('NAM218: from_utc=%s to_utc=%s fc [%d,%d] cycle %s' % (str(from_utc), str(to_utc), fc_start, fc_hours, str(cycle_start)))

        # computes the relative paths of the desired files (the manifest)
        manifest = self.compute_manifest(cycle_start, fc_start, fc_hours)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't are available remotely
        url_base = self.remote_url
        unavailables = filter(lambda x: requests.head(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs we need
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def compute_manifest(self, cycle_start, fc_start, fc_hours):
        """
        Computes the relative paths of required GRIB files.

        Note, the system is built so that relative paths are the same in local cache
        and in remote system w.r.t. URL base.

        :param cycle_start: UTC time of cycle start
        :param fc_start: index of first file we need
        :param fc_hours: final forecast hour 
        """
        fc_seq = range(fc_start, 36) + range(max(fc_start, 39), 85, 3)
        # get all time points up to fc_hours plus one (we must cover entire timespan)
        fc_list = [x for x in fc_seq if x < fc_hours]
        fc_list.append(fc_seq[len(fc_list)])
        
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: self.path_tmpl % (year, mon, day, hour, x), fc_list)

    # instance variables
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.awphys%02d.grb2.tm00'



class NARR(GribSource):
    """
    The NARR (North American Regional Reanalysis) grib source as provided by NOMADS.
    """

    def __init__(self, ingest_dir):
        super(NARR, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NARR data source.

        :return: a dictionary of variable tables
        """
        return {'ungrib_vtable': 'Vtable.NARR'}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NARR.

        NAM 218 requires that ''num_metgrid_soil_levels'' is set to 8.
        """
        return {}

    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned, or if retrieval fails, an error is raised.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        # select cycle (at least two hours behind real-time), out of the four NAM 218 cycles per day
        # which occurr at [0, 6, 12, 18] hours
        cycle_start = min(from_utc, ref_utc - timedelta(hours=2))
        cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % 6)

        # check if the request is even satisfiable
        delta = to_utc - cycle_start
        fc_hours = delta.days * 24 + delta.seconds / 3600

        if fc_hours > 21:
            raise GribError('Unsatisfiable: NARR initial and boundary conditions are only available for 21 hours in one cycle.')

        # computes the relative paths of the desired files (the manifest)
        manifest = self.compute_manifest(cycle_start, fc_hours)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't are available remotely
        url_base = self.remote_url
        unavailables = filter(lambda x: requests.head(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs we need
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def compute_manifest(self, cycle_start, fc_hours):
        """
        Computes the relative paths of required GRIB files.

        Note, the system is built so that relative paths are the same in local cache
        and in remote system w.r.t. URL base.

        :param cycle_start: UTC time of cycle start
        :param fc_hours: final forecast hour 
        """
        fc_seq = range(0, 22, 3)
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        # get all time points up to fc_hours plus one (we must cover entire timespan)
        fc_list = [x for x in fc_seq if x < fc_hours]
        fc_list.append(fc_seq[len(fc_list)])
        return map(lambda x: self.path_tmpl % (year, mon, year, mon, day, year, mon, day, x), fc_list)

    # instance variables
    remote_url = 'http://nomads.ncdc.noaa.gov/data/narr'
    path_tmpl = '%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb'




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


