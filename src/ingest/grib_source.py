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
from downloader import download_url, DownloadError

from datetime import datetime, timedelta
import pytz
import requests
import os
import os.path as osp
import sys
import logging
from ingest_utils import readhead


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

    def __init__(self, ingest_dir):
        """
        Initialize grib source with ingest directory (where GRIB files are stored).

        :param ingest_dir: root of GRIB storage
        """
        self.ingest_dir = osp.abspath(ingest_dir)

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

    def retrieve_gribs(self, from_utc, to_utc, ref_utc = None):
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

    def download_grib(self, url_base, rel_path, max_retries=3):
        """
        Download a GRIB file from a GRIB service and stream to <rel_path> in ingest_dir.

        :param url_base: the base URL part of the GRIB service
        :param rel_path: the relative path of the file (w.r.t GRIB base url and w.r.t self.ingest_dir)
        :param max_retries: how many times we may retry to download the file
        """
        url = url_base + '/' + rel_path
        grib_path = osp.join(self.ingest_dir, rel_path)
        try:
            download_url(url, grib_path, max_retries)
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
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.HRRR',
                'metgrid_vtable': 'METGRID.TBL.HRRR'}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with HRRR.

        HRRR requires that ''num_metgrid_soil_levels'' is set to 8.
        """
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
        logging.info('Retrieving GRIBs from %s' % url_base)
        unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs we need
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def compute_manifest(self, cycle_start, fc_hours):
        """
        Computes the relative paths of required GRIB2 files.

        HRRR provides 16 GRIB2 files, one per hour and performs a cycle every hour.

        :param cycle_start: UTC time of cycle start
        :param fc_hours: final forecast hour 
        """
        path_tmpl = 'hrrr.%04d%02d%02d/hrrr.t%02dz.wrfprsf%02d.grib2'
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: path_tmpl % (year, mon, day, hour, x), range(fc_hours+1))

    # instance variables
    # remote_url = 'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod'
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/'
    period_hours = 1

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
        return {'geogrid_vtable': 'GEOGRID.TBL',
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

        # it is possible that a cycle output is delayed and unavailable when we expect it (3 hours after cycle time)
        # in this case, the NAM grib source supports using previous cycles (up to 2)
        cycle_shift = 0
        while cycle_shift < 3:
    
            # select cycle (at least two hours behind real-time), out of the four NAM 218 cycles per day
            # which occurr at [0, 6, 12, 18] hours
            ref_utc_2 = ref_utc - timedelta(hours=3)
            ref_utc_2 = ref_utc_2.replace(minute=0,second=0,microsecond=0)
            cycle_start = min(from_utc, ref_utc_2)
            cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % 6)
            cycle_start -= timedelta(hours=6*cycle_shift)

            # check if the request is even satisfiable
            delta_end = to_utc - cycle_start
            fc_hours = delta_end.days * 24 + int(delta_end.seconds / 3600)

            if fc_hours > 84:
                raise GribError('Unsatisfiable: NAM 218 initial and boundary conditions are only available for 84 hours.')

            delta_start = from_utc - cycle_start
            fc_start = delta_start.days * 24 + int(delta_start.seconds / 3600)

            logging.info('NAM218: from_utc=%s to_utc=%s attempting retrieval for fc [%d,%d] of cycle %s' %
                          (str(from_utc), str(to_utc), fc_start, fc_hours, str(cycle_start)))

            # computes the relative paths of the desired files (the manifest)
            manifest = self.compute_manifest(cycle_start, fc_start, fc_hours)

            # check what's available locally
            nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

            # check if GRIBs we don't are available remotely
            url_base = self.remote_url
            logging.info('Retrieving GRIBs from %s' % url_base)
            unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
            if len(unavailables) > 0:
                logging.warning('NAM218: failed retrieving cycle data for cycle %s, unavailables %s' % (cycle_start, repr(unavailables)))
                cycle_shift += 1
                continue

            # download all gribs we need
            map(lambda x: self.download_grib(url_base, x), nonlocals)

            # return manifest
            return manifest

        raise GribError('Unsatisfiable: failed to retrieve GRIB2 files in eligible cycles %s' % repr(unavailables))
    
    def compute_manifest(self, cycle_start, fc_start, fc_hours):
        """
        Computes the relative paths of required GRIB files.

        NAM218 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
        every 3 hours, starting with 39 and ending with 84.

        :param cycle_start: UTC time of cycle start
        :param fc_start: index of first file we need
        :param fc_hours: final forecast hour 
        """
        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.awphys%02d.tm00.grib2'
        fc_seq = range(fc_start, 36) + range(max(fc_start, 39), 85, 3)
        # get all time points up to fc_hours plus one (we must cover entire timespan)
        fc_list = [x for x in fc_seq if x < fc_hours]
        fc_list.append(fc_seq[len(fc_list)])
        
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: path_tmpl % (year, mon, day, hour, x), fc_list)

    # instance variables
    #remote_url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    period_hours = 3


class NAM227(GribSource):
    """
    The NAM (North American Mesoscale) 227 grib source as provided by NOMADS (5km resolution).
    """

    def __init__(self, ingest_dir):
        super(NAM227, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NAM data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable':'Vtable.NAM',
                'metgrid_vtable':'METGRID.TBL.NAM'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NAM.

        NAM 227 requires that ''num_metgrid_soil_levels'' is set to 8.
        """
        return { 'domains' : { 'num_metgrid_levels': 43, 'num_metgrid_soil_levels' : 4 }}

    
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

        # it is possible that a cycle output is delayed and unavailable when we expect it (3 hours after cycle time)
        # in this case, the NAM grib source supports using previous cycles (up to 2)
        cycle_shift = 0
        while cycle_shift < 3:
    
            # select cycle (at least two hours behind real-time), out of the four NAM 227 cycles per day
            # which occurr at [0, 6, 12, 18] hours
            ref_utc_2 = ref_utc - timedelta(hours=3)
            ref_utc_2 = ref_utc_2.replace(minute=0,second=0,microsecond=0)
            cycle_start = min(from_utc, ref_utc_2)
            cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % 6)
            cycle_start -= timedelta(hours=6*cycle_shift)

            # check if the request is even satisfiable
            delta_end = to_utc - cycle_start
            fc_hours = delta_end.days * 24 + int(delta_end.seconds / 3600)

            if fc_hours > 60:
                raise GribError('Unsatisfiable: NAM 227 initial and boundary conditions are only available for 60 hours.')

            delta_start = from_utc - cycle_start
            fc_start = delta_start.days * 24 + int(delta_start.seconds / 3600)

            logging.info('NAM227: from_utc=%s to_utc=%s attempting retrieval for fc [%d,%d] of cycle %s' %
                          (str(from_utc), str(to_utc), fc_start, fc_hours, str(cycle_start)))

            # computes the relative paths of the desired files (the manifest)
            manifest = self.compute_manifest(cycle_start, fc_start, fc_hours)

            # check what's available locally
            nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

            # check if GRIBs we don't are available remotely
            url_base = self.remote_url
            logging.info('Attempting to retrieve GRIBs from %s' % url_base)
            unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
            if len(unavailables) > 0:
                logging.warning('NAM227: failed retrieving cycle data for cycle %s, unavailables %s' % (cycle_start, repr(unavailables)))
                cycle_shift += 1
                continue

            # download all gribs we need
            map(lambda x: self.download_grib(url_base, x), nonlocals)

            # return manifest
            return manifest

        raise GribError('Unsatisfiable: failed to retrieve GRIB2 files in eligible cycles %s' % repr(unavailables))

    def compute_manifest(self, cycle_start, fc_start, fc_hours):
        """
        Computes the relative paths of required GRIB files.

        NAM227 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
        every 3 hours, starting with 39 and ending with 60.

        :param cycle_start: UTC time of cycle start
        :param fc_start: index of first file we need
        :param fc_hours: final forecast hour 
        """
        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.conusnest.hiresf%02d.tm00.grib2'
        fc_seq = range(fc_start, 36) + range(max(fc_start, 39), 60, 3)
        # get all time points up to fc_hours plus one (we must cover entire timespan)
        fc_list = [x for x in fc_seq if x < fc_hours]
        fc_list.append(fc_seq[len(fc_list)])
        
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: path_tmpl % (year, mon, day, hour, x), fc_list)

    # instance variables
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    period_hours = 3
# CFSR
class CFSR_P(GribSource):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, ingest_dir):
        super(CFSR_P, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the CFSRv2 data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.CFSR_press_pgbh06',
                'metgrid_vtable': 'METGRID.TBL.CFSR'}

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.wps with CFSR_P
        return: a dictionary of namelist entries
        """
        return { 'ungrib' : {'prefix': 'COLMET_P'},
                 'metgrid': {'COLMET_S':'COLMET_P'} 
               }

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with CFSRv2.

        CFSRv2 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels' : 38,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

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

        # CFSRv2 is only available in 6 hour increments
        start_utc = from_utc.replace(hour = from_utc.hour - from_utc.hour % 6)
        end_utc = to_utc + timedelta(hours=5,minutes=59,seconds=59)
        end_utc = end_utc.replace(hour=end_utc.hour - end_utc.hour % 6)

        if (start_utc < datetime(2011,4,1,tzinfo=pytz.UTC)) | (end_utc > datetime.now(pytz.UTC)):
            logging.error('CFSRv2 is available after 04/01/2011 only')
            logging.info('Check https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-regional-reanalysis-narr')
            raise GribError('Unsatisfiable: CFSRv2 not available for the requested dates')

        # compute the manifest here
        at_time = start_utc
        manifest = []
        while at_time <= end_utc:
            manifest.append(self.make_relative_url(at_time))
            logging.info('Adding to manifest input file %s' % self.make_relative_url(at_time))
            at_time += timedelta(hours=6)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't have are available remotely
        url_base = self.remote_url
        logging.info('Retrieving GRIBs from %s' % url_base)
        unavailables = filter(lambda x: requests.head(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs not available remotely
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the CFSRv2 GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        path_tmpl = '%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.pgrbh00.grib2'
        print "path_tmpl=",path_tmpl
        
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, year, mon, year, mon, day, hour)
        print "path_tmpl",path_tmpl
    # instance variables
    remote_url = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_pgbh'
    period_hours = 6

class CFSR_S(GribSource):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, ingest_dir):
        super(CFSR_S, self).__init__(ingest_dir)

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.wps with CFSR_S
        return: a dictionary of namelist entries
        """
        return { 'ungrib' : {'prefix': 'COLMET_S'},
                 'metgrid': {'COLMET_S':'COLMET_P'} 
               }

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the CFSRv2 data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.CFSR_sfc_flxf06',
                'metgrid_vtable': 'METGRID.TBL.CFSR'}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with CFSRv2.

        CFSRv2 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels' : 38,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

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

        # CFSRv2 is only available in 6 hour increments
        start_utc = from_utc.replace(hour = from_utc.hour - from_utc.hour % 6)
        end_utc = to_utc + timedelta(hours=5,minutes=59,seconds=59)
        end_utc = end_utc.replace(hour=end_utc.hour - end_utc.hour % 6)

        if (start_utc < datetime(2011,4,1,tzinfo=pytz.UTC)) | (end_utc > datetime.now(pytz.UTC)):
            logging.error('CFSRv2 is available after 04/01/2011 only')
            logging.info('Check https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-regional-reanalysis-narr')
            raise GribError('Unsatisfiable: CFSRv2 not available for the requested dates')

        # compute the manifest here
        at_time = start_utc
        manifest = []
        while at_time <= end_utc:
            manifest.append(self.make_relative_url(at_time))
            logging.info('Adding to manifest input file %s' % self.make_relative_url(at_time))
            at_time += timedelta(hours=6)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't have are available remotely
        url_base = self.remote_url
        logging.info('Retrieving GRIBs from %s' % url_base)
        unavailables = filter(lambda x: requests.head(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs not available remotely
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the CFSRv2 GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        path_tmpl = '%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.sfluxgrbf00.grib2'
        print "path_tmpl=",path_tmpl

        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, year, mon, year, mon, day, hour)
        print "path_tmpl",path_tmpl
    # instance variables
    remote_url = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_flxf'
    period_hours = 6


#end CFSR
class NARR(GribSource):
    """
    The NARR (North American Regional Reanalysis) grib source as provided by NOMADS.

    NARR is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product computed with a large delay (18 months at this time) and
    the conditions are encoded every 3 hours [0, 3, 6, 9, ..., 21] every day.
    """

    def __init__(self, ingest_dir):
        super(NARR, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NARR data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.NARR',
                'metgrid_vtable': 'METGRID.TBL.NARR'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NARR.

        NARR requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels' : 30,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

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

        # NARR is only available with a large delay and is available in 3 hour increments
        # current delay is about 18 months, so ref_utc is not even used here.
        start_utc = from_utc.replace(hour = from_utc.hour - from_utc.hour % 3)
        end_utc = to_utc + timedelta(hours=2,minutes=59,seconds=59)
        end_utc = end_utc.replace(hour=end_utc.hour - end_utc.hour % 3)

        if (start_utc < datetime(1979,1,1,tzinfo=pytz.UTC)) | (end_utc > datetime(2014,10,2,tzinfo=pytz.UTC)):
            logging.error('NARR is available 01Jan1979 - 02Oct2014 only')
            logging.info('Check https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-regional-reanalysis-narr')
            raise GribError('Unsatisfiable: NARR not available for the requested dates')

        # compute the manifest here
        at_time = start_utc
        manifest = []
        while at_time <= end_utc:
            manifest.append(self.make_relative_url(at_time))
            logging.info('Adding to manifest input file %s' % self.make_relative_url(at_time))
            at_time += timedelta(hours=3)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't have are available remotely
        url_base = self.remote_url
        logging.info('Retrieving GRIBs from %s' % url_base)
        unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs not available remotely
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the NARR GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        path_tmpl = '%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb'
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, mon, year, mon, day, year, mon, day, hour)

    # instance variables
    remote_url = 'http://nomads.ncdc.noaa.gov/data/narr'
    period_hours = 3

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


