from ingest.grib_source import GribForecast, GribError
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead

class NAM227(GribForecast):
    """
    The NAM (North American Mesoscale) 227 grib source as provided by NOMADS (5km resolution).
    """

    def __init__(self, arg):
        super(NAM227, self).__init__(arg)

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

    
#    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None, cycle_start_utc = None, download_all_gribs = False):
#        """
#        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.
#
#        Starts with the most recent cycle available an hour ago, then moves further
#        into the past.  For each candidate cycle, the filenames are computed, the local cache is
#        checked for files that are already there.  The presence of remaining files is checked
#        on server, if not available, we try an older cycle, if yes, download is attempted.
#        Once all files are downloaded, the manifest is returned, or if retrieval fails, an error is raised.
#
#        :param from_utc: forecast start time
#        :param to_utc: forecast end time
#        :return: a list of paths to local GRIB files
#        """
#        # ensure minutes and seconds are zero, simplifies arithmetic later
#        from_utc = from_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)
#        to_utc = to_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)
#
#        if ref_utc is None:
#            ref_utc = datetime.now(pytz.UTC)
#
#        # it is possible that a cycle output is delayed and unavailable when we expect it (3 hours after cycle time)
#        # in this case, the NAM grib source supports using previous cycles (up to 2)
#        cycle_shift = 0
#        while cycle_shift < 3:
#    
#            # select cycle (at least two hours behind real-time), out of the four NAM 227 cycles per day
#            # which occurr at [0, 6, 12, 18] hours
#            ref_utc_2 = ref_utc - timedelta(hours=3)
#            ref_utc_2 = ref_utc_2.replace(minute=0,second=0,microsecond=0)
#            cycle_start = min(from_utc, ref_utc_2)
#            cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % 6)
#            cycle_start -= timedelta(hours=6*cycle_shift)
#
#            # check if the request is even satisfiable
#            delta_end = to_utc - cycle_start
#            fc_hours = delta_end.days * 24 + int(delta_end.seconds / 3600)
#
#            if fc_hours > 60:
#                raise GribError('Unsatisfiable: NAM 227 initial and boundary conditions are only available for 60 hours.')
#
#            delta_start = from_utc - cycle_start
#            fc_start = delta_start.days * 24 + int(delta_start.seconds / 3600)
#
#            # computes the relative paths of the desired files (the manifest)
#            
#            manifest = self.compute_manifest(cycle_start, fc_start, fc_hours)
#
#            # check what's available locally
#            nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)
#
#            # check if GRIBs we don't are available remotely
#            url_base = self.remote_url
#            logging.info('Attempting to retrieve GRIBs from %s' % url_base)
#            unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
#            if len(unavailables) > 0:
#                logging.warning('NAM227: failed retrieving cycle data for cycle %s, unavailables %s' % (cycle_start, repr(unavailables)))
#                cycle_shift += 1
#                continue
#
#            # download all gribs we need
#            map(lambda x: self.download_grib(url_base, x), nonlocals)
#
#            # return manifest
#            return Dict({'grib_files': manifest})
#
#        raise GribError('Unsatisfiable: failed to retrieve GRIB2 files in eligible cycles %s' % repr(unavailables))
#
#    def compute_manifest(self, cycle_start, fc_start, fc_hours):
#        """
#        Computes the relative paths of required GRIB files.
#
#        
#        :param cycle_start: UTC time of cycle start
#        :param fc_start: index of first file we need
#        :param fc_hours: final forecast hour 
#        """
#        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.conusnest.hiresf%02d.tm00.grib2'
#        fc_seq = range(fc_start, 36) + range(max(fc_start, 39), 60, 3)
#        # get all time points up to fc_hours plus one (we must cover entire timespan)
#        fc_list = [x for x in fc_seq if x < fc_hours]
#        fc_list.append(fc_seq[len(fc_list)])
#        
#        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
#        return map(lambda x: path_tmpl % (year, mon, day, hour, x), fc_list)

    def file_names(self, cycle_start, fc_list, colmet_files_utc):
        """
        Computes the relative paths of required GRIB and COLMET files.
        Dependent on the grib source.
    
        :param cycle_start: UTC time of cycle start
        :param fc_list: list of hours in the cycle when forecast will be donwloaded
        :param colmet_files_utc: 
        """

        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.conusnest.hiresf%02d.tm00.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, x) for x in fc_list]

        colmet_prefix_tmpl = '%s.%04d%02d%02dt%02d'
        colmet_files_tmpl ='COLMET:%04d-%02d-%02d_%02d'

        colmet_prefix = colmet_prefix_tmpl % (self.id, cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour)
        colmet_files = [colmet_files_tmpl % (x.year, x.month, x.day, x.hour) for x in colmet_files_utc]
        
        return grib_files, colmet_prefix, colmet_files 


    # instance variables
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    id = "NAM227"
    info_url="http://www.nco.ncep.noaa.gov/pmb/products/nam"
    info_text="NAM NEST over CONUS (5 km Resolution - Grid 227)"
    #    NAM227 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
    #    every 3 hours, starting with 39 and ending with 60.
    # grib_forecast_hours_periods = [{'hours':36,'period':1},{'hours':60,'period':3}]
    grib_forecast_hours_periods = [{'hours':60,'period':3}]
    cycle_hours = 6
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    period_hours = 3    # for METGRID and WRF
 
