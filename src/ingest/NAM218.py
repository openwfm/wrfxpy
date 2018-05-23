from ingest.grib_source import GribSource, GribError
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead


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

    def forecast_times(self,cycle_start, from_utc, to_utc):  
        """
        Compute the span of hours to be used in a forecast cycle
        This should be common to all forecast data sources

        :param cycle_start: UTC time of cycle start
        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return fc_start, fc_hours: first and last hour in the forecast to be used
        """

        logging.info('%s cycle %s forecast from %s to %s UTC' % (self.id, str( cycle_start), str( from_utc), str( to_utc)))

        # check if the request is even satisfiable
        if (from_utc - cycle_start).total_seconds() < 0:
            raise GribError('cycle start %s is larger than forecast start %s' % (str(cycle_start), str(from_utc)))
        if (to_utc - from_utc).total_seconds() < 3600:
            raise GribError('forecast from %s to %s is less than one hour' % (str(from_utc), str(to_utc)))

        fc_hours = timedelta_hours(to_utc - cycle_start)

        if fc_hours > self.max_forecast_hours :
            logging.error('cycle start %s to forecast end %s is more than %s hours' % (str(cycle_start), str(to_utc),self.max_forecast_hours))
            raise GribError('Unsatisfiable: %s forecast is only available for %s hours.' % (self.id, self.max_forecast_hours))

        fc_start = timedelta_hours(from_utc - cycle_start, False) # rounding down

        logging.info('%s using cycle %s hours %d to %d' % (self.id, str(cycle_start), fc_start, fc_hours))

        return fc_start, fc_hours


    def file_times(self, cycle_start, fc_start, fc_hours):
        """
        Computes the file times of required GRIB and COLMET files from start and end forecast hour
        This may depend on the forecast source
         
        NAM218 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
        every 3 hours, starting with 39 and ending with 84.

        :param cycle_start: UTC time of cycle start
        :param fc_start: index of first file we need
        :param fc_hours: final forecast hour 
        :return fc_list: hours from cycle_start for which is forecast required
        :return colmet_files_utc: utc time of files after ungrib
        """

        logging.info('period_hours = %d' % self.period_hours)
        if self.period_hours not in [1, 3]:
            raise GribError('period_hours = %d must be 1 or 3' % self.period_hours) 

        fc_seq = range(fc_start, 37, self.period_hours) + range(max(fc_start, 39), 85, 3)
        # get all time points up to fc_hours plus one (we must cover entire timespan)
        fc_list = [x for x in fc_seq if x < fc_hours]
        fc_list.append(fc_seq[len(fc_list)])

        colmet_files_utc = [cycle_start + timedelta(hours = x) for x in range(fc_start, fc_list[-1] +1, self.period_hours)]
  
        return fc_list, colmet_files_utc


    def file_names(self, cycle_start, fc_list, colmet_files_utc):
        """
        Computes the relative paths of required GRIB and COLMET files.
        Defintely dependent on the grib source.

        NAM218 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
        every 3 hours, starting with 39 and ending with 84.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of hours in the cycle when forecast will be donwloaded
        :param colmet_files_utc: 
        """

        # grib path: nam.YYYYMMDD/nam.tccz.awphysfh.tm00.grib2 
        #cc is the model cycle runtime (i.e. 00, 06, 12, 18) 
        #YYYYMMDD is the Year Month Day Hour of model runtime
        #fh is the forecast hour (i.e. 00, 03, 06, ..., 84) 

        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.awphys%02d.tm00.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, x) for x in fc_list]

        # met path: nam.YYYYMMDDtcc/COLMET:YYYY-MM-DD_hh
        # YYYYMMDD is the Year Month Day Hour of the cycle
        # YYYY-MM-DD_hh the Year Month Day Hour of the forecast

        colmet_prefix_tmpl = '%s.%04d%02d%02dt%02d'
        colmet_files_tmpl ='COLMET:%04d-%02d-%02d_%02d'

        colmet_prefix = colmet_prefix_tmpl % (self.id, cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour)
        colmet_files = [colmet_files_tmpl % (x.year, x.month, x.day, x.hour) for x in colmet_files_utc]
        
        return grib_files, colmet_prefix, colmet_files 

    # instance variables
    id = "NAM218"
    max_forecast_hours = 84
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    #period_hours = 1
    period_hours = 3
    info_text = "NAM 218 AWIPS Grid - CONUS (12-km Resolution; full complement of pressure level fields and some surface-based fields)"
    info_url = "http://www.nco.ncep.noaa.gov/pmb/products/nam/"
 

