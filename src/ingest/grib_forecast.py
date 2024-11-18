from ingest.grib_source import GribError, GribSource
from utils import timedelta_hours, readhead, Dict
from datetime import datetime, timedelta
import pytz
import os.path as osp
import logging

class GribForecast(GribSource):
    """
    Common part for all grib forecast products.
    """

    def __init__(self, arg):
        super(GribForecast, self).__init__(arg)
        self.max_forecast_hours = self.grib_forecast_hours_periods[-1]['hours']
        
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
        :return: dictionary with
       	    'grib_files': list of grib files available, 
            'colmet_files_utc': list of datetimes for the colmet files, 
            'colmet_prefix': string as colmet file prefix, 
            'colmet_files': list of all colmet files, 
            'colmet_missing': list of colmet files that need to be created
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)
 
        logging.info('retrieve_gribs %s from_utc=%s to_utc=%s ref_utc=%s cycle_start=%s download_whole_cycle=%s' %
            (self.id, from_utc, to_utc, ref_utc, cycle_start, download_whole_cycle ))

        # it is possible that a cycle output is delayed and unavailable when we expect it (3 hours after cycle time)
        # in this case, the grib source supports using previous cycles (up to 2)
        cycle_shift = 0
        while cycle_shift < 3:
    
            if cycle_start is not None:
                cycle_start = cycle_start.replace(minute=0,second=0,microsecond=0)
                logging.info('forecast cycle start given as %s' % cycle_start)
            else:
                # select cycle (at least hours_behind_real_time behind)
                # for NAM218 which occurr at [0, 6, 12, 18] hours
                ref_utc_2 = ref_utc - timedelta(hours=self.hours_behind_real_time)
                ref_utc_2 = ref_utc_2.replace(minute=0,second=0,microsecond=0)
                cycle_start = min(from_utc, ref_utc_2)
                cycle_start = cycle_start.replace(hour = cycle_start.hour - cycle_start.hour % self.cycle_hours)
                cycle_start -= timedelta(hours=self.cycle_hours*cycle_shift)
                logging.info('forecast cycle start selected as %s' % cycle_start)

            if download_whole_cycle:
                logging.info('%s downloading whole cycle' % self.id)
                fc_start, fc_hours = 0, self.max_forecast_hours
            else:
                logging.info('%s downloading from %s to %s' % (self.id, from_utc, to_utc))
                fc_start, fc_hours = self.forecast_times(cycle_start, from_utc, to_utc)

            logging.info('%s downloading cycle %s forecast hours %d to %d' % (self.id, cycle_start, fc_start, fc_hours))

            # computes the relative paths of the desired files (the manifest)
            fc_list, colmet_files_utc = self.file_times(cycle_start,fc_start, fc_hours)
            grib_files = self.file_names(cycle_start, fc_list)
            colmet_prefix, colmet_files = self.colmet_names(cycle_start, colmet_files_utc)

            for f in grib_files:
               logging.info('%s will retrive %s' % (self.id, f)) 

            colmet_missing = self.colmet_missing(colmet_prefix,colmet_files)
            if len(colmet_missing) > 0:

                # check what's available locally
                nonlocals = [x for x in grib_files if not self.grib_available_locally(osp.join(self.ingest_dir, x))]
    
                # check if GRIBs we don't are available remotely
                url_bases = self.remote_url
                if isinstance(url_bases,str):
                    url_bases = [url_bases]
                for url_base in url_bases:
                    logging.info('Retrieving %s GRIBs from %s' % (self.id, url_base))
                    if url_base[:5] == 's3://':
                        unavailables = [x for x in nonlocals if readhead(osp.join(osp.dirname(self.browse_aws), x), msg_level=0).status_code != 200]
                        logging.debug(f"from s3: {unavailables}")
                    else:
                        unavailables = [x for x in nonlocals if readhead(osp.join(url_base, x), msg_level=0).status_code != 200]
                    if len(unavailables) == 0:
                        break
                if len(unavailables) > 0:
                    logging.warning('%s failed retrieving cycle data for cycle %s, unavailables %s'
                                         % (self.id, cycle_start, repr(unavailables)))
                    cycle_shift += 1
                    continue
    
                # download all gribs we need
                list(map(lambda x: self.download_grib(url_base, x), nonlocals))

            # return manifest

            return Dict({'grib_files': [osp.join(self.ingest_dir, x) for x in grib_files], 
                'colmet_files_utc': colmet_files_utc, 
                'colmet_prefix': colmet_prefix, 
                'colmet_files': colmet_files,
                'colmet_missing': colmet_missing})

        raise GribError('Unsatisfiable: failed to retrieve GRIB2 files in eligible cycles %s' % repr(unavailables))

    # GribForecast instance variables
    hours_behind_real_time = 3     # choose forecast cycle at least this much behind
    

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
            raise GribError('cycle start %s is after forecast start %s' % (str(cycle_start), str(from_utc)))
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

        
        g=self.grib_forecast_hours_periods
        fc_seq = [] 
        for i in range(0, len(g)):
            fc_seq += list(range(max(int(fc_start), 0 if i == 0 else g[i-1]['hours'] + g[i]['period']), 
            g[i]['hours'] + g[i]['period'], g[i]['period']))
        # get all time points up to fc_hours plus one (we must cover entire timespan)
        fc_list = [x for x in fc_seq if x < fc_hours]
        fc_list.append(fc_seq[len(fc_list)])

        colmet_files_utc = [cycle_start + timedelta(hours = x) for x in range(int(fc_start), fc_list[-1] +1, self.period_hours)]
  
        return fc_list, colmet_files_utc



    def colmet_names(self, cycle_start, colmet_files_utc):
        """
        Computes the relative paths of cached COLMET files.
    
        :param cycle_start: UTC time of cycle start
        :param colmet_files_utc: 
        """

        # met path: nam.YYYYMMDDtcc/COLMET:YYYY-MM-DD_hh
        # YYYYMMDD is the Year Month Day Hour of the cycle
        # YYYY-MM-DD_hh the Year Month Day Hour of the forecast

        colmet_prefix_tmpl = '%s.%04d%02d%02dt%02d'
        colmet_prefix = colmet_prefix_tmpl % (self.id, cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour)
        colmet_files = self.colmet_files(colmet_files_utc)
        
        return colmet_prefix, colmet_files 
