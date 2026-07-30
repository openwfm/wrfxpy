from ingest.grib_source import GribError, GribSource
from utils import timedelta_hours, readhead, Dict
from datetime import datetime, timedelta, timezone
import os.path as osp
import logging

class GribForecast(GribSource):
    """
    Common part for all grib forecast products.
    """

    def minimum_forecast_lead_hours(self):
        """First forecast hour allowed for a historical run."""
        return 0

    def cycle_forecast_hours(self, cycle_start):
        """Last forecast hour available from this cycle."""
        return self.max_forecast_hours

    def __init__(self, arg):
        super(GribForecast, self).__init__(arg)
        self.max_forecast_hours = self.grib_forecast_hours_periods[-1]['hours']
        
    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None, cycle_start = None, download_whole_cycle=False):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the newest eligible cycle.  For each candidate cycle, the
        filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on every server.  If every server returns 404, we try an older cycle.
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
        from_utc = from_utc.replace(minute=0, second=0, microsecond=0, tzinfo=timezone.utc)
        to_utc = to_utc.replace(minute=0, second=0, microsecond=0, tzinfo=timezone.utc)

        if ref_utc is None:
            ref_utc = datetime.now(timezone.utc)
 
        logging.info('retrieve_gribs %s from_utc=%s to_utc=%s ref_utc=%s cycle_start=%s download_whole_cycle=%s' %
            (self.id, from_utc, to_utc, ref_utc, cycle_start, download_whole_cycle ))

        explicit_cycle = cycle_start is not None
        if explicit_cycle:
            first_cycle = cycle_start.replace(minute=0, second=0, microsecond=0)
            minimum_lead = 0 if download_whole_cycle else self.minimum_forecast_lead_hours()
            if timedelta_hours(from_utc - first_cycle, False) < minimum_lead:
                raise GribError('%s cycle %s starts before f%02d'
                                % (self.id, first_cycle, minimum_lead))
        else:
            minimum_lead = 0 if download_whole_cycle else self.minimum_forecast_lead_hours()
            ref_utc_2 = ref_utc - timedelta(hours=self.hours_behind_real_time)
            ref_utc_2 = ref_utc_2.replace(minute=0, second=0, microsecond=0)
            first_cycle = min(from_utc - timedelta(hours=minimum_lead), ref_utc_2)
            first_cycle = first_cycle.replace(
                hour=first_cycle.hour - first_cycle.hour % self.cycle_hours)

        attempts = 1 if explicit_cycle else (
            3 if download_whole_cycle else self.cycle_search_attempts)
        last_missing = {}
        for cycle_shift in range(attempts):
            cycle_start = first_cycle - timedelta(hours=self.cycle_hours * cycle_shift)
            cycle_forecast_hours = self.cycle_forecast_hours(cycle_start)
            if not download_whole_cycle:
                requested_hours = timedelta_hours(to_utc - cycle_start)
                if not explicit_cycle and requested_hours > self.max_forecast_hours:
                    break
                if requested_hours > cycle_forecast_hours:
                    if explicit_cycle:
                        raise GribError(
                            '%s cycle %s is only available through f%02d'
                            % (self.id, cycle_start, cycle_forecast_hours))
                    logging.info(
                        '%s cycle %s ends at f%02d; trying an older cycle'
                        % (self.id, cycle_start, cycle_forecast_hours))
                    continue

            if explicit_cycle:
                logging.info('forecast cycle start given as %s' % cycle_start)
            else:
                logging.info('forecast cycle start selected as %s' % cycle_start)

            if download_whole_cycle:
                logging.info('%s downloading whole cycle' % self.id)
                fc_start, fc_hours = 0, cycle_forecast_hours
            else:
                logging.info('%s downloading from %s to %s' % (self.id, from_utc, to_utc))
                fc_start, fc_hours = self.forecast_times(cycle_start, from_utc, to_utc)

            logging.info('%s downloading cycle %s forecast hours %d to %d' % (self.id, cycle_start, fc_start, fc_hours))

            # computes the relative paths of the desired files (the manifest)
            fc_list, colmet_files_utc = self.file_times(cycle_start, fc_start, fc_hours)
            grib_files = self.file_names(cycle_start, fc_list)
            colmet_prefix, colmet_files = self.colmet_names(cycle_start, colmet_files_utc)

            for f in grib_files:
               logging.info('%s will retrive %s' % (self.id, f)) 

            colmet_missing = self.colmet_missing(colmet_prefix,colmet_files)
            if len(colmet_missing) > 0:

                # check what's available locally
                nonlocals = [x for x in grib_files if not self.grib_available_locally(osp.join(self.ingest_dir, x))]
    
                # Use one complete server copy. Only 404 permits an older cycle.
                url_bases = self.remote_url
                if isinstance(url_bases,str):
                    url_bases = [url_bases]
                selected_url = None
                missing = {}
                check_order = nonlocals[-1:] + nonlocals[:-1]
                for url_base in url_bases:
                    logging.info('Checking %s GRIBs at %s' % (self.id, url_base))
                    check_base = osp.dirname(self.browse_aws) if url_base[:5] == 's3://' else url_base
                    for path in check_order:
                        status = readhead(osp.join(check_base, path), msg_level=0).status_code
                        if status == 404:
                            missing[url_base] = path
                            break
                        if status != 200:
                            raise GribError('%s availability check returned %s for %s'
                                            % (self.id, status, osp.join(check_base, path)))
                    else:
                        selected_url = url_base
                        break

                if selected_url is None:
                    last_missing = missing
                    if explicit_cycle:
                        raise GribError('%s cycle %s is unavailable: %s'
                                        % (self.id, cycle_start, repr(missing)))
                    logging.warning('%s cycle %s is unavailable: %s'
                                    % (self.id, cycle_start, repr(missing)))
                    continue
    
                # download all gribs not available remotely
                if selected_url[:5] == 's3://':
                    self.download_grib_many(selected_url, nonlocals, workers=32)
                else:
                    list(map(lambda x: self.download_grib(selected_url, x), nonlocals))

            # return manifest
            return Dict({'grib_files': [osp.join(self.ingest_dir, x) for x in grib_files], 
                'colmet_prefix': colmet_prefix, 
                'colmet_files_utc': colmet_files_utc,
                'colmet_files': [osp.join(self.cache_dir, colmet_prefix, f) for f in colmet_files],
                'colmet_missing': [osp.join(self.cache_dir, colmet_prefix, f) for f in colmet_missing]})

        raise GribError('Unsatisfiable: no complete %s cycle: %s'
                        % (self.id, repr(last_missing)))

    # GribForecast instance variables
    hours_behind_real_time = 3     # choose forecast cycle at least this much behind
    cycle_search_attempts = 3      # try this cycle and two older cycles
    

    def forecast_times(self, cycle_start, from_utc, to_utc):  
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
            logging.warning('forecast from %s to %s is less than one hour' % (str(from_utc), str(to_utc)))

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
