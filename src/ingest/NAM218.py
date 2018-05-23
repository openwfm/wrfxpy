from ingest.grib_source import GribForecast, GribError
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead


class NAM218(GribForecast):
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
 

