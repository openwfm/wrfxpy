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

    

    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB files.
        Dependent on the grib source.

 
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

        
        return grib_files

    # instance variables
    id = "NAM218"
    #    NAM218 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
    #    every 3 hours, starting with 39 and ending with 84.
    # grib_forecast_hours_periods = [{'hours':36,'period':3} , {'hours':84,'period':3}]
    grib_forecast_hours_periods = [{'hours':84,'period':3}]
    cycle_hours = 6
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod'
    #period_hours = 1
    period_hours = 3    # for METGRID and WRF
    info_text = "NAM 218 AWIPS Grid - CONUS (12-km Resolution; full complement of pressure level fields and some surface-based fields)"
    info_url = "http://www.nco.ncep.noaa.gov/pmb/products/nam/"
 

