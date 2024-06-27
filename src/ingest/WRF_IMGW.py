from __future__ import absolute_import
from ingest.grib_source import GribError
from ingest.grib_forecast import GribForecast
from datetime import datetime
import pytz



class WRF_IMGW(GribForecast):
    """
    The NAM (North American Mesoscale) 218 grib source as provided by NOMADS.
    """

    def __init__(self, arg):
        super(WRF_IMGW, self).__init__(arg)


    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NAM data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable':'Vtable.WRF_IMGW',
                'metgrid_vtable':'METGRID.TBL.WRF_IMGW'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NAM.

        NAM 218 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels': 47, 'num_metgrid_soil_levels' : 4 }}



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


        path_tmpl = 'METEOPG_GFS_d02_%04d-%02d-%02d_%02d.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, x) for x in fc_list]

        return grib_files

    # instance variables
    id = "WRF_IMGW"
    info_url = "W.I.P"
    info_text = "W.I.P"
    info = "W.I.P"
    remote_url = 'https://home.chpc.utah.edu/~u0631741/wrfgrib/' #Temporary path
    cycle_hours = 6
    period_hours = 1    # for METGRID and WRF
    #    NAM218 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
    #    every 3 hours, starting with 39 and ending with 84.
    grib_forecast_hours_periods = [{'hours':60,'period':1}]


