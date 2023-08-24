from __future__ import absolute_import
from ingest.grib_source import GribError
from ingest.grib_forecast import GribForecast
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

        NAM 227 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels': 43, 'num_metgrid_soil_levels' : 4 }}


    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB files.
        Dependent on the grib source.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of hours in the cycle when forecast will be donwloaded
        """

        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.conusnest.hiresf%02d.tm00.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, x) for x in fc_list]

        return grib_files


    # instance variables
    id = "NAM227"
    info_url="https://www.nco.ncep.noaa.gov/pmb/products/nam"
    info_text="NAM NEST over CONUS (5 km Resolution - Grid 227)"
    info = "North American Mesoscale (NAM) Forecast System Grid 227"
    remote_url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/'
    cycle_hours = 6
    period_hours = 3    # for METGRID and WRF
    #    NAM227 provides hourly GRIB2 files up to 60.
    grib_forecast_hours_periods = [{'hours':60,'period':1}]

