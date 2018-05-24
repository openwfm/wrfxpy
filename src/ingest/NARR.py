from ingest.grib_source import GribError
from ingest.grib_reanalysis import GribReanalysis
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead


class NARR(GribReanalysis):
    """
    The NARR (North American Regional Reanalysis) grib source as provided by NOMADS.

    NARR a reanalysis product computed with a large delay (18 months at this time) and
    the conditions are encoded every 3 hours [0, 3, 6, 9, ..., 21] every day.
    """

    def __init__(self, arg):
        super(NARR, self).__init__(arg)


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
        :return: a list of paths to local GRIB files
        """

        #NARR requires that ''num_metgrid_soil_levels'' is set to 4.
        return { 'domains' : { 'num_metgrid_levels' : 30,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

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
    info_url = 'https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-regional-reanalysis-narr'
    info = "North American Regional Reanalysis (NARR)"
    remote_url = 'https://nomads.ncdc.noaa.gov/data/narr'
    period_hours = 3
    cycle_hours = 3
    id = "NARR"
    available_from_utc = datetime(1979,1,1,tzinfo=pytz.UTC)
    available_to_utc = datetime(2014,10,2,tzinfo=pytz.UTC)
