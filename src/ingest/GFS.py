from __future__ import absolute_import
from ingest.grib_source import GribError
from ingest.grib_reanalysis import GribReanalysis
from datetime import datetime
import pytz

class GFS(GribReanalysis):
    """
    The GFS (Global Forecast System) grib source as provided by NOMADS.

    GFS a reanalysis product computed with a large delay (18 months at this time) and
    the conditions are encoded every 3 hours [0, 3, 6, 9, ..., 21] every day.
    """

    def __init__(self, arg):
        super(GFS, self).__init__(arg)


    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NARR data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.GFS',
                'metgrid_vtable': 'METGRID.TBL.GFS'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NARR.
        :return: a list of paths to local GRIB files
        """

        #GFS requires that ''num_metgrid_soil_levels'' is set to 4.
        return { 'domains' : { 'num_metgrid_levels' : 49,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the GFS GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        path_tmpl = '%04d%02d/%04d%02d%02d/gfsanl_4_%04d%02d%02d_%02d00_000.grb2'
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, mon, year, mon, day, year, mon, day, hour)

    # instance variables
    info_url = 'https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ncdc:C00634'
    info = "Global Forecast System (GFS)"
    remote_url = 'https://nomads.ncdc.noaa.gov/data/gfsanl'
    period_hours = 6
    id = "GFS"
    available_from_utc = datetime(2004,3,1,tzinfo=pytz.UTC)
    available_to_utc = datetime.now(pytz.UTC)

    # see also
    # https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs
    # https://developers.google.com/earth-engine/datasets/catalog/NOAA_GFS0P25
    # https://catalog.data.gov/dataset/noaa-ncep-global-forecast-system-gfs-atmospheric-model