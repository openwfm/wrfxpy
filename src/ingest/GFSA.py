from __future__ import absolute_import
from ingest.grib_source import GribError
from ingest.grib_reanalysis import GribReanalysis
from datetime import datetime
import pytz


class GFSA(GribReanalysis):
    """
    The GFS (Global Forecast System) grib source as provided by NOMADS.
    
    The NCEP operational Global Forecast System analysis grids are on a 0.5 global latitude longitude grid.
    Model analysis runs occur at 00, 06, 12, and 18 UTC daily. 
    Grids include forecast time steps at a 3 hourly interval from 0 to 6.
    """

    def __init__(self, arg):
        super(GFSA, self).__init__(arg)


    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the GFS data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.GFS',
                'metgrid_vtable': 'METGRID.TBL.GFS'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with GFS.
        :return: a list of paths to local GRIB files
        """

        #GFS requires that ''num_metgrid_soil_levels'' is set to 4.
        return { 'domains' : { 'num_metgrid_levels' : 34,
                               'num_metgrid_soil_levels': 4,
                               'p_top_requested': 10000 }}

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the GFS GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour

        path_tmpls = ['%04d%02d/%04d%02d%02d/gfsanl_4_%04d%02d%02d_%02d00_000.grb2' % (year, mon, year, mon, day, year, mon, day, hour),
                      '%04d%02d/%04d%02d%02d/gfs_4_%04d%02d%02d_%02d00_000.grb2' % (year, mon, year, mon, day, year, mon, day, hour)]
        
        return self.available_online(path_tmpls)

    # instance variables
    id = "GFSA"
    info_url = 'https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ncdc:C00634'
    info_text = "Global Forecast System (GFS) Analysis"
    info = "Global Forecast System (GFS) Analysis"
    remote_url = ['https://www.ncei.noaa.gov/data/global-forecast-system/access/historical/analysis', 'https://www.ncei.noaa.gov/data/global-forecast-system/access/grid-004-0.5-degree/analysis']
    period_hours = 6
    available_from_utc = datetime(2004,3,1,tzinfo=pytz.UTC)
    available_to_utc = datetime.now(pytz.UTC)

    # see also
    # https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs
    # https://developers.google.com/earth-engine/datasets/catalog/NOAA_GFS0P25
    # https://catalog.data.gov/dataset/noaa-ncep-global-forecast-system-gfs-atmospheric-model
