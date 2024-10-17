from ingest.grib_reanalysis import GribReanalysis
from datetime import datetime, timedelta
import pytz


class HRRRA(GribReanalysis):
    """
    The HRRR (High Resolution Rapid Refresh) grib source as provided by AWS.
    """

    def __init__(self, arg):
        super(HRRRA, self).__init__(arg)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the HRRR data source.
        :return:
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.HRRR',
                'metgrid_vtable': 'METGRID.TBL.HRRR'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with HRRR.

        HRRR requires that ''num_metgrid_soil_levels'' is set to 9.
        """
        return { 'domains' : { 'num_metgrid_levels': 41, 'num_metgrid_soil_levels' : 9 }}

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the GFS GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        utc_time -= timedelta(hours=3)
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        path_tmpls = [
            'hrrr.{:04d}{:02d}{:02d}/conus/hrrr.t{:02d}z.wrfprsf03.grib2'.format(
                year, mon, day, hour
            )
        ]

        return self.available_online(path_tmpls)

    # instance variables
    id = "HRRRA"
    info_url = "https://rapidrefresh.noaa.gov/hrrr"
    info_aws = "https://registry.opendata.aws/noaa-hrrr-pds/"
    info_text = "NOAA HRRR 3-km CONUS High-Resolution Rapid Refresh Forecast"
    info = "The High-Resolution Rapid Refresh (HRRR)"
    remote_url = "s3://noaa-hrrr-bdp-pds/"
    browse_aws = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/"
    period_hours = 1
    # HRRR provides hourly GRIB2 files up to hour 48.
    grib_forecast_hours_periods = [{'hours':48, 'period':1}]
    # more general info: https://rapidrefresh.noaa.gov/internal/pdfs/RAPX_HRRRX_NWS-13sep2016-pub.pdf
    # file content: http://www.nco.ncep.noaa.gov/pmb/products/hrrr/hrrr.t00z.wrfprsf00.grib2.shtml
    available_from_utc = datetime(2011,4,1,tzinfo=pytz.UTC)
    available_to_utc = datetime.now(pytz.UTC) + timedelta(hours=3)

class HRRRA_S(HRRRA):
    """
    The HRRR (High Resolution Rapid Refresh) grib source as provided by AWS or NOMADS.
    The 2D surface product.
    """

    def __init__(self, js):
        super(HRRRA_S, self).__init__(js)

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the GFS GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        utc_time -= timedelta(hours=3)
        
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        path_tmpls = [
            'hrrr.{:04d}{:02d}{:02d}/conus/hrrr.t{:02d}z.wrfsfcf03.grib2'.format(
                year, mon, day, hour
            )
        ]

        return self.available_online(path_tmpls)
    