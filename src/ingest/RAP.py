from ingest.grib_forecast import GribForecast


class RAP(GribForecast):
    """
    The RAP (Rapid Refresh) forecast grib source as provided by AWS or NOMADS.

    This class uses the full-domain pressure-level RAP product:
        rap.tHHz.wrfprsfFF.grib2
    """

    def __init__(self, arg):
        super(RAP, self).__init__(arg)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the RAP data source.
        :return:
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.RAP',
                'metgrid_vtable': 'METGRID.TBL.RAP'}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with RAP.

        RAP wrfprs is a pressure-level WRF-oriented product with 38 metgrid levels and
        9 soil levels.
        """
        return {'domains': {'num_metgrid_levels': 41, 'num_metgrid_soil_levels': 9}}

    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB2 files.

        RAP provides one GRIB2 file per forecast hour and performs a cycle every hour.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of forecast hours
        """
        path_tmpl = 'rap.%04d%02d%02d/rap.t%02dz.wrfprsf%02d.grib2'
        grib_files = [path_tmpl % (
                cycle_start.year, cycle_start.month,
                cycle_start.day, cycle_start.hour, x
            ) for x in fc_list
        ]

        return grib_files

    # instance variables
    id = "RAP"
    info_url = "https://rapidrefresh.noaa.gov/"
    info_aws = "https://registry.opendata.aws/noaa-rap/"
    info_text = "NOAA RAP 13-km North America Rapid Refresh Forecast"
    info = "The Rapid Refresh (RAP)"
    remote_url = [
        "s3://noaa-rap-pds/",
        "https://nomads.ncep.noaa.gov/pub/data/nccf/com/rap/prod/"
    ]
    browse_aws = "https://noaa-rap-pds.s3.amazonaws.com/index.html"
    cycle_hours = 1
    period_hours = 1
    hours_behind_real_time = 1

    # RAP provides hourly forecasts.
    # Standard cycles go to FH21.
    # Extended cycles 03, 09, 15, and 21 UTC go to FH51.
    grib_forecast_hours_periods = [{'hours': 51, 'period': 1}]