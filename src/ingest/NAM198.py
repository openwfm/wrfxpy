from ingest.grib_forecast import GribForecast


class NAM198(GribForecast):
    """
    The NAM (North American Mesoscale for Alaska) 198 grib source as provided by NOMADS (6km resolution).
    """

    def __init__(self, arg):
        super(NAM198, self).__init__(arg)

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

        NAM 198 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels': 43, 'num_metgrid_soil_levels' : 4 }}


    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB files.
        Dependent on the grib source.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of hours in the cycle when forecast will be donwloaded
        """

        path_tmpl = 'nam.%04d%02d%02d/nam.t%02dz.alaskanest.hiresf%02d.tm00.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, x) for x in fc_list]

        return grib_files


    # instance variables
    id = "NAM198"
    info_url = "https://www.nco.ncep.noaa.gov/pmb/products/nam"
    info_aws = "https://registry.opendata.aws/noaa-nam/"
    info_text = "NAM NEST over Alaska (6 km Resolution - Grid 198)"
    info = "North American Mesoscale (NAM) Forecast System Grid 198"
    remote_url = ["https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/"]
    #browse_aws = "https://noaa-nam-pds.s3.amazonaws.com/index.html"
    cycle_hours = 6
    period_hours = 3    # for METGRID and WRF
    #    NAM198 provides hourly GRIB2 files up to 60.
    grib_forecast_hours_periods = [{'hours':60, 'period':1}]

