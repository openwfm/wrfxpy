from ingest.grib_forecast import GribForecast


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

        NAM 218 requires that ''num_metgrid_soil_levels'' is set to 4.
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
            
        grib_files = [
            [
                'nam.{0:04d}{1:02d}{2:02d}/nam.t{3:02d}z.awphys{4:02d}.tm00.grib2'.format(
                    cycle_start.year, cycle_start.month, 
                    cycle_start.day, cycle_start.hour, x
                ),
                '{0:04d}{1:02d}/{0:04d}{1:02d}{2:02d}/nam_218_{0:04d}{1:02d}{2:02d}_{3:02d}00_{4:03d}.grb2'.format(
                    cycle_start.year, cycle_start.month, 
                    cycle_start.day, cycle_start.hour, x
                )
            ] 
            for x in fc_list
        ]

        return [self.available_online(grib_file) for grib_file in grib_files]

    # instance variables
    id = "NAM218"
    info_url = "https://www.nco.ncep.noaa.gov/pmb/products/nam/"
    info_aws = "https://registry.opendata.aws/noaa-nam/"
    info_text = "NAM 218 AWIPS Grid - CONUS (12-km Resolution; full complement of pressure level fields and some surface-based fields)"
    info = "North American Mesoscale (NAM) Forecast System Grid 218"
    remote_url = [
        "https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/", 
        "s3://noaa-nam-pds/", 
        "https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/forecast/",
        "https://www.ncei.noaa.gov/thredds/fileServer/model-nam218-old/"
    ]
    browse_aws = "https://noaa-nam-pds.s3.amazonaws.com/index.html"
    cycle_hours = 6
    period_hours = 3    # for METGRID and WRF
    #    NAM218 provides hourly GRIB2 files up to hour 36 and then one GRIB2 file
    #    every 3 hours, starting with 39 and ending with 84.
    grib_forecast_hours_periods = [{'hours':84, 'period':3}]


