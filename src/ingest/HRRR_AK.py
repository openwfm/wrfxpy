from ingest.grib_forecast import GribForecast


class HRRR_AK(GribForecast):
    """
    The HRRR (High Resolution Rapid Refresh) grib source as provided by NOMADS.
    """

    def __init__(self, arg):
        super(HRRR_AK, self).__init__(arg)

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

    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB2 files.

        HRRR_AK provides 16 GRIB2 files, one per hour and performs a cycle every hour.

        :param cycle_start: UTC time of cycle start
        :param fc_hours: final forecast hour 
        """
        path_tmpl = 'hrrr.%04d%02d%02d/alaska/hrrr.t%02dz.wrfprsf%02d.ak.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, x) for x in fc_list]

        return grib_files

    # instance variables
    id = "HRRR_AK"
    info_url = "https://rapidrefresh.noaa.gov/hrrr"
    info_aws = "https://registry.opendata.aws/noaa-hrrr-pds/"
    info_text = "NOAA HRRR Alaska 3-km CONUS High-Resolution Rapid Refresh Forecast"
    info = "The High-Resolution Rapid Refresh (HRRR)"
    #S3 first, NOMADS as fallback -- the pattern HRRR and NAM218 already use.  Both
    #lines were commented out together in 835a893 (2024-05-26), the commit that added
    #this source, and the reason was almost certainly the coupling rather than an S3
    #outage: grib_source tests s3:// availability through browse_aws, so a remote_url
    #added without it dereferences a missing attribute.  The data was on S3 that very
    #day.  NOMADS alone retains about 1 day for HRRR Alaska, which made anything older irreproducible;
    #verified 2026-10-05 by retrieving 233 MB per lead through this code path for a date 30 days
    #back, where NOMADS has nothing.
    remote_url = ["s3://noaa-hrrr-bdp-pds/",
                  "https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/"]
    browse_aws = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/index.html"
    #cycle_hours = 6, not 1, and for two reasons measured against the bucket on
    #2026-10-05.  HRRR Alaska does not run hourly -- it runs every 3 h, so t01z is a 404
    #outright -- and of those, only 00/06/12/18 reach the f48 declared below, while
    #03/09/15/21 stop at f18 (f18 is 200, f20 is 404).  With cycle_hours = 1 the selector
    #asked for cycles that do not exist two hours in three, and for leads the remaining
    #ones do not carry.  Same defect as RAP and RRFS_CONUS, fixed the same way; HRRR
    #CONUS already uses 6 for the f48 half of it.  Drop to 3 only together with
    #grib_forecast_hours_periods hours = 18.
    cycle_hours = 6
    period_hours = 1
    hours_behind_real_time = 1     # choose forecast cycle at least one hour behind
    # HRRR provides hourly GRIB2 files up to hour 48.
    grib_forecast_hours_periods = [{'hours':48, 'period':1}]
    # more general info: https://rapidrefresh.noaa.gov/internal/pdfs/RAPX_HRRRX_NWS-13sep2016-pub.pdf
    # file content: http://www.nco.ncep.noaa.gov/pmb/products/hrrr/hrrr.t00z.wrfprsf00.grib2.shtml

