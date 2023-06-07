from __future__ import absolute_import
from ingest.grib_reanalysis import GribReanalysis
from datetime import datetime
import pytz

class CFSR(GribReanalysis):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    CFSR consists if two products, pressure and surface, this is container class for them to avoid copying
    Methods return value if common, otherwise None and must be specified in contained classes
    """

    def __init__(self, js):
        super(CFSR, self).__init__(js)

    def namelist_wps_keys(self):
        return None

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the CFSRv2 data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.CFSR',
                'metgrid_vtable': 'METGRID.TBL.CFSR'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with CFSRv2.

        CFSRv2 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels' : 38,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

 
    id = "CFSR"
    info_url = "https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2"
    info_aws = "https://registry.opendata.aws/noaa-cfs/"
    info_text = "The CFSRv2 (Climate Forecast System Reanalysis v2)"
    info = "The CFSRv2 (Climate Forecast System Reanalysis v2)"
    browse_aws = "https://noaa-cfs-pds.s3.amazonaws.com/"
    available_from_utc = datetime(2011,4,1,tzinfo=pytz.UTC)
    available_to_utc = datetime.now(pytz.UTC)
    period_hours = 6


class CFSR_P(CFSR):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, js):
        super(CFSR_P, self).__init__(js)

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.wps with CFSR_P
        return: a dictionary of namelist entries
        """
        return { 'ungrib' : {'prefix': 'COLMET_P'},
                 'metgrid': {'fg_name': ['COLMET_S','COLMET_P']} 
               }

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the CFSRv2 GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour

        path_tmpls = ['cdas.%04d%02d%02d/cdas1.t%02dz.pgrbh00.grib2' % (year, mon, day, hour),
                    '%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.pgrbh00.grib2' % (year, year, mon, year, mon, day, hour),
                    '%04d/%04d%02d/%04d%02d%02d/%04d%02d%02dcdas1.t%02dz.pgrbh00.grib2' % (year, year, mon, year, mon, day, year, mon, day, hour)]

        return self.available_online(path_tmpls)

    # instance variables
    id = "CFSR_P"
    prefix = 'COLMET_P'
    remote_url = ["s3://noaa-cfs-pds/", "https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-analysis/6-hourly-by-pressure"]

class CFSR_S(CFSR):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, js):
        super(CFSR_S, self).__init__(js)

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.wps with CFSR_S
        return: a dictionary of namelist entries
        """
        return { 'ungrib' : {'prefix': 'COLMET_S'},
                 'metgrid': {'fg_name': ['COLMET_S','COLMET_P']} 
               }

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the CFSRv2 GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour

        path_tmpls = ['cdas.%04d%02d%02d/cdas1.t%02dz.sfluxgrbf00.grib2' % (year, mon, day, hour),
                    '%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.sfluxgrbf00.grib2' % (year, year, mon, year, mon, day, hour),
                    '%04d/%04d%02d/%04d%02d%02d/%04d%02d%02dcdas1.t%02dz.sfluxgrbf00.grib2' % (year, year, mon, year, mon, day, year, mon, day, hour)]

        return self.available_online(path_tmpls)

    # instance variables
    id = "CFSR_S"
    prefix = 'COLMET_S'
    remote_url = ["s3://noaa-cfs-pds/", "https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-analysis/6-hourly-flux"]


