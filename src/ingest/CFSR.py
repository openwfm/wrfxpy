from ingest.grib_source import GribError
from grib_reanalysis import GribReanalysis
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead

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

    def vtables(self):
        return None

    def namelist_wps_keys(self):
        return None

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with CFSRv2.

        CFSRv2 requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels' : 38,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

 
    period_hours = 6
    id = "CFSR"
    info_text = "The CFSRv2 (Climate Forecast System Reanalysis v2)"
    info_url = "https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2"
    available_to_utc = datetime.now(pytz.UTC)
    available_from_utc = datetime(2011,4,1,tzinfo=pytz.UTC)


class CFSR_P(CFSR):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, js):
        super(CFSR_P, self).__init__(js)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the CFSRv2 data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.CFSR_press_pgbh06',
                'metgrid_vtable': 'METGRID.TBL.CFSR'}

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

        path_tmpls = ['%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.pgrbh00.grib2' % (year, year, mon, year, mon, day, hour),
                    '%04d/%04d%02d/%04d%02d%02d/%04d%02d%02dcdas1.t%02dz.pgrbh00.grib2' % (year, year, mon, year, mon, day, year, mon, day, hour)]

        return self.available_online(path_tmpls)

    # instance variables
    remote_url = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_pgbh'
    id = "CFSR_P"
    prefix = 'COLMET_P'

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

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the CFSRv2 data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.CFSR_sfc_flxf06',
                'metgrid_vtable': 'METGRID.TBL.CFSR'}

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the CFSRv2 GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour

        path_tmpls = ['%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.sfluxgrbf00.grib2' % (year, year, mon, year, mon, day, hour),
                    '%04d/%04d%02d/%04d%02d%02d/%04d%02d%02dcdas1.t%02dz.sfluxgrbf00.grib2' % (year, year, mon, year, mon, day, year, mon, day, hour)]

        return self.available_online(path_tmpls)

    # instance variables
    remote_url = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_flxf'
    id = "CFSR_S"
    prefix = 'COLMET_S'


