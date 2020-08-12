from __future__ import absolute_import
from ingest.grib_source import GribError
from ingest.grib_forecast import GribForecast
from datetime import datetime
import pytz


class GFSF(GribForecast):
    """
    The GFS (Global Forecast System) grib source as provided by NOMADS.
    
    The NCEP operational Global Forecast System forecast grids are on a 0.25 global latitude longitude grid.
    Model forecast runs occur at 00, 06, 12, and 18 UTC daily. 
    Grids include forecast time steps at a hourly interval from 0 to 120, and a 3 hourly interval from 120 to 384. 
    """

    def __init__(self, arg):
        super(GFSF, self).__init__(arg)

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

        return { 'domains' : { 'num_metgrid_levels' : 34,
                               'num_metgrid_soil_levels' : 4,
                               'p_top_requested' : 10000 }}

    def file_names(self, cycle_start, fc_list):
        return None

    # instance variables
    id = "GFSF"
    info_url = "https://www.nco.ncep.noaa.gov/pmb/products/gfs/#GFS"
    info_text = "NCEP GFS 0.25 Degree Global Forecast Grids Historical Archive"
    info = "Global Forecast System (GFS) Forecast"
    remote_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
    cycle_hours = 6
    period_hours = 3 # for METGRID and WRF
    #    GFS provides hourly GRIB2 files up to hour 120 and then one GRIB2 file
    #    every 3 hours, starting with 123 and ending with 384.
    grib_forecast_hours_periods = [{'hours':384,'period':3}]
    # more information: https://www.nco.ncep.noaa.gov/pmb/products/gfs/nomads/


class GFSF_P(GFSF):
    """
    The GFS (Global Forecast System) grib source as provided by NOMADS. Pressure fields.
    
    The NCEP operational Global Forecast System forecast grids are on a 0.25 global latitude longitude grid.
    Model forecast runs occur at 00, 06, 12, and 18 UTC daily. 
    Grids include forecast time steps at a hourly interval from 0 to 120, and a 3 hourly interval from 120 to 384. 
    """

    def __init__(self, js):
        super(GFSF_P, self).__init__(js)

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.wps with GFS_P
        return: a dictionary of namelist entries
        """
        return { 'ungrib' : {'prefix': 'COLMET_P'},
                 'metgrid': {'fg_name': ['COLMET_S','COLMET_P']} 
               }

    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB files.
        Dependent on the grib source.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of hours in the cycle when forecast will be donwloaded
        """

        path_tmpl = 'gfs.%04d%02d%02d/%02d/gfs.t%02dz.pgrb2.0p25.f%03d'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, cycle_start.hour, x) for x in fc_list]

        return grib_files

    # instance variables
    id = "GFSF_P"
    prefix = 'COLMET_P'


class GFSF_S(GFSF):
    """
    The GFS (Global Forecast System) grib source as provided by NOMADS. Surface fields.
    
    The NCEP operational Global Forecast System forecast grids are on a 0.25 global latitude longitude grid.
    Model forecast runs occur at 00, 06, 12, and 18 UTC daily. 
    Grids include forecast time steps at a hourly interval from 0 to 120, and a 3 hourly interval from 120 to 384. 
    """

    def __init__(self, js):
        super(GFSF_S, self).__init__(js)

    def namelist_wps_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.wps with GFS_S
        return: a dictionary of namelist entries
        """
        return { 'ungrib' : {'prefix': 'COLMET_S'},
                 'metgrid': {'fg_name': ['COLMET_S','COLMET_P']} 
               }

    def file_names(self, cycle_start, fc_list):
        """
        Computes the relative paths of required GRIB files.
        Dependent on the grib source.

        :param cycle_start: UTC time of cycle start
        :param fc_list: list of hours in the cycle when forecast will be donwloaded
        """
        path_tmpl = 'gfs.%04d%02d%02d/%02d/gfs.t%02dz.sfluxgrbf%03d.grib2'
        grib_files = [path_tmpl % (cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour, cycle_start.hour, x) for x in fc_list]

        return grib_files

    # instance variables
    id = "GFSF_S"
    prefix = 'COLMET_S'
