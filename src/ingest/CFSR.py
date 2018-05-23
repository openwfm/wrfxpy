from ingest.grib_source import GribSource, GribError
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead

class CFSR(GribSource):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    CFSR consists if two products, pressure and surface, this is container class for them to avoid copying
    Methods return value if common, otherwise None and must be specified in contained classes
    """

    def __init__(self, ingest_dir):
        super(CFSR, self).__init__(ingest_dir)

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

    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None, cycle_start_utc = None, download_all_gribs = False):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned, or if retrieval fails, an error is raised.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """

        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)

        # CFSRv2 is only available in 6 hour increments
        start_utc = from_utc.replace(hour = from_utc.hour - from_utc.hour % 6)
        end_utc = to_utc + timedelta(hours=5,minutes=59,seconds=59)
        end_utc = end_utc.replace(hour=end_utc.hour - end_utc.hour % 6)

        if (start_utc < datetime(2011,4,1,tzinfo=pytz.UTC)) | (end_utc > datetime.now(pytz.UTC)):
            logging.error('CFSRv2 is available after 04/01/2011 only')
            logging.info('Check https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-regional-reanalysis-narr')
            raise GribError('Unsatisfiable: CFSRv2 not available for the requested dates')

        # compute the manifest here
        at_time = start_utc
        manifest = []
        while at_time <= end_utc:
            manifest.append(self.make_relative_url(at_time))
            logging.info('Adding to manifest input file %s' % self.make_relative_url(at_time))
            at_time += timedelta(hours=6)

        print 'manifest = ' + str(manifest)
        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)
        print 'nonlocals = ' + str(nonlocals)

        # check if GRIBs we don't have are available remotely
        url_base = self.remote_url
        logging.info('Retrieving CFSR GRIBs from %s' % url_base)
        unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs not available remotely
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return Dict({'grib_files': manifest})

    period_hours = 6
    id = "CFSR"

class CFSR_P(CFSR):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, ingest_dir):
        super(CFSR_P, self).__init__(ingest_dir)

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
        path_tmpl = '%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.pgrbh00.grib2'
        
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, year, mon, year, mon, day, hour)
    # instance variables
    remote_url = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_pgbh'
    id = "CFSR_P"

class CFSR_S(CFSR):
    """
    The CFSRv2 (Climate Forecast System Reanalysis v2) grib source as provided by NOMADS.

    CFSRv2 is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product and the conditions are encoded every 6 hours [0, 6, 12, 18] every day.
    """

    def __init__(self, ingest_dir):
        super(CFSR_S, self).__init__(ingest_dir)


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
        path_tmpl = '%04d/%04d%02d/%04d%02d%02d/cdas1.t%02dz.sfluxgrbf00.grib2'

        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, year, mon, year, mon, day, hour)
    # instance variables
    remote_url = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_flxf'
    id = "CFSR_S"

