from ingest.grib_source import GribSource, GribError
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead



class NARR(GribSource):
    """
    The NARR (North American Regional Reanalysis) grib source as provided by NOMADS.

    NARR is different from other GRIB2 sources since it does not really have cycles.
    It's a reanalysis product computed with a large delay (18 months at this time) and
    the conditions are encoded every 3 hours [0, 3, 6, 9, ..., 21] every day.
    """

    def __init__(self, ingest_dir):
        super(NARR, self).__init__(ingest_dir)


    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the NARR data source.

        :return: a dictionary of variable tables
        """
        return {'geogrid_vtable': 'GEOGRID.TBL',
                'ungrib_vtable': 'Vtable.NARR',
                'metgrid_vtable': 'METGRID.TBL.NARR'}


    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with NARR.

        NARR requires that ''num_metgrid_soil_levels'' is set to 4.
        """
        return { 'domains' : { 'num_metgrid_levels' : 30,
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

        # NARR is only available with a large delay and is available in 3 hour increments
        # current delay is about 18 months, so ref_utc is not even used here.
        start_utc = from_utc.replace(hour = from_utc.hour - from_utc.hour % 3)
        end_utc = to_utc + timedelta(hours=2,minutes=59,seconds=59)
        end_utc = end_utc.replace(hour=end_utc.hour - end_utc.hour % 3)

        if (start_utc < datetime(1979,1,1,tzinfo=pytz.UTC)) | (end_utc > datetime(2014,10,2,tzinfo=pytz.UTC)):
            logging.error('NARR is available 01Jan1979 - 02Oct2014 only')
            logging.info('Check https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-regional-reanalysis-narr')
            raise GribError('Unsatisfiable: NARR not available for the requested dates')

        # compute the manifest here
        at_time = start_utc
        manifest = []
        while at_time <= end_utc:
            manifest.append(self.make_relative_url(at_time))
            logging.info('Adding to manifest input file %s' % self.make_relative_url(at_time))
            at_time += timedelta(hours=3)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't have are available remotely
        url_base = self.remote_url
        logging.info('Retrieving NARR GRIBs from %s' % url_base)
        unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs not available remotely
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return Dict({'grib_files': manifest})

    def make_relative_url(self, utc_time):
        """
        Build the relative URL of the NARR GRIB2 file, which is based on the UTC time.

        :param utc_time: the UTC time
        :return: the relative URL
        """
        path_tmpl = '%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb'
        year, mon, day, hour = utc_time.year, utc_time.month, utc_time.day, utc_time.hour
        return path_tmpl % (year, mon, year, mon, day, year, mon, day, hour)

    # instance variables
    remote_url = 'http://nomads.ncdc.noaa.gov/data/narr'
    period_hours = 3
    id = "NARR"
