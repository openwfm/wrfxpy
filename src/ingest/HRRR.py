from ingest.grib_source import GribSource, GribError
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp
from utils import Dict, timedelta_hours, readhead


class HRRR(GribSource):
    """
    The HRRR (High Resolution Rapid Refresh) grib source as provided by NOMADS.
    """

    def __init__(self, ingest_dir):
        super(HRRR, self).__init__(ingest_dir)

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

        HRRR requires that ''num_metgrid_soil_levels'' is set to 8.
        """
        return {}

    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None, cycle_start_utc = None, download_all_gribs = False):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)

        if ref_utc is None:
            ref_utc = datetime.now(pytz.UTC)

        # select cycle (at least one hour behind)
        cycle_start = min(from_utc, ref_utc - timedelta(hours=1))

        # check if the request is even satisfiable
        delta = to_utc - cycle_start
        fc_hours = delta.days * 24 + delta.seconds / 3600

        if fc_hours > 15:
            raise GribError('Unsatisfiable: HRRR only forecasts 15 hours ahead.')

        # computes the relative paths of the desired files (the manifest)
        manifest = self.compute_manifest(cycle_start, fc_hours)

        # check what's available locally
        nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't are available remotely
        url_base = self.remote_url
        logging.info('Retrieving HRRR GRIBs from %s' % url_base)
        unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs we need
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return Dict({'grib_files': manifest})

    def compute_manifest(self, cycle_start, fc_hours):
        """
        Computes the relative paths of required GRIB2 files.

        HRRR provides 16 GRIB2 files, one per hour and performs a cycle every hour.

        :param cycle_start: UTC time of cycle start
        :param fc_hours: final forecast hour 
        """
        path_tmpl = 'hrrr.%04d%02d%02d/hrrr.t%02dz.wrfprsf%02d.grib2'
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: path_tmpl % (year, mon, day, hour, x), range(fc_hours+1))

    # instance variables
    # remote_url = 'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod'
    id = "HRRR"
    remote_url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/'
    period_hours = 1

