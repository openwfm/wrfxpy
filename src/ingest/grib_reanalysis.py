from ingest.grib_source import GribError, GribSource
from utils import ensure_dir, symlink_unless_exists, timedelta_hours, readhead, Dict
from downloader import download_url, DownloadError
from datetime import datetime, timedelta
import pytz
import requests
import os
import os.path as osp
import sys
import logging


class GribReanalysis(GribSource):
    """
    Common part for all grib reanalysis products.
    """

    def __init__(self, arg):
        super(GribReanalysis, self).__init__(arg)
    
    def retrieve_gribs(self, from_utc, to_utc, ref_utc=None, cycle_start_utc = None, download_all_gribs = False):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: dictionary with
        	'grib_files': list of grib files available, 
                'colmet_files_utc': list of datetimes for the colmet files, 
                'colmet_prefix': string as colmet file prefix (directory names)
                'colmet_files': list of all colmet files needed, 
                'colmet_missing': those missing in the cache
        """

        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)
        to_utc = to_utc.replace(minute=0, second=0, tzinfo=pytz.UTC)

        # round start_utc down and end_utc up to period - reanalysis has no forecast cycles
        start_utc = from_utc.replace(hour = from_utc.hour - from_utc.hour % self.period_hours)
        end_utc = to_utc + timedelta(hours=self.period_hours)-timedelta(seconds=1)
        end_utc = end_utc.replace(hour=end_utc.hour - end_utc.hour % self.period_hours)

        if (start_utc < self.available_from_utc) | (end_utc > self.available_to_utc):
            logging.error('%s is available from %s to %s only' % (self.id, self.available_from_utc, self.available_to_utc))
            logging.info('Check %s for %s' % (self.info_url, self.info))
            raise GribError('Unsatisfiable: %s not available for the requested dates' % self.id)

        # compute the manifest here
        at_time = start_utc
        grib_files = []
        colmet_files_utc=[]
        while at_time <= end_utc:
            grib_files.append(self.make_relative_url(at_time))
            # logging.info('Adding to manifest input file %s' % self.make_relative_url(at_time))
            colmet_files_utc.append(at_time)
            at_time += timedelta(hours=self.period_hours)
        colmet_prefix = self.id
        colmet_files = self.colmet_files(colmet_files_utc)
        colmet_missing = self.colmet_missing(colmet_prefix,colmet_files)

        # if no missing cached files, we do not care about gribs
        if len(colmet_missing) > 0:
            # print 'grib_files = ' + str(grib_files)
            # check what's available locally
            nonlocals = filter(lambda x: not self.grib_available_locally(osp.join(self.ingest_dir, x)), grib_files)
            #print 'nonlocals = ' + str(nonlocals)
            # check if GRIBs we don't have are available remotely
            url_base = self.remote_url
            logging.info('Retrieving CFSR GRIBs from %s' % url_base)
            unavailables = filter(lambda x: readhead(url_base + '/' + x).status_code != 200, nonlocals)
            if len(unavailables) > 0:
                raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

            # download all gribs not available remotely
            map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return Dict({'grib_files': grib_files, 
            'colmet_files_utc': colmet_files_utc, 
            'colmet_prefix': colmet_prefix, 
            'colmet_files': colmet_files,
            'colmet_missing': colmet_missing}) 

# instance variables - need to be defined in the subclasses
    info = None
    remote_url = None
    period_hours = None
    cycle_hours = None
    info_url = None
    available_from_utc = None
    available_to_utc = None
