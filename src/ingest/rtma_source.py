# Copyright (C) 2013-2016 Martin Vejmelka, CU Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from downloader import download_url, DownloadError
from utils import esmf_to_utc

import requests
from datetime import datetime, timedelta
import pytz
import logging
import os.path as osp


class RTMA(object):
    """
    This class supports the ingest of NOAA RTMA (Real-time Mesoscale Analysis) data.

    RTMA is different since there is only one set of files for 24 hours (1 fileset
    per hour), which are continuously overwritten.  This RTMA class thus checks the
    filetime for each file to verify it belongs to 'today' instead of the day before.

    Note that RTMA operates in the UTC time zone.
    """

    def __init__(self, ingest_path, var_list = ['temp', 'td', 'precipa', 'wspd', 'wdir']):
        """
        Initialize the RTMA object with variables to track.

        :param var_list: variables to download
        :param ingest_path: path to the ingest directory
        """
        self.var_list = var_list
        self.ingest_path = ingest_path


    def retrieve_rtma(self, cycle):
        """
        Attempts to retrieve the variables passed in during initialization.
        Any files already downloaded are not modified.
        Returns a list of variables that have not been downloaded.

        :param cycle: the cycle (UTC) for which to retrieve the RTMA
        :return: tuple with list of all variables that are not ready yet
                 and dictonary with path to stored files
        """
        ts = cycle.replace(minute=0, second=0, microsecond=0)
        logging.info('RTMA retrieving variables %s for cycle %s.' % (self.var_list, str(ts)))

        vars_paths = map(lambda x: (x, self._local_var_path(ts, x)), self.var_list)
        ready = dict(filter(lambda x: self._is_var_cached(x[1]), vars_paths))
        nonlocals = filter(lambda x: not self._is_var_cached(x[1]), vars_paths)
        if nonlocals:
            nl_vars = [x[0] for x in nonlocals]
            logging.info('RTMA variables %s are not available locally, trying to download.' % nl_vars)

        not_ready = []
        for var, local_path in nonlocals:
            if self._is_var_ready(ts, var):
                download_url(self._remote_var_url(cycle.hour, var), local_path)
                ready[var] = local_path
            else:
                not_ready.append(var)

        if not_ready:
            logging.info('RTMA the variables %s for hour %d are not ready.' % (not_ready, cycle.hour))
            # unless a file was downloaded, it makes no sense to check the server immediately again
        else:
            # if all files are available, return
            logging.info('RTMA success obtaining variables %s for hour %d.' % (self.var_list, cycle.hour))

        return not_ready, ready


   
    def _local_var_path(self, ts, var):
        """
        Build the path to the file storing variable var at timestamp ts (only hours taken into account).
        This path is relative to the ingest_path specified at init time.

        :param ts: the timestamp
        :param var: the variable name
        """
        rel_path = 'rtma/%04d%02d%02d/%02d/%s.grib' % (ts.year, ts.month, ts.day, ts.hour, var)
        return osp.join(self.ingest_path, rel_path)


    def _is_var_ready(self, cycle, var):
        """
        Checks if the variable var is ready for the given forecast hour by comparing its
        filetime to the timestamp given by the forecast hour.  If the filetime is newer
        (later) then the variable is ready.

        :param cycle: which cycle are we working with (UTC)
        :param var: the variable identifier
        :return: true if the variable is ready
        """
        # find last-modified time of file in UTC timezone
        url = self._remote_var_url(cycle.hour, var)
        r = requests.head(url)
        if r.status_code != 200:
            raise ValueError('Cannot find variable %s for hour %d at url %s' % (var, cycle.hour, url))
        last_modif = self._parse_header_timestamp(r.headers['Last-Modified'])

        return last_modif > cycle


    @staticmethod
    def _remote_var_url(hour, var):
        """
        Build the URL pointing to the file to be retrieved.

        :param hour: the forecast hour
        :param var: the variable to download
        """
        # rtma_base = 'http://weather.noaa.gov/pub/SL.us008001/ST.opnl/DF.gr2/DC.ndgd/GT.rtma/AR.conus/'
        rtma_base = 'http://tgftp.nws.noaa.gov/SL.us008001/ST.opnl/DF.gr2/DC.ndgd/GT.rtma/AR.conus/'
        return rtma_base + '/RT.%02d/' % hour + 'ds.%s.bin' % var


    @staticmethod
    def _is_var_cached(path):
        """
        Check if the file for the variable is available in the ingest directory
        already and it's filesize is correct.

        :param path: the (absolute) path to the variable file
        :param var: the variable name
        :return: True if the file is available and verified (via file size), False otherwise
        """
        info_path = path + '.size'
        if osp.exists(path) and osp.exists(info_path):
            content_size = int(open(info_path).read())
            return osp.getsize(path) == content_size
        else:
            return False


    @staticmethod
    def _parse_header_timestamp(ts):
        """
        Parse a timestamp in the header, example 'Tue, 12 Apr 2016 18:51:18 GMT'.

        :param ts: the timestamp
        :return: a datetime object in UTC timezone
        """
        return datetime.strptime(ts, '%a, %d %b %Y %H:%M:%S GMT').replace(tzinfo=pytz.UTC)


if __name__ == '__main__':
    import sys

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if len(sys.argv) != 2:
        print('usage: %s <esmf_time>' % sys.argv[0])
        sys.exit(1)

    # initialize the RTMA object with standard variables
    rtma = RTMA('ingest', ['utd', 'utemp', 'precipa'])

    # try to download them
    rtma.retrieve_rtma(esmf_to_utc(sys.argv[1]))

