#
# Dalton Burke, CU Denver
#

#from utils import ensure_dir, symlink_unless_exists
from .downloader import download_url, DownloadError, get_dList

from datetime import datetime, timedelta

import os
import os.path as osp
import sys
import logging


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Downloads the latest firms MODIS and VIIRS data from US and Alaska zones')
        print('Usage: %s <target_directory>' % sys.argv[0])
        sys.exit(-1)

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    ingest_dir = sys.argv[1]

    # today's data is not available yet, we want yesterdays which has recently become available
    yesterday = datetime.utcnow() - timedelta(days=1)
    julian_day = (yesterday - datetime(yesterday.year, 1, 1)).days + 1
    year = yesterday.year

    urls = ['ftp://dmtburke:WRFfire6@nrt3.modaps.eosdis.nasa.gov/FIRMS/c6/USA_contiguous_and_Hawaii/',
            'ftp://dmtburke:WRFfire6@nrt3.modaps.eosdis.nasa.gov/FIRMS/c6/Alaska/',
            'ftp://dmtburke:WRFfire6@nrt3.modaps.eosdis.nasa.gov/FIRMS/viirs/USA_contiguous_and_Hawaii/',
            'ftp://dmtburke:WRFfire6@nrt3.modaps.eosdis.nasa.gov/FIRMS/viirs/Alaska/']

    filenames = ['MODIS_C6_USA_contiguous_and_Hawaii_MCD14DL_NRT_%04d%03d.txt' % (year, julian_day),
                 'MODIS_C6_Alaska_MCD14DL_NRT_%04d%03d.txt' % (year, julian_day),
                 'VIIRS_I_USA_contiguous_and_Hawaii_VNP14IMGTDL_NRT_%04d%03d.txt' % (year, julian_day),
                 'VIIRS_I_Alaska_VNP14IMGTDL_NRT_%04d%03d.txt' % (year, julian_day)]

    for i in range(len(urls)):
        download_url(urls[i] + filenames[i], ingest_dir + '/' + filenames[i])

    logging.info('SUCCESS, the following files are now available:')

    print('')
    for f in filenames:
        print((osp.join(ingest_dir, f)))
