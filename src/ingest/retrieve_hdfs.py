#
# Dalton Burke, UC Denver
#

from ingest.level0_source import MODIS_AQUA, MODIS_TERRA, VIIRS_NPP
from utils import esmf_to_utc

import logging
import sys
import os.path as osp


## Standalone script that can be used to download files
if __name__ == '__main__':
    if len(sys.argv) != 9:
        print('Usage: %s <hdf_source_name> <esmf_from_utc> <esmf_to_utc> <low_long> <high_long> <low_lat> <high_lat> <target_directory>' % sys.argv[0])
        print('        supported HDF sources: MODIS_AQUA, MODIS_TERRA, VIIRS_NPP')
        print('        time of format YYYY.MM.DD-HH.MM.SS')
        sys.exit(-1)


    # configure basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    data_src_name = sys.argv[1]
    from_utc = esmf_to_utc(sys.argv[2])
    to_utc = esmf_to_utc(sys.argv[3])
    lonlat = [float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])]
    ingest_dir = osp.abspath(sys.argv[8])

    hdf_src = None
    if data_src_name == 'MODIS_AQUA':
        data_src = MODIS_AQUA(ingest_dir)
    if data_src_name == 'MODIS_TERRA':
        data_src = MODIS_TERRA(ingest_dir)
    elif data_src_name == 'VIIRS_NPP':
        data_src = VIIRS_NPP(ingest_dir)
    else:
        raise ValueError('Invalid HDF source %s' % data_src_name)

    logging.info('Initiating download of files from HDF source %s' % grib_src_name)

    hdfs = data_src.retrieve_data(from_utc, to_utc, lonlat)

    logging.info('SUCCESS, the following files are now available:')
    print('')

    for h in hdfs:
        print(osp.join(ingest_dir, h))
