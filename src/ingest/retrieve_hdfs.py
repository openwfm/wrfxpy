#
# Dalton Burke, UC Denver
#

from ingest.hdf_source import MODIS, VIIRS
from utils import esmf_to_utc

import logging
import sys
import os.path as osp


## Standalone script that can be used to download files
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: %s <hdf_source_name> <esmf_from_utc> <esmf_to_utc> <target_directory>' % sys.argv[0])
        print('        supported HDF sources: MODIS, VIIRS')
        sys.exit(-1)


    # configure basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    hdf_src_name = sys.argv[1]
    from_utc = esmf_to_utc(sys.argv[2])
    to_utc = esmf_to_utc(sys.argv[3])
    ingest_dir = sys.argv[4]

    hdf_src = None
    if grib_src_name == 'MODIS':
        hdf_src = MODIS(ingest_dir)
    elif hdf_src_name == 'VIIRS':
        hdf_src = VIIRS(ingest_dir)
    else:
        raise ValueError('Invalid HDF source %s' % grib_src_name)

    logging.info('Initiating download of files from HDF source %s' % grib_src_name)

    hdfs = hdf_src.retrieve_hdfs(from_utc, to_utc)

    logging.info('SUCCESS, the following files are now available:')
    print('')

    for h in hdfs:
        print(osp.join(ingest_dir, h))

    print('\n** NOTE **')
    print('The following variable tables must be used with this hdf source:')
    print(repr(hdf_src.vtables()))
    print('The following keys must be set in the namelists:')
    print(repr(hdf_src.namelist_keys()))


