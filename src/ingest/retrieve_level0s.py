#
# Dalton Burke
#

from ingest.level0_source import MODIS_TERRA, MODIS_AQUA, VIIRS_NPP
from utils import esmf_to_utc

import logging
import sys
import os.path as osp

## Standalone script to simply download files
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print(('Usage: %s <level0_source_name> <esmf_from_utc> <esmf_to_utc> <target_directory>' % sys.argv[0]))
        print('\tsupported level0 sources: MODIS_TERRA, MODIS_AQUA, VIIRS_NPP')
        sys.exit(-1)

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    level0_src_name = sys.argv[1]
    from_utc = esmf_to_utc(sys.argv[2])
    to_utc = esmf_to_utc(sys.argv[3])
    ingest_dir = sys.argv[4]

    level0_src = None
    if level0_src_name == 'MODIS_TERRA':
        level0_src = MODIS_TERRA(ingest_dir)
    elif level0_src_name == 'MODIS_AQUA':
        level0_src = MODIS_AQUA(ingest_dir)
    elif level0_src_name == 'VIIRS_NPP':
        level0_src = VIIRS_NPP(ingest_dir)
    else:
        raise ValueError('Invalid level0 source %s' % grib_src_name)

    logging.info('Initiating download of files from level0 source %s' % level0_src_name)

    level0s = level0_src.retrieve_level0s(from_utc, to_utc)

    logging.info('SUCCESS, the following files are now available:')
    print('')
    for l in level0s:
        print((osp.join(ingest_dir, l)))
