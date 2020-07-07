#
# Dalton Burke, CU Denver
#

# from ingest.level0_source import MODIS_AQUA, MODIS_TERRA, VIIRS_NPP
from __future__ import absolute_import
from __future__ import print_function
from . import var_wisdom as vw
from utils import esmf_to_utc

import logging
import sys
import os.path as osp



if __name__ == '__main__':
    if len(sys.argv) != 9  and len(sys.argv) != 5:
        print(('Usage: %s <sat_name> <esmf_from_utc> <esmf_to_utc> <target_directory> <low_long> <high_long> <low_lat> <high_lat>' % sys.argv[0]))
        print('        supported HDF sources: MODIS_AQUA, MODIS_TERRA')
        print('        time of format YYYY.MM.DD-HH.MM.SS')
        print('        (lonlat parameters optional, CONUS = -124.7844079 -66.9513812 24.7433195 49.3457868)')
        print('example: ./retrieve_hdfs.sh MODIS_AQUA 2017.06.20-10.00.00 2017.06.20-15.00.00 ~/hdf_test -124.7844079 -66.9513812 24.7433195 49.3457868')
        sys.exit(-1)

    # configure basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    sat_name = sys.argv[1]
    if sat_name not in vw.get_sat_wisdom_variables():
        print('Invalid instrument/satellite name!')
        print('supported instrument/satellite pairs:')
        for sat in vw.get_sat_wisdom_variables():
            print(sat)
        sys.exit(-1)

    from_utc = esmf_to_utc(sys.argv[2])
    from_utc = from_utc.replace(tzinfo=None)
    to_utc = esmf_to_utc(sys.argv[3])
    to_utc = to_utc.replace(tzinfo=None)

    ingest_dir = osp.abspath(osp.expanduser(sys.argv[4]))

    lonlat = []
    if len(sys.argv) == 9:
        lonlat = [float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])]

    logging.info('sat_name: %s' % sat_name)
    logging.info('from_utc:      %s' % from_utc)
    logging.info('to_utc:        %s' % to_utc)
    logging.info('lonlat:        %s:%s %s:%s' % tuple(lonlat))
    logging.info('ingest_dir:    %s' % ingest_dir)

### Standalone script that can be used to download files
#if __name__ == '__main__':
#    if len(sys.argv) != 9:
#        print('Usage: %s <data_src_name> <esmf_from_utc> <esmf_to_utc> <low_long> <high_long> <low_lat> <high_lat> <target_directory>' % sys.argv[0])
#        print('        supported HDF sources: MODIS_AQUA, MODIS_TERRA, VIIRS_NPP')
#        print('        time of format YYYY.MM.DD-HH.MM.SS')
#        sys.exit(-1)
#
#
#    # configure basic logger
#    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
#
#    data_src_name = sys.argv[1]
#    from_utc = esmf_to_utc(sys.argv[2])
#    from_utc = from_utc.replace(tzinfo=None)
#    to_utc = esmf_to_utc(sys.argv[3])
#    to_utc = to_utc.replace(tzinfo=None)
#    lonlat = [float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])]
#    ingest_dir = osp.abspath(osp.expanduser(sys.argv[8]))
#
#    logging.info('data_src_name: %s' % data_src_name)
#    logging.info('from_utc:      %s' % from_utc)
#    logging.info('to_utc:        %s' % to_utc)
#    logging.info('lonlat:        %s:%s %s:%s' % tuple(lonlat))
#    logging.info('ingest_dir:    %s' % ingest_dir)
#
#    if data_src_name == 'MODIS_AQUA':
#        data_src = MODIS_AQUA(ingest_dir)
#    elif data_src_name == 'MODIS_TERRA':
#        data_src = MODIS_TERRA(ingest_dir)
#    elif data_src_name == 'VIIRS_NPP':
#        data_src = VIIRS_NPP(ingest_dir)
#    else:
#        raise ValueError('Invalid HDF source %s' % data_src_name)

    logging.info('Initiating download of files from HDF source %s' % sat_name)

    hdfs = vw.retr_sat(ingest_dir, sat_name, from_utc, to_utc, lonlat)

    logging.info('SUCCESS, the following files are now available:')
    print('')

    for h in hdfs:
        print((osp.join(ingest_dir, h)))
