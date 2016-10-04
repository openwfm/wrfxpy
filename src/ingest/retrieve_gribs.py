#!/usr/bin/env python
# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
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


from ingest.grib_source import HRRR, NAM218, NAM227, NARR
from utils import esmf_to_utc

import logging
import sys
import os.path as osp

## Standalone script that can be used to simply download files
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: %s <grib_source_name> <esmf_from_utc> <esmf_to_utc> <target_directory>' % sys.argv[0])
        print('       supported GRIB sources: HRRR, NAM, NARR')
        sys.exit(-1)

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    grib_src_name = sys.argv[1]
    from_utc = esmf_to_utc(sys.argv[2])
    to_utc = esmf_to_utc(sys.argv[3])
    ingest_dir = sys.argv[4]

    grib_src = None
    if grib_src_name == 'HRRR':
        grib_src = HRRR(ingest_dir)
    elif grib_src_name == 'NAM':
        grib_src = NAM218(ingest_dir)
    elif grib_src_name == 'NAM227':
        grib_src = NAM227(ingest_dir)
    elif grib_src_name == 'NARR':
        grib_src = NARR(ingest_dir)
    else:
        raise ValueError('Invalid GRIB source %s' % grib_src_name)

    logging.info('Initiating download of files from GRIB source %s' % grib_src_name)

    gribs = grib_src.retrieve_gribs(from_utc, to_utc)
    
    logging.info('SUCCESS, the following files are now available:')
    print('')
    for g in gribs:
        print(osp.join(ingest_dir, g))
    
    print('\n** NOTE **')
    print('The following variable tables must be used with this grib source:')
    print(repr(grib_src.vtables()))
    print('The following keys must be set in namelists:')
    print(repr(grib_src.namelist_keys()))


