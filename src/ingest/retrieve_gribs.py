#!/usr/bin/env python
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


from __future__ import absolute_import
from __future__ import print_function
from ingest.HRRR import HRRR
from ingest.HRRR_AK import HRRR_AK
from ingest.NARR import NARR
from ingest.NAM218 import NAM218
from ingest.NAM227 import NAM227
from ingest.NAM198 import NAM198
from ingest.NAM196 import NAM196
from ingest.CFSR import CFSR_P, CFSR_S
from ingest.GFSA import GFSA
from ingest.GFSF import GFSF_P, GFSF_S
from utils import esmf_to_utc, load_sys_cfg

import logging
import sys
import os.path as osp

## Standalone script that can be used to simply download files
if __name__ == '__main__':
    if len(sys.argv) not in (5, 6):
        print(('Usage: %s <grib_source_name> <esmf_from_utc> <esmf_to_utc> <target_directory> [esmf_cycle_start]' % sys.argv[0]))
        print('       supported GRIB sources: HRRR, HRRR_AK, NAM, NAM227, NAM196, NAM198, CFSR_P, CFSR_S, NARR, GFSA, GFSF_P, GFSF_S')
        print('       cycle_start pins which forecast cycle is used instead of letting')
        print('       the source pick the most recent one.  Needed for HRRR, whose')
        print('       cycle_hours is 1: without it a caller asking for the 00/06/12/18')
        print('       cycles silently gets whatever hourly cycle happens to be newest.')
        sys.exit(-1)


    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    js = load_sys_cfg()
    grib_src_name = sys.argv[1]
    from_utc = esmf_to_utc(sys.argv[2])
    to_utc = esmf_to_utc(sys.argv[3])
    ingest_dir = sys.argv[4]
    #optional: pin the cycle rather than taking the source's most recent.  The
    #parameter has always existed on retrieve_gribs -- hrrr_cycler passes it -- but was
    #not reachable from the command line, which is how cache_grib_files drives this.
    cycle_start = esmf_to_utc(sys.argv[5]) if len(sys.argv) == 6 else None
    js.ingest_dir = ingest_dir

    grib_src = None
    #'elif', not 'if'.  With a plain 'if' the HRRR branch set grib_src and then fell
    #into the HRRR_AK/NAM/... chain below, matched nothing, and hit its else -- so
    #'retrieve_gribs.sh HRRR' always died with "Invalid GRIB source HRRR" even though
    #the source is fully supported.  HRRR_AK was unaffected because it is the chain's
    #own first test.
    if grib_src_name == 'HRRR':
        grib_src = HRRR(js)
    elif grib_src_name == 'HRRR_AK':
        grib_src = HRRR_AK(js)
    elif grib_src_name == 'NAM':
        grib_src = NAM218(js)
    elif grib_src_name == 'NAM227':
        grib_src = NAM227(js)
    elif grib_src_name == 'NAM198':
        grib_src = NAM198(js)
    elif grib_src_name == 'NAM196':
        grib_src = NAM196(js)
    elif grib_src_name == 'CFSR_P':
        grib_src = CFSR_P(js)
    elif grib_src_name == 'CFSR_S':
        grib_src = CFSR_S(js)
    elif grib_src_name == 'NARR':
        grib_src = NARR(js)
    elif grib_src_name == 'GFSA':
        grib_src = GFSA(js)
    elif grib_src_name == 'GFSF_P':
        grib_src = GFSF_P(js)
    elif grib_src_name == 'GFSF_S':
        grib_src = GFSF_S(js)
    else:
        raise ValueError('Invalid GRIB source %s' % grib_src_name)

    logging.info('Initiating download of files from GRIB source %s' % grib_src_name)

    if cycle_start is not None:
        logging.info('pinning forecast cycle to %s' % cycle_start)
        gribs = grib_src.retrieve_gribs(from_utc, to_utc, cycle_start=cycle_start)
    else:
        gribs = grib_src.retrieve_gribs(from_utc, to_utc)

    logging.info('SUCCESS, the following files are now available:')
    print('')
    #retrieve_gribs returns a manifest dict, not a list, so iterating it printed the
    #dict's keys -- 'ingest/grib_files', 'ingest/colmet_prefix' and so on, which read
    #exactly like paths and made an empty download look like a successful one.
    for g in (gribs['grib_files'] if isinstance(gribs, dict) else gribs):
        print((osp.join(ingest_dir, g)))

    print('\n** NOTE **')
    print('The following variable tables must be used with this grib source:')
    print((repr(grib_src.vtables())))
    print('The following keys must be set in namelists:')
    print((repr(grib_src.namelist_keys())))


