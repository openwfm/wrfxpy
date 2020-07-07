from __future__ import absolute_import
from forecast import process_sat_output
import sys, logging

if len(sys.argv) < 2: 
    raise SystemExit('usage: ./process_sat_output.sh job_id')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
process_sat_output(sys.argv[1])
