from forecast import process_output
import sys, logging

if len(sys.argv) < 2: 
    raise SystemExit('usage: ./process_output.sh job_id')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
process_output(sys.argv[1])
