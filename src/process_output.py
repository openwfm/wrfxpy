from forecast import process_output
import sys, logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
process_output(sys.argv[1])
