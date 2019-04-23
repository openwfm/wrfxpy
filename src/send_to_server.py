from forecast import send_products_to_server
import sys, logging

if len(sys.argv) < 2: 
    raise SystemExit('usage: ./send_to_server.sh job_id')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
send_products_to_server(sys.argv[1])
