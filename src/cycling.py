from forecast import execute, process_arguments
from utils import load_sys_cfg
import sys, json, logging

if __name__ == '__main__':
	# configure the basic logger
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        # logging.basicConfig(level=logging.DEBUG)

        # load configuration JSON
        sys_cfg = load_sys_cfg()
        # logging.info('sys_cfg = %s' % json.dumps(sys_cfg, indent=4, separators=(',', ': ')))

        # load job JSON
        job_args = json.load(open(sys.argv[1]), 'ascii')
        # logging.info('job_args = %s' % json.dumps(job_args, indent=4, separators=(',', ': ')))

        # process arguments
        args = process_arguments(job_args,sys_cfg)
        # logging.info('processed args = %s' % str(args))

        # execute the job
        logging.info('calling execute')
        execute(args,job_args)

        logging.info('forecast.py done')
