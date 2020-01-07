from __future__ import absolute_import
from forecast import make_job_file, JobState, process_arguments, load_sys_cfg
import json
import logging
import sys

if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    # load configuration JSON
    sys_cfg = load_sys_cfg()
    # note: the execution flow allows us to override anything in the etc/conf.json file
    job_args = json.load(open(sys.argv[1]), 'ascii')
    args = sys_cfg 
    args.update(job_args)
    process_arguments(args)
    js = JobState(args)
    jsub = make_job_file(js)
    jsub.pid = None
    jsub.process_create_time = None
    logging.debug('job.json:\n%s' % json.dumps(jsub, indent=4, separators=(',', ': ')))
    try:
        jsub_old=json.load(open(jsub.jobfile,'r'))
        logging.info('Found existing job.json file, updating.')
        logging.debug('%s' % json.dumps(jsub_old, indent=4, separators=(',', ': ')))
        jsub.update(jsub_old)
        logging.debug('updated:\n%s' % json.dumps(jsub, indent=4, separators=(',', ': ')))
    except:
        logging.info('No existing job.json file found.')

    #json.dump(jsub, open(jsub.jobfile,'w'), indent=4, separators=(',', ': '))
    # logging.info('Written job state to file %s' % sub.jobfile)
