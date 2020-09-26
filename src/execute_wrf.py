from forecast import wrf_execute, process_output
import logging, sys

if __name__ == '__main__':

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # get job_id
    job_id = sys.argv[1]    

    # execute the job
    logging.info('calling wrf_execute')
    wrf_execute(job_id)
    logging.info('calling process_output')
    process_output(job_id)

    logging.info('execute_wrf.py done')
