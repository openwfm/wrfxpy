from forecast import restart_execute
import logging, sys

if __name__ == '__main__':

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # get job_id
    job_id = sys.argv[1]    

    # execute the job
    logging.info('calling restart_execute')
    restart_execute(job_id)

    logging.info('restart_forecast.py done')
