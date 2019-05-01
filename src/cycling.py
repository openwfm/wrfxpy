from forecast import execute, process_arguments
from utils import load_sys_cfg, matching_files, make_dir, copy, move
import sys, json, logging


def cycling_verify_inputs(args):
	required_inputs = [('ncycles', 'Non-existent number of cycles ncycles in %s'),
			('dcycle', 'Non-existent length of cycle in days dcycle in %s')]
	for key, err in required_inputs:
		if key not in args.keys():
			raise OSError(err % args.keys()) 	


def cycling_step(args,jobfile,cycle):

	logging.info('cycling_step: loading job description from %s' % jobfile)
        try:
                js = Dict(json.load(open(jobfile,'r')))
        except Exception as e:
                logging.error('Cannot load the job description file %s' % jobfile)
                logging.error('%s' % e)
                sys.exit(1)
	args['restart'] = True	

	domains = args['domains'].keys()
	cycle_dir = osp.join(js.wrf_dir,'cycle%01d' % cycle-1)	
	make_dir(cycle_dir)
	reg_cycle = ('wrfrst*', 'wrfout*', 'rsl.error.*', 'rsl.out.*')
	for reg in reg_cycle:
		for f in matching_files(js.wrf_dir,reg):
			move(f,cycle_dir)

	for domain in domains:
		dom_id = int(domain)
		rst_dir = sorted(glob.glob(osp.join(cycle_dir,'wrfrst_d%02d*' % dom_id)))[-1]  	
		copy(rst_dir,js.wrf_dir)
		esmf_time = re.match(r'.*wrfrst_d%02d_([0-9_\-:]{19})' % dom_id, rst_dir).groups()[0]
		logging.info('domain %s, at time %s' % (domain, esmf_time))


if __name__ == '__main__':
	# configure the basic logger
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        
	# load configuration JSON
        sys_cfg = load_sys_cfg()

        # load job JSON
        job_args = json.load(open(sys.argv[1]), 'ascii')
	
        # process arguments
        args = process_arguments(job_args,sys_cfg)
	
	# process cycling arguments
	cycling_verify_inputs(args)
		
	logging.info('first simulation for cycling')
        # execute the job
        logging.info('calling execute')
        jobfile = execute(args,job_args)
        logging.info('execute done')
	logging.info('returned job file path %s' % jobfile)
	
	for cycle in range(1,args['ncycles']+1):
		logging.info('cycling step %d of %d' % (cycle,args['ncycles']))
		cycling_step(args,job_file,cycle)
		logging.info('cycling step %d done' % cycle)
        
	logging.info('cycling.py done')
        
        
