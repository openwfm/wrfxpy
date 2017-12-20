# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
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

from ssh_shuttle import SSHShuttle
import json
import os
import os.path as osp
import sys
import logging
import shutil
import glob
import signal
import subprocess
from utils import kill_process, Dict, process_create_time, load_sys_cfg

cfg = load_sys_cfg()

def remote_rmdir(s,dirname):
    logging.info('SHUTTLE removing remote directory %s' % dirname)
    try:
        s.rmdir(dirname)
        logging.info('SHUTTLE delete successful.')
        return 0
    except:
        logging.error('SHUTTLE cannot delete directory %s' % dirname)
        return 'Error'

def local_rmdir(dirname):
    work_dir = osp.abspath(osp.join(cfg['workspace_path'], dirname))
    logging.info('Deleting directory %s' % work_dir)
    try:
        shutil.rmtree(work_dir)
        logging.info('Deleted directory %s' % work_dir)
        return 0
    except:
        logging.error('Deleting directory %s failed' % work_dir)
        return 'Error'

def cancel_job(job_num,qsys):
        if job_num is not None:
            logging.info('Deleting parallel job %s on %s.' % (job_num, qsys))
            cluster = load_cluster_file(qsys)
            try:
                ret = subprocess.check_output([cluster['qdel_cmd'], job_num])
                logging.info(ret)
            except:
                logging.error('Deleting parallel job %s failed.' % job_num)

def load_cluster_file(qsys):
    clusters = json.load(open('etc/clusters.json','r'))
    if qsys in clusters:
        return clusters[qsys]
    else:
        logging.warning('Cluster %s not known' % qsys)
        return None

def load_job_file(job_id):
    jobfile = osp.join(cfg['workspace_path'], job_id, 'job.json')
    logging.info('Loading job state from %s' % jobfile)
    try:
        js = Dict(json.load(open(jobfile,'r')))
    except IOError:
        logging.error('Cannot open %s' % jobfile)
        js = None
    except:
        logging.error('Cannot load %s' % jobfile)
        js = None
    # print(json.dumps(js, indent=4, separators=(',', ': ')))
    return js, jobfile

def forecast_process_running(js):
    if 'pid' in js and 'process_create_time' in js:
        fs =  js.process_create_time == process_create_time(js.pid)
        if fs:
           logging.info('Forecast process is running. pid=%s' % (js.pid))
        else:
            logging.info('Forecast process is not running')
        return fs
    else:
        logging.warning('Cannot determine if forecast process is running - old job file?')
        return None 

def parallel_job_running(js):
    logging.debug('Checking for WRF job %s' % js.job_num)
    cluster = load_cluster_file(js.qsys)
    if cluster is None:
        logging.warning('No cluster system, cannot check for WRF job')
        return False
    ret = subprocess.check_output(cluster['qstat_cmd'],stderr=subprocess.STDOUT)
    for line in ret.split('\n'):
        ls=line.split()
        if len(ls) >0 and ls[0] == js.job_num:
             status = ls[4]
             logging.debug('WRF job %s status is %s' % (js.job_num, status))
             return True
    logging.info('WRF job %s is not running.' % js.job_num)
    return False 

# the cleanup functions called as argv[2]:
# connecting to remote host if needed

def list(s):
    s.connect()
    cat = s.retrieve_catalog()
    s.disconnect()
    print('%-60s desc' % 'id')
    print('-' * 70)
    for k in sorted(cat):
        print('%-60s %s' % (k, cat[k]['description']))

def workspace(s):
    logging.info('Deleting all directories in local workspace that are not in the remote catalog.')
    s.connect()
    cat = s.retrieve_catalog(s)
    s.disconnect()
    for f in glob.glob(osp.join(cfg['workspace_path'],'*')):
        if osp.isdir(f):
            ff = osp.basename(f)
            if ff not in cat:
                logging.error('%s not in the catalog' % ff)
                local_rmdir(cfg, f)

def output(s,name):
    logging.info('Trying to delete WRF output and visualization of job %s' % name)
    s.connect()
    remote_rmdir(s, name)
    s.disconnect()
    local_rmdir(cfg,osp.join(name,'products'))
    wrf_dir = osp.join(cfg['workspace_path'], name,'wrf')
    files = glob.glob(osp.join(wrf_dir,'rsl.*'))+glob.glob(osp.join(wrf_dir,'wrfout_*'))
    logging.info('Deleting %s WRF output files in  %s' % (len(files),wrf_dir ))
    for f in files:
        os.remove(f)

def cancel(name):
    logging.info('Trying to cancel job %s' % name)
    js, jobfile = load_job_file(name)
    if js is not None:
        if 'state' in js and js.state == 'Cancelled':
            logging.warning('Job %s was already cancelled, proceeding anyway' % name)
        if 'job_num' in js and js.job_num is not None:
            try:
                logging.info('Trying to delete job %s from the queue scheduler' % js.job_num)
                cancel_job(js.job_num,js.qsys)
            except:
                logging.error('Delete failed, please check your queue scheduler for job %s' % js.job_num)
            js.old_job_num = js.job_num
        js.job_num = None
        if 'pid' in js and js.pid is not None:
            logging.info('Trying to kill forecast process %s' % js.pid)
            kill_process(js.pid)
        js.pid = None
        js.state = "Cancelled"
        json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

def delete(s,name):
    s.connect()
    logging.info('Trying to delete all files of job %s' % name)
    remote_rmdir(s,name)
    local_rmdir(name)
    s.simple_command('wrfxweb/join_catalog.sh')
    s.disconnect()

def update(name):
    logging.info('Updating state of job %s' % name)
    js, jobfile = load_job_file(name)
    if js is not None:
        logging.debug('Old state is: %s' % js.state)
        if not forecast_process_running(js):
            logging.debug('Forecast process is not running')
            js.pid = None
        if not parallel_job_running(js):    
            logging.debug('WRF is not running')
            js.old_job_num = js.job_num
            js.job_num = None
        if js.state == 'Completed' or js.state == 'Cancelled':
            if js.pid is not None:
                 js.state = 'Forecast runaway'
            if js.job_num is not None:
                 js.state = 'WRF runaway'
        elif js.state == 'Preparing':
                 if js.pid is None:
                      js.state = 'Stopped'
                 if js.job_num is not None:
                      js.state = 'Cannot happen'
        elif js.state == 'Processing':
                 if js.pid is None and js.job_num is None:
                      js.state = 'Stopped'
                 if js.pid is not None and js.job_num is None:
                      js.state = 'WRF stopped'
                 if js.pid is None and js.job_num is not None:
                      js.state = 'Forecast stopped'

        logging.info('State is: %s' % js.state)

        json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    commands = [ 'list', 'cancel', 'output', 'delete', 'workspace', 'update']

    if len(sys.argv) < 2 or sys.argv[1] not in commands: 
        print('usage: ./cleanup.sh ' + '|'.join(commands) +' [job_id]')
        print('list            : show list of current simulations with their job_id and description')
        print('cancel <job_id> : kill all processes and the WRF parallel job, do not delete any files')
        print('output <job id> : cancel, and delete all WRF output and visualization files only')
        print('delete <job_id> : cancel, and delete all files')
        print('workspace       : delete jobs that are not on the visulalization server')
        print('update <job_id> : check if the job is running and update its job state file')
        sys.exit(1)

    cmd = sys.argv[1]
    if cmd in ['delete' , 'cancel', 'output', 'update']: 
        if len(sys.argv) == 3 and not sys.argv[2] == "" :
            job_id = sys.argv[2]
        else:
            print('%s function requires one job id' % cmd)
            sys.exit(1)
    else:
        job_id = None

    logging.info('cleanup: command=%s job_id=%s' % (cmd,job_id))
    s = SSHShuttle(cfg)

    
    if cmd == 'list':
        list(s)

    if cmd == 'cancel':
        cancel(job_id)

    if cmd == 'output':
        cancel(job_id)
        output(s,job_id)

    if cmd == 'delete':
        cancel(job_id)
        delete(s,job_id)

    if cmd == 'workspace':
        workspace(s)

    if cmd == 'update':
        update(job_id)

