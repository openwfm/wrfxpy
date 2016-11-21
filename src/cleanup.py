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
from utils import kill_process, Dict, process_create_time


def retrieve_catalog(cfg):
    """
    Retrieve the catalog from the remote visualization server.
    
    :param cfg: the configuration of the remote visualization host
    """
    s = SSHShuttle(cfg)
    logging.info('SHUTTLE connecting to remote host %s' % s.host)
    s.connect()
    logging.info('SHUTTLE retrieving catalog file.')
    s.get('catalog.json', '_catalog.json')

    cat = json.load(open('_catalog.json'))
    s.disconnect()
    os.remove('_catalog.json')
    logging.info('SHUTTLE retrieve complete.')
    return cat


def store_catalog(cfg, cat):
    """
    Store a new version of the catalog on the remote server given
    by cfg.

    :param cfg: the configuration of the remote visualization host
    :param cat: the JSON object representing the new catalog
    """
    s = SSHShuttle(cfg)
    logging.info('SHUTTLE connecting to remote host %s' % s.host)
    s.connect()
    logging.info('SHUTTLE sending catalog file.')
    json.dump(cat, open('_catalog.json', 'w'), indent=4, separators=(',', ': '))
    s.put('_catalog.json', 'catalog.json')
    os.remove('_catalog.json')
    s.disconnect()
    logging.info('SHUTTLE send complete.')

def remote_rmdir(cfg, dirname):
    s = SSHShuttle(cfg)
    logging.info('SHUTTLE connecting to remote host %s' % s.host)
    s.connect()
    logging.info('SHUTTLE removing remote directory %s' % dirname)
    try:
        s.rmdir(dirname)
        logging.info('SHUTTLE delete successful.')
        s.disconnect()
        return 0
    except:
        logging.error('SHUTTLE cannot delete directory %s' % dirname)
        s.disconnect()
        return 'Error'

def local_rmdir(cfg, dirname):
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
    return clusters[qsys]

def load_job_file(job_id,cfg):
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
        logging.info('Forecast process is running: %s pid=%s' % (fs, js.pid))
        return fs
    else:
        logging.warning('Cannot determine if forecast process is running - old job file?')
        return None 

def parallel_job_running(js):
    cluster = load_cluster_file(js.qsys)
    ret = subprocess.check_output(cluster['qstat_cmd'])
    ws = False
    for line in ret.split('\n'):
        ls=line.split()
        if len(ls) >0 and ls[0] == js.job_num:
             ws = True
             break
    logging.info('WRF job is running: %s name=%s' % (ws, js.job_num))
    return ws

# the cleanup functions called as argv[2]:

def list(cfg):
        cat = retrieve_catalog(cfg)
        print('%-60s desc' % 'id')
        print('-' * 70)
        for k in sorted(cat):
            print('%-60s %s' % (k, cat[k]['description']))

def workspace(cfg):
        logging.info('Deleting all directories in local workspace that are not in the remote catalog.')
        cat = retrieve_catalog(cfg)
        for f in glob.glob(osp.join(cfg['workspace_path'],'*')):
            if osp.isdir(f):
                ff = osp.basename(f)
                if ff not in cat:
                    logging.error('%s not in the catalog' % ff)
                    local_rmdir(cfg, f)

def output(name,cfg):
    logging.info('Trying to delete WRF output and visualization of job %s' % name)
    remote_rmdir(cfg, name)
    local_rmdir(cfg,osp.join(name,'products'))
    wrf_dir = osp.join(cfg['workspace_path'], name,'wrf')
    files = glob.glob(osp.join(wrf_dir,'rsl.*'))+glob.glob(osp.join(wrf_dir,'wrfout_*'))
    logging.info('Deleting %s WRF output files in  %s' % (len(files),wrf_dir ))
    cat = retrieve_catalog(cfg)
    for f in files:
        os.remove(f)
    if name not in cat:
        logging.error('Simulation %s not in the catalog' % name)
    else:
        logging.info('Deleting simulation %s from the catalog' % name)
        del cat[name]
    store_catalog(cfg, cat)

def cancel(name,cfg):
    logging.info('Trying to cancel job %s' % name)
    js, jobfile = load_job_file(name,cfg)
    if 'state' in js and js.state == 'Cancelled':
        logging.warning('Job %s was already cancelled, proceeding anyway' % name)
    if js is not None:
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

def delete(name,cfg):
        logging.info('Trying to delete all files of job %s' % name)
        remote_rmdir(cfg, name)
        local_rmdir(cfg,name)
        cat = retrieve_catalog(cfg)
        if name not in cat:
            logging.error('Simulation %s not in the catalog' % name)
        else:
            logging.info('Deleting simulation %s from the catalog' % name)
            del cat[name]
            store_catalog(cfg, cat)

def update(name,cfg):
    logging.info('Updating state of job %s' % name)
    js, jobfile = load_job_file(name,cfg)
    if js is not None:
        logging.info('Old state is: %s' % js.state)
        if not forecast_process_running(js):
            #logging.info('Forecast process is not running')
            js.pid = None
        if not parallel_job_running(js):    
            #logging.info('WRF is not running')
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

        logging.info('New state is: %s' % js.state)


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

    cfg = Dict(json.load(open('etc/conf.json')))

    if cmd == 'list':
        list(cfg)

    if cmd == 'cancel':
        cancel(job_id,cfg)

    if cmd == 'output':
        cancel(job_id,cfg)
        output(job_id,cfg)

    if cmd == 'delete':
        cancel(job_id,cfg)
        delete(job_id,cfg)

    if cmd == 'workspace':
        workspace(cfg)

    if cmd == 'update':
        update(job_id,cfg)
