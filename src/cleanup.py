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
from utils import kill_process, Dict


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
            clusters = json.load(open('etc/clusters.json','r'))
            cluster = clusters[qsys]
            try:
                ret = subprocess.check_output([cluster['qdel_cmd'], job_num])
                logging.info(ret)
            except:
                logging.error('Deleting parallel job %s failed.' % job_num)

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
    jobfile = osp.join(cfg['workspace_path'], name, 'job.json')
    logging.info('Loading job state from %s' % jobfile)
    try:
        js = Dict(json.load(open(jobfile,'r')))
    except IOError:
        logging.error('Cannot open %s' % jobfile)
    except:
        logging.error('Cannot load %s' % jobfile)
    else:
        print(json.dumps(js, indent=4, separators=(',', ': ')))
        if 'job_num' in js:
            cancel_job(js.job_num,js.qsys)
            js.old_job_num = js.job_num
        js.job_num = None
        if 'pid' in js:
            kill_process(js.pid)
            js.old_pid = js.pid
        js.pid = None
        js.status = "Cancelled"
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

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    commands = [ 'list', 'cancel', 'output', 'delete', 'workspace']

    if len(sys.argv) < 2 or sys.argv[1] not in commands: 
        print('usage: ./cleanup.sh ' + '|'.join(commands) +' [job_id]')
        print('list            : show list of current simulations with their job_id and description')
        print('cancel <job_id> : kill all processes and the WRF parallel job, do not delete any files')
        print('output <job id> : cancel, and delete all WRF output and visualization files only')
        print('delete <job_id> : cancel, and delete all files')
        print('workspace       : delete jobs that are not on the visulalization server')
        sys.exit(1)

    cmd = sys.argv[1]
    if cmd in ['delete' , 'cancel', 'output']: 
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
