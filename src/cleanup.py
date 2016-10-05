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
    logging.info('SHUTTLE removing directory %s' % dirname)
    try:
        s.rmdir(dirname)
        logging.info('SHUTTLE delete successful.')
        s.disconnect()
        return 0
    except:
        logging.info('SHUTTLE cannot delete directory %s' % dirname)
        s.disconnect()
        return 'Error'

def local_rmdir(cfg, dirname):
    work_dir = osp.abspath(osp.join(cfg['workspace_path'], dirname))
    logging.info('Deleting directory %s' % work_dir)
    try:
        shutil.rmtree(work_dir)
        logging.debug('Directory %s deleted OK' % work_dir)
        return 0
    except:
        logging.error('Deleting directory %s failed' % work_dir)
        return 'Error'

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    commands = [ 'list', 'delete' , 'workspace' ]
    
    if len(sys.argv) < 2:
        print('usage: %s [list|delete <name>|workspace check|workspace delete]' % sys.argv[0])
        sys.exit(1)

    cfg = json.load(open('etc/conf.json'))

    if sys.argv[1] == 'list':
        cat = retrieve_catalog(cfg)
        print('%-60s desc' % 'id')
        print('-' * 70)
        for k in sorted(cat):
            print('%-60s %s' % (k, cat[k]['description']))
    elif sys.argv[1] == 'workspace':
        cat = retrieve_catalog(cfg)
        for f in glob.glob(osp.join(cfg['workspace_path'],'*')):
            if osp.isdir(f):
                ff = osp.basename(f)
                if ff not in cat:
                    logging.error('%s not found in the catalog' % ff)
                    if sys.argv[2] == 'delete':
                        local_rmdir(cfg, f)
    elif sys.argv[1] == 'delete':
        name = sys.argv[2]
        cat = retrieve_catalog(cfg)
        if name not in cat:
            logging.error('Cannot find simulation %s' % name)
        else:
            logging.info('Deleting simulation %s' % name)
            remote_rmdir(cfg, name)
            local_rmdir(cfg,name)
            del cat[name]
            store_catalog(cfg, cat)
            
    else:
        logging.error('command line not understood %s' % sys.argv)
        

