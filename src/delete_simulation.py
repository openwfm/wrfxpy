from ssh_shuttle import delete_product_from_server
import json
import logging
from utils import Dict
import os
import os.path as osp
import glob
import shutil

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print len(sys.argv) 
        print sys.argv
        print('usage: %s <sim-name>' % sys.argv[0])
        sys.exit(1)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    sim_name = sys.argv[1]
    logging.info('Removing %s' % sim_name)
    cfg = Dict(json.load(open('etc/conf.json')))
    try:
        delete_product_from_server(cfg, sim_name)
        logging.info('Remote delete OK')
    except:
        logging.error('Remote delete failed,')
    work_dir = osp.abspath(osp.join(cfg.workspace_path, sim_name))
    logging.info('Removing work tree %s' % work_dir)
    try:
        shutil.rmtree(work_dir)
        logging.info('Work tree delete OK')
    except:
        logging.error('Work tree delete failed,')

