from __future__ import absolute_import
from clamp2mesh import fill_subgrid
import sys, logging

if len(sys.argv) < 2: 
    raise SystemExit('usage: ./fill_subgrid.sh wrf_path')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
fill_subgrid(sys.argv[1])
