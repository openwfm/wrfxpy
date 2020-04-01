# convert_geotiff.py
# Angel Farguell, March 2020

from __future__ import absolute_import
from geo.geodriver import GeoDriver
from geo.var_wisdom import get_wisdom_variables
from utils import file_exists
import logging, sys, os

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    if len(sys.argv) != 4:
        print('Usage: ./convert_geotiff.sh geotiff_file geogrid_folder var_name')
        print('Example: ./convert_geotiff.sh ./fuel.tif ./geo_data NFUEL_CAT')
        print('Available var_name options are %s' % get_wisdom_variables())
        exit(1)
    elif not file_exists(sys.argv[1]):
        print('File %s does not exist or is not readable' % sys.argv[1])
        exit(1)
    elif not sys.argv[3] in get_wisdom_variables():
        print('Variable %s not in available variables %s' % (sys.argv[3],get_wisdom_variables()))
        exit(1)

    _,file,path,var = sys.argv
    print('Generating geogrid folder %s from GeoTIFF file %s as variable %s' % (path,file,var))
    GeoDriver.from_file(file).to_geogrid(path,var)
