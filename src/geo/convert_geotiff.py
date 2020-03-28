# convert_geotiff.py
# Angel Farguell, March 2020

from __future__ import absolute_import
from geo.geotiff import GeoTIFF
import logging, sys

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    if len(sys.argv) != 4:
        print('Usage: %s geotiff_file geogrid_folder var_name' % sys.argv[0])
        print('Example: %s ./fuel.tif ./fuel NFUEL_CAT' % sys.argv[0])
        exit(1)

    _,file,path,var = sys.argv
    print('Generating geogrid folder %s from GeoTIFF file %s as variable %s' % (path,file,var))
    tif = GeoTIFF(file)
    tif.to_geogrid(path,var)
