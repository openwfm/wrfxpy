# clamp2mesh.py
# Jan Mandel March 2020

from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os.path as osp
import os
import netCDF4 as nc4
import numpy as np

def clamp2mesh(nc_path,x,y):
    """
    Return the closes node on the mesh 
    :param nc_path: netcdf file with fxlong and fxlat arrays
    :param x, y: coordinates of the point
    :return x1,x2: nearby coordinates in the mesh 
    """

    d=nc4.Dataset(nc_path) 
    lats = d.variables['FXLAT']
    lons = d.variables['FXLONG']
    lats = np.array(lats)
    lons = np.array(lons)
    lons = lons.flatten()
    lats = lats.flatten()
    xd = lats - x 
    yd = lons - y
    idx = (xd*xd + yd*yd).argmin()
    return lats[idx], lons[idx] 
    

if __name__ == '__main__':

    # process arguments
    self, path, x, y  = sys.argv
    nlon, nlat =  clamp2mesh(path,float(x),float(y))
    print('%f %f\n' % (nlon, nlat))

