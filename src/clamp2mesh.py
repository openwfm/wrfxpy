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

def interpolate_coords(lons,lats,srx,sry):
    """
    Return the closes node on the mesh 
    :param lons,lats: longitude and latitude coordinate arrays for the coarse grid
    :param srx,sry: refinement
    :return lonsr,latsr: longitude and latitude coordinate arrays for the refinement grid
    """
    # compute mesh ratios
    rx = 1./srx
    ry = 1./sry
    # coarse indices
    isc,iec,jsc,jec = 0,lons.shape[1]-1,0,lons.shape[0]-1
    # fine indices
    isf,ief,jsf,jef = isc,(iec-isc+2)*srx-1,jsc,(jec-jsc+2)*sry-1
    # both meshes lined up by lower left corner of domain
    ripc,rjpc,ripf,rjpf = float(isc),float(jsc),isc+(srx-1)*.5,jsc+(sry-1)*.5
    # boundary lenghts
    ih = int(np.ceil(ripf))
    jh = int(np.ceil(rjpf))
    # initialize fire mesh
    lonsr = np.zeros((jef+1,ief+1))
    latsr = np.zeros((jef+1,ief+1))

    # bilinear interpolation in the inside
    for j_c in range(jsc,jec):
        rjo = rjpf+sry*(j_c-rjpc)
        j_s = max(jsf,int(np.ceil(rjo)))
        j_e = min(jef,int(np.floor(rjo)+sry))
        for i_c in range(isc,iec):
            rio = ripf+srx*(i_c-ripc)
            i_s = max(isf,int(np.ceil(rio)))
            i_e = min(ief,int(np.floor(rio)+srx))
            for j_f in range(j_s,j_e+1):
                ty = (j_f-rjo)*ry
                for i_f in range(i_s,i_e+1):
                    tx = (i_f-rio)*rx
                    lonsr[j_f,i_f] = (1-tx)*(1-ty)*lons[j_c,i_c]+ \
                                    (1-tx)*ty*lons[j_c+1,i_c]+ \
                                    tx*(1-ty)*lons[j_c,i_c+1]+ \
                                    tx*ty*lons[j_c+1,i_c+1]
                    latsr[j_f,i_f] = (1-tx)*(1-ty)*lats[j_c,i_c]+ \
                                    (1-tx)*ty*lats[j_c+1,i_c]+ \
                                    tx*(1-ty)*lats[j_c,i_c+1]+ \
                                    tx*ty*lats[j_c+1,i_c+1]

    return lonsr,latsr

def clamp2mesh(nc_path,x,y):
    """
    Return the closes node on the fire mesh
    :param nc_path: netcdf file with WRF coordinate arrays
    :param x, y: coordinates of the point
    :return x1,x2: nearby coordinates in the mesh 
    """

    d=nc4.Dataset(nc_path) 
    varis = d.variables
    attrs = d.ncattrs()
    if 'FXLONG' in varis and 'FXLAT' in varis:
        print('> fxlong and fxlat exist')
        lons = np.array(d.variables['FXLONG'][0])
        lats = np.array(d.variables['FXLAT'][0])
    elif 'XLONG_M' in varis and 'XLAT_M' in varis and 'sr_x' in attrs and 'sr_y' in attrs:
        print('> fxlong and fxlat does not exist')
        lons_atm = np.array(d.variables['XLONG_M'][0])
        lats_atm = np.array(d.variables['XLAT_M'][0])
        srx = d.getncattr('sr_x')
        sry = d.getncattr('sr_y')
        print('> interpolating xlong to fxlong and xlat to fxlat...')
        lons,lats = interpolate_coords(lons_atm,lats_atm,srx,sry)
    elif 'XLONG' in varis and 'XLAT' in varis and 'sr_x' in attrs and 'sr_y' in attrs:
        print('> fxlong and fxlat does not exist')
        lons_atm = np.array(d.variables['XLONG'][0])
        lats_atm = np.array(d.variables['XLAT'][0])
        srx = d.getncattr('sr_x')
        sry = d.getncattr('sr_y')
        print('> interpolating xlong to fxlong and xlat to fxlat...')
        lats,lons = interpolate_coords(lons_atm,lats_atm,srx,sry)
    else:
        print('Error: %s NetCDF file specifiedc has not coordinates specified')

    lons = lons.flatten()
    lats = lats.flatten()
    xd = lons - x
    yd = lats - y
    idx = (xd*xd + yd*yd).argmin()
    return lats[idx], lons[idx] 
    

if __name__ == '__main__':

    # process arguments
    self, path, x, y  = sys.argv
    nlon, nlat =  clamp2mesh(path,float(x),float(y))
    print('%f %f\n' % (nlon, nlat))

