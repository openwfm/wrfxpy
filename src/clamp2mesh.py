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

def nearest_idx(lons,lats,x,y):
    flons = lons.flatten()
    flats = lats.flatten()
    xd = flons - x
    yd = flats - y
    idx = (xd*xd + yd*yd).argmin()
    return np.unravel_index(idx,lons.shape)

def array_filled(array):
    return bool(np.array(array).sum())

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
    # coarse sizes 
    nyc,nxc = lons.shape 
    # refined sizes
    nyr,nxr = (nyc+1)*sry,(nxc+1)*srx
    # real index of fire mesh of first atmosphere point
    i0 = (srx-1)*.5
    j0 = (sry-1)*.5
    # first fire index
    i0f = int(np.ceil(i0))
    j0f = int(np.ceil(j0))
    # initialize refined grid
    lonsr = np.zeros((nyr,nxr))
    latsr = np.zeros((nyr,nxr))
    # create bilinear coefficients
    tx,ty = np.meshgrid(np.arange(i0f-i0,srx,1)*rx,np.arange(j0f-j0,sry,1)*ry)
    t0 = np.multiply((1-tx),(1-ty))
    t1 = np.multiply((1-tx),ty)
    t2 = np.multiply(tx,(1-ty))
    t3 = np.multiply(tx,ty)
    # create bilinear function components
    lon0 = lons[:-1,:-1]
    lon1 = lons[1:,:-1]
    lon2 = lons[:-1,1:]
    lon3 = lons[1:,1:]
    lat0 = lats[:-1,:-1]
    lat1 = lats[1:,:-1]
    lat2 = lats[:-1,1:]
    lat3 = lats[1:,1:]
    # loop over refined position
    for y in range(t0.shape[0]):
        for x in range(t0.shape[1]):
            lonsr[j0f+y:nyr-2*sry+y:sry,i0f+x:nxr-2*srx+x:srx] = t0[y,x]*lon0+t1[y,x]*lon1+t2[y,x]*lon2+t3[y,x]*lon3
            latsr[j0f+y:nyr-2*sry+y:sry,i0f+x:nxr-2*srx+x:srx] = t0[y,x]*lat0+t1[y,x]*lat1+t2[y,x]*lat2+t3[y,x]*lat3

    return lonsr,latsr

def clamp2mesh(nc_path,x,y):
    """
    Return the closes node on the fire mesh
    :param nc_path: netcdf file with WRF coordinate arrays
    :param x, y: coordinates of the point
    :return x1,x2: nearby coordinates in the mesh 
    """

    d=nc4.Dataset(nc_path,'r') 
    varis = d.variables
    attrs = d.ncattrs()

    if 'sr_x' in attrs and 'sr_y' in attrs:
        srx = d.getncattr('sr_x')
        sry = d.getncattr('sr_y')
    else:
        srx = int(d.dimensions['west_east_subgrid'].size/(d.dimensions['west_east'].size+1))
        sry = int(d.dimensions['south_north_subgrid'].size/(d.dimensions['south_north'].size+1))

    if 'FXLONG' in varis and 'FXLAT' in varis and array_filled(d.variables['FXLONG']) and array_filled(d.variables['FXLAT']):
        print('> fxlong and fxlat exist')
        lons = np.array(d.variables['FXLONG'][0])
        lats = np.array(d.variables['FXLAT'][0])
    elif 'XLONG_M' in varis and 'XLAT_M' in varis:
        print('> fxlong and fxlat does not exist')
        lons_atm = np.array(d.variables['XLONG_M'][0])
        lats_atm = np.array(d.variables['XLAT_M'][0])
        print('> interpolating xlong_m to fxlong and xlat_m to fxlat...')
        lons,lats = interpolate_coords(lons_atm,lats_atm,srx,sry)
    elif 'XLONG' in varis and 'XLAT' in varis:
        print('> fxlong and fxlat does not exist')
        lons_atm = np.array(d.variables['XLONG'][0])
        lats_atm = np.array(d.variables['XLAT'][0])
        print('> interpolating xlong to fxlong and xlat to fxlat...')
        lons,lats = interpolate_coords(lons_atm,lats_atm,srx,sry)
    else:
        print('Error: %s NetCDF file specifiedc has not coordinates specified' % nc_path)
        sys.exit(1)

    idx = nearest_idx(lons,lats,x,y)
    return lons[idx], lats[idx] 
    

if __name__ == '__main__':

    # process arguments
    self, path, x, y  = sys.argv
    nlon, nlat =  clamp2mesh(path,float(x),float(y))
    print('%f %f\n' % (nlon, nlat))

