# clamp2mesh.py
# Jan Mandel March 2020

import netCDF4 as nc4
import numpy as np
import numpy.typing as npt
import sys
import logging


def nearest_idx(
        lons: np.array,
        lats: np.array,
        x: float,
        y: float) -> tuple:
    """Compute nearest index in the array to a point.

    Args:
        lons (np.array): longitude grid
        lats (np.array): latitude grid
        x (float): longitude point coordinate
        y (float): latitude point coordinate

    Returns:
        tuple: index tuple to slice grid
    """    
    flons = lons.flatten()
    flats = lats.flatten()
    xd = flons - x
    yd = flats - y
    idx = (xd*xd + yd*yd).argmin()
    return np.unravel_index(idx, lons.shape)


def array_filled(array: npt.ArrayLike) -> bool:
    """Assess if array is filled (variable with all 0s).

    Args:
        array (npt.ArrayLike): array to check for all 0s

    Returns:
        bool: True if the array was filled (some value different than 0)
    """
    return bool(np.array(array).sum())


def interpolate_coords(lons,lats,srx,sry,extrap=True):
    """
    Return the closes node on the mesh 
    :param lons,lats: longitude and latitude coordinate arrays for the coarse grid
    :param srx,sry: refinement
    :param extrap: extrapolation
    :return lonsr,latsr: longitude and latitude coordinate arrays for the refinement grid
    """
    def continue_at_boundary(lons,lats):
        # extrapolation, max quarded
        def EX(a,b,bias=0.):
            a = np.array(a)
            b = np.array(b)
            return np.squeeze((1.-bias)*(2.*a-b)+bias*np.max(np.c_[2.*a-b,a,b],axis=1))
        # continue at boundary
        clons = np.zeros(np.array(lons.shape)+2)
        clats = np.zeros(np.array(lats.shape)+2)
        clons[1:-1,1:-1] = lons
        clats[1:-1,1:-1] = lats
        clons[0,1:-1] = EX(lons[0],lons[1])
        clats[0,1:-1] = EX(lats[0],lats[1])
        clons[-1,1:-1] = EX(lons[-1],lons[-2])
        clats[-1,1:-1] = EX(lats[-1],lats[-2])
        clons[1:-1,0] = EX(lons[:,0],lons[:,1])
        clats[1:-1,0] = EX(lats[:,0],lats[:,1])
        clons[1:-1,-1] = EX(lons[:,-1],lons[:,-2])
        clats[1:-1,-1] = EX(lats[:,-1],lats[:,-2])
        clons[0,0] = EX(lons[0,0],lons[1,1])
        clats[0,0] = EX(lats[0,0],lats[1,1])
        clons[0,-1] = EX(lons[0,-1],lons[1,-2])
        clats[0,-1] = EX(lats[0,-1],lats[1,-2])
        clons[-1,0] = EX(lons[-1,0],lons[-2,1])
        clats[-1,0] = EX(lats[-1,0],lats[-2,1])
        clons[-1,-1] = EX(lons[-1,-1],lons[-2,-2])
        clats[-1,-1] = EX(lats[-1,-1],lats[-2,-2])
        return clons,clats

    # apply extrapolation if requiered
    if extrap:
        # continue at boundary for atmospheric mesh (extrapolation)
        lons,lats = continue_at_boundary(lons,lats)

    # coarse sizes 
    nyc,nxc = lons.shape 
    # refined sizes
    nyr,nxr = (nyc+1)*sry,(nxc+1)*srx
    # initialize refined grid
    lonsr = np.zeros((nyr,nxr))
    latsr = np.zeros((nyr,nxr))
    # real index of fire mesh of first atmosphere point
    i0 = (srx-1)*.5
    j0 = (sry-1)*.5
    # first fire index
    i0f = int(np.ceil(i0))
    j0f = int(np.ceil(j0))
    # last fire index
    i1f = nxr-2*srx
    j1f = nyr-2*sry
    # compute mesh ratios
    rx = 1./srx
    ry = 1./sry
    # create bilinear coefficients
    tx,ty = np.meshgrid(np.arange(.5,srx,1)*rx,np.arange(.5,sry,1)*ry)
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
            lonsr[j0f+y:j1f+y:sry,i0f+x:i1f+x:srx] = t0[y,x]*lon0+t1[y,x]*lon1+t2[y,x]*lon2+t3[y,x]*lon3
            latsr[j0f+y:j1f+y:sry,i0f+x:i1f+x:srx] = t0[y,x]*lat0+t1[y,x]*lat1+t2[y,x]*lat2+t3[y,x]*lat3

    if extrap:
        return lonsr[sry:-sry,srx:-srx],latsr[sry:-sry,srx:-srx]
    else:
        return lonsr,latsr


def get_subgrid_coordinates(nc_path):
    """
    Get subgrid coordinates from NetCDF file

    :param nc_path: path to NetCDF file with WRF coordinate arrays
    """

    with nc4.Dataset(nc_path) as d:
        varis = d.variables
        attrs = d.ncattrs()
        if 'sr_x' in attrs and 'sr_y' in attrs:
            srx = d.getncattr('sr_x')
            sry = d.getncattr('sr_y')
        else:
            srx = int(d.dimensions['west_east_subgrid'].size/(d.dimensions['west_east'].size+1))
            sry = int(d.dimensions['south_north_subgrid'].size/(d.dimensions['south_north'].size+1))        
        if 'FXLONG' in varis and 'FXLAT' in varis and array_filled(d.variables['FXLONG']) and array_filled(d.variables['FXLAT']):
            lons = np.array(d.variables['FXLONG'][0])
            lats = np.array(d.variables['FXLAT'][0])
            return lons,lats
        else:
            logging.info('get subgrid coordinates from netCDF file...')
            if 'XLONG_M' in varis and 'XLAT_M' in varis:
                lons_atm = np.array(d.variables['XLONG_M'][0])
                lats_atm = np.array(d.variables['XLAT_M'][0])
            elif 'XLONG' in varis and 'XLAT' in varis:
                lons_atm = np.array(d.variables['XLONG'][0])
                lats_atm = np.array(d.variables['XLAT'][0])
            else:
                raise TypeError('atmospheric coordinates not found in file, skipping')
            lons,lats = interpolate_coords(lons_atm,lats_atm,srx,sry)

    return lons,lats


def fill_subgrid(nc_path):
    """
    Fill netCDF file with refined coordinates, if necessary
    :param nc_path: NetCDF file with WRF coordinate arrays
    """
    
    logging.info(f'filling subgrid coordinates in NetCDF file {nc_path}...')
    with nc4.Dataset(nc_path, 'a') as d:
        varis = d.variables
        if 'FXLONG' in varis and 'FXLAT' in varis and array_filled(d.variables['FXLONG']) and array_filled(d.variables['FXLAT']):
            logging.info('subgrid coordinates already defined')
            return
        lons,lats = get_subgrid_coordinates(nc_path)
        if 'FXLONG' in d.variables.keys():
            d['FXLONG'][:] = lons[np.newaxis,:,:]
            d['FXLAT'][:] = lats[np.newaxis,:,:]
        else:
            fxlong = d.createVariable("FXLONG","f4",("Time","south_north_subgrid","west_east_subgrid")) 
            fxlat = d.createVariable("FXLAT","f4",("Time","south_north_subgrid","west_east_subgrid")) 
            fxlong[:] = lons[np.newaxis,:,:]
            fxlat[:] = lats[np.newaxis,:,:]
    return


def clamp2mesh(nc_path,x,y):
    """
    Return the closes node on the fire mesh
    :param nc_path: netcdf file with WRF coordinate arrays
    :param x, y: coordinates of the point
    :return x1,x2: nearby coordinates in the mesh 
    """

    lons,lats = get_subgrid_coordinates(nc_path)
    idx = nearest_idx(lons,lats,x,y)
    return lons[idx], lats[idx] 
    

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # process arguments
    if len(sys.argv) == 4:
        self, path, x, y  = sys.argv
        nlon, nlat =  clamp2mesh(path,float(x),float(y))
        logging.info('%f %f\n' % (nlon, nlat))
    elif len(sys.argv) == 2:
        self, path = sys.argv
        fill_subgrid(path)
        logging.info('filled subgrid coordinates in {}'.format(path))
