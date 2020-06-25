
from __future__ import absolute_import
from __future__ import print_function
from wrf.wps_domains import WPSDomainLCC
import json, netCDF4, os, sys, glob
import numpy as np
import os.path as osp
from six.moves import range

if len(sys.argv) != 2:
    print('Usage: python %s path' % sys.argv[0])
    print('            path must be a folder with geo_em.d0?.nc files OR')
    print('                           JSON file with domains information on it')
    sys.exit(1)
else:     
    PATH = sys.argv[1]
    if osp.isfile(PATH) and os.access(PATH, os.R_OK):
        pass
    else:
        files = glob.glob(osp.join(PATH,'geo_em.d0?.nc'))
        if osp.isdir(PATH) and len(files):
            files = sorted(glob.glob(osp.join(PATH,'geo_em.d0?.nc')))
            domains = []
            for d,file in enumerate(files):
                if d == 0:
                    domains.append(WPSDomainLCC(1, { 'precomputed' : file, 'time_step' : 50 }, None ))
                else:
                    domains.append(WPSDomainLCC(d+1, { 'precomputed' : file, 'parent_time_step_ratio' : 4 }, domains[d-1]))
        else:
            print('Error: something wrong with path provided %s' % sys.argv[1])
            print('        path must be a folder with geo_em.d0?.nc files OR')
            print('                       JSON file with domains information on it')
            sys.exit(1)

def get_latlon(path):
    d = netCDF4.Dataset(path)
    lat = d.variables['XLAT_M'][0,:,:]
    lon = d.variables['XLONG_M'][0,:,:]
    d.close()
    return lat, lon

coords = []
for file in files:
    coords.append(get_latlon(file))

err = 0
for k,d in enumerate(domains):
    print('>> Domain 0{} <<'.format(k+1))
    ilist = [0,d.domain_size[0]-2,d.domain_size[0]-2,0,5]
    jlist = [0,0,d.domain_size[1]-2,d.domain_size[1]-2,10]
    print('ij -> ll')
    for c in range(len(ilist)):
        i = ilist[c]
        j = jlist[c]
        latlon = d.ij_to_latlon(i,j)
        print('i={0},j={1}:  latlon=({2}, {3}) - ij_to_latlon=({4}, {5})'.format(i, j, coords[k][0][j,i], coords[k][1][j,i], latlon[0], latlon[1]))
        err += np.sqrt((coords[k][0][j,i]-latlon[0])**2+(coords[k][1][j,i]-latlon[1])**2)
    print('ll -> ij')
    for c in range(len(ilist)):
        i = ilist[c]
        j = jlist[c]
        ij = d.latlon_to_ij(coords[k][0][j,i],coords[k][1][j,i])
        print('i={0},j={1}:  latlon_to_ij=({2},{3})'.format(i, j, round(ij[0]), round(ij[1])))    

print('>> Summary coordinate errors <<')
print('sum(sqrt(lon_diff**2+lat_diff**2))={}'.format(err))

print('>> Bounding box <<')
for d in domains:
    bbox = d.bounding_box()
    lons = [b[1] for b in bbox]
    lats = [b[0] for b in bbox]
    print('{}'.format(bbox))
    print('bbox={0},{1},{2},{3}'.format(min(lons),max(lons),min(lats),max(lats)))

