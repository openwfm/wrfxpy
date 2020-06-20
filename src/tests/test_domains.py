
from __future__ import absolute_import
from __future__ import print_function
from wrf.wps_domains import WPSDomainLCC
import json, netCDF4, sys
import numpy as np
from six.moves import range

d1 = WPSDomainLCC(1, { 'precomputed' : 'dom_test/geo_em.d01.nc', 'time_step' : 50 }, None )
d2 = WPSDomainLCC(2, { 'precomputed' : 'dom_test/geo_em.d02.nc', 'parent_time_step_ratio' : 4 }, d1)
d3 = WPSDomainLCC(3, { 'precomputed' : 'dom_test/geo_em.d03.nc', 'parent_time_step_ratio' : 4 }, d2)
domains = [d1,d2,d3]

def get_latlon(path):
    d = netCDF4.Dataset(path)
    lat = d.variables['XLAT_M'][0,:,:]
    lon = d.variables['XLONG_M'][0,:,:]
    d.close()
    return lat, lon

lat1, lon1 = get_latlon('dom_test/geo_em.d01.nc')
lat2, lon2 = get_latlon('dom_test/geo_em.d02.nc')
lat3, lon3 = get_latlon('dom_test/geo_em.d03.nc')
coords = [(lat1,lon1),(lat2,lon2),(lat3,lon3)]

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
