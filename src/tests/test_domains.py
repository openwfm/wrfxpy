
from __future__ import absolute_import
from __future__ import print_function
from wrf.wps_domains import WPSDomainLCC
import json, netCDF4
import numpy as np
from six.moves import range

d1 = WPSDomainLCC(1, { 'precomputed' : 'dom_test/geo_em.d01.nc', 'time_step' : 50 }, None )
d2 = WPSDomainLCC(2, { 'precomputed' : 'dom_test/geo_em.d02.nc', 'parent_time_step_ratio' : 4 }, d1)

def get_latlon(path):
    d = netCDF4.Dataset(path)
    lat = d.variables['XLAT_M'][0,:,:]
    lon = d.variables['XLONG_M'][0,:,:]
    d.close()
    return lat, lon

lat1, lon1 = get_latlon('dom_test/geo_em.d01.nc')
lat2, lon2 = get_latlon('dom_test/geo_em.d02.nc')

pari, parj = d2.parent_start[0]-1, d2.parent_start[1]-1
pcsr = d2.parent_cell_size_ratio
delta = (pcsr - 1) / 2

print(('Child domain starts at %d, %d pcsr=%d delta=%d' % (pari, parj, pcsr, delta)))

for i in range(3):
    for j in range(3):
        pi, pj = pari + i, parj + j
        ci, cj = i * pcsr + delta, j * pcsr + delta
        print(('i=%d j=%d parent: %g, %g child: %g, %g' % (i, j, lat1[pi,pj], lon1[pi,pj], lat2[ci,cj], lon2[ci,cj])))
        print((' ll->ij   projection from parent: %g, %g' % d1.latlon_to_ij(lat1[pi,pj], lon1[pi,pj])))
        print((' ll->ij   projection from child : %g, %g' % d2.latlon_to_ij(lat2[ci,cj], lon2[ci,cj])))
        print((' ij->ll   projection from parent: %g, %g' % (d1.ij_to_latlon(pi, pj))))
        print((' ij->ll   projection from child : %g, %g' % (d2.ij_to_latlon(ci, cj))))
        print('')
