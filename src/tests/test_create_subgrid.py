from geo.geodriver import GeoDriver as GD
import pyproj,osr,sys
from netCDF4 import Dataset
import numpy as np

var = ['NFUEL_CAT', 'ZSF']
nc = Dataset(sys.argv[1])
lat1 = nc.getncattr('TRUELAT1')
lat2 = nc.getncattr('TRUELAT2')
lat0 = nc.getncattr('MOAD_CEN_LAT')
lon0 = nc.getncattr('STAND_LON')
lcc_proj4 = '+proj=lcc +lat_1=%.10f +lat_2=%.10f +lat_0=%.10f +lon_0=%.10f +a=6370000. +b=6370000. +towgs84=0,0,0 +no_defs' % (lat1,lat2,lat0,lon0)
lcc_proj = pyproj.Proj(lcc_proj4)
ref_proj = lcc_proj.to_latlong()
crs = osr.SpatialReference()
crs.ImportFromProj4(lcc_proj4)
dx = nc.getncattr('DX')
dy = nc.getncattr('DY')
srx = nc.getncattr('sr_x')
sry = nc.getncattr('sr_y')
dxs = dx/srx
dys = dy/sry
lon00 = nc.getncattr('corner_lons')[13]
lat00 = nc.getncattr('corner_lats')[13]
x00 = pyproj.transform(ref_proj,lcc_proj,lon00,lat00)
gt = (x00[0],dxs,0,x00[1],0,-dys)
for v in var:
    if v in nc.variables:
        array = np.flipud(nc[v][0,:-sry,:-srx])
        GD.from_elements(array,crs,gt).to_geotiff(v.lower()+'.tif')
