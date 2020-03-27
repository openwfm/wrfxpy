# convert_geotiff.py
# Angel Farguell, March 2020

from __future__ import absolute_import
import gdal,osr,pyproj,rasterio
import logging,sys
from utils import Dict
from write_geogrid import write_geogrid,addquotes

# WPS projections from proj4 attribute +proj
proj4_to_projwrf = {
    'lcc': 'lambert',
    'stere': 'polar',
    'merc': 'mercator',
    'latlong': 'regular_ll',
    'aea': 'albers_nad83',
    'stere': 'polar_wgs84'
}

def get_proj_str(rasterio_obj):
    """
    Get projection string for WPS from pyproj object

    :param rasterio_obj: rasterio object (projection information)
    :return proj_str: WRF project string 
    """
    proj = rasterio_obj.to_dict()['proj']
    proj_str = proj4_to_projwrf.get(proj,None)
    return proj_str

def create_index(ds,var):
    """
    Create index from gdal object and variable from WPS

    :param ds: geotiff file opened using gdal 
    :param var: variable name from WPS (NFUEL_CAT or ZSF)
    :return index: python dictionary with the entries of the index file
    """
    # get WKT string projection
    wkt = ds.GetProjectionRef() 
    # get Spatial Reference (projection)
    crs = osr.SpatialReference()
    crs.ImportFromWkt(wkt)
    # get Geo Transform
    gt = ds.GetGeoTransform()
    # get proj4 string
    proj4 = crs.ExportToProj4()
    # get rasterio object 
    rast_proj = rasterio.crs.CRS.from_proj4(proj4)
    # get pyproj element for tif file
    tif_proj = pyproj.Proj(proj4)
    # get pyproj element for WGS84
    ref_proj = pyproj.Proj(proj='lonlat',ellps='WGS84',datum='WGS84',no_defs=True)
    # define index dictionary
    index = Dict({})
    # projection string <- get_proj_str
    index.projection = get_proj_str(rast_proj)
    # truelat1 <- standard_parallel_1
    index.truelat1 = crs.GetProjParm("standard_parallel_1")
    # truelat2 <- standard_parallel_2
    index.truelat2 = crs.GetProjParm("standard_parallel_2")
    # stdlon <- longitude_of_center
    index.stdlon = crs.GetProjParm("longitude_of_center")
    # known_x, known_y <- size_x/2, size_y/2
    index.known_x = ds.RasterXSize//2
    index.known_y = ds.RasterYSize//2       
    # known_lon, known_lat <- trans(known_x), trans(known_y)
    posX, posY = gdal.ApplyGeoTransform(gt,index.known_x-1,index.known_y-1)
    index.known_lon,index.known_lat = pyproj.transform(tif_proj,ref_proj,posX,posY)
    # dx, dy <- from geotransform string
    index.dx = gt[1]
    index.dy = abs(gt[5])
    # tile_x, tile_y, tile_z <- RasterXSize, RasterYSize, 1
    index.tile_x = ds.RasterXSize
    index.tile_y = ds.RasterYSize
    index.tile_z = 1
    # units, description, missing_value, and others <- depending on var
    if var == 'NFUEL_CAT':
        index.units = addquotes("fuel category")
        index.description = addquotes("Anderson 13 fire behavior categories")
        index.category_min = 0
        index.category_max = 99
        index.missing_value = 99
    elif var == 'ZSF':
        index.units = addquotes("meters")
        index.description = addquotes("National Elevation Dataset 1/3 arcsecond resolution")
        index.missing_value = 0.0
    # row_order <- from sign of dy
    if index.dy > 0:
        index.row_order = 'top_bottom'
    else:
        index.row_order = 'bottom_top'
    return index

def read_geotiff(path,var):
    """
    Read GeoTIFF file using GDAL library

    :param path: path to the GeoTIFF file
    :param var: variable name from WPS (NFUEL_CAT or ZSF)
    :return array: numpy array with raster from GeoTIFF file
    :return index: python dictionary with the entries of the index file
    """
    ds = gdal.Open(path)
    array = ds.ReadAsArray()
    index = create_index(ds,var)
    return array,index

def geotiff2geogrid(path_out,path_file,var):
    """
    Transform GeoTIFF file into geogrid folder

    :param path_out: path to the geogrid folder to write
    :param path_file: path to the GeoTIFF file
    :param var: variable name from WPS (NFUEL_CAT or ZSF)
    """
    array,index = read_geotiff(path_file,var)
    if var == 'NFUEL_CAT':
        write_geogrid(path_out,array,index,bits=16,scale=1.,data_type='categorical')
    elif var == 'ZSF':
        write_geogrid(path_out,array,index,bits=16,scale=1.)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    if len(sys.argv) != 4:
        print('Usage: %s geotiff_file geogrid_folder var_name' % sys.argv[0])
        print('Example: %s ./fuel.tif ./fuel NFUEL_CAT' % sys.argv[0])
        exit(1)

    _,file,path,var = sys.argv
    print('Generating geogrid folder %s from GeoTIFF file %s' % (path,file))
    geotiff2geogrid(path,file,var)
