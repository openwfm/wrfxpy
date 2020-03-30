# geotiff.py
# Angel Farguell, March 2020

import gdal, osr, pyproj, rasterio, gdalconst
import logging

from utils import Dict
from geo.write_geogrid import write_geogrid_var
from geo.geo_utils import coord2str

class GeoTIFF(object):
    """
    Represents the content of one GeoTIFF file.
    """

    def __init__(self, path):
        """
        Initializes metadata information.

        :param path: path to the file
        """
        self.path = path
        logging.info('Reading GeoTIFF file %s' % path)
        # open GeoTIFF file
        self.dgdal = gdal.Open(path,gdalconst.GA_ReadOnly)
        # raster size
        self.nx = self.dgdal.RasterXSize
        self.ny = self.dgdal.RasterYSize
        # define projection objects
        self.init_proj()

    def close(self):
        """
        Close the current file.
        """
        del self.dgdal

    def init_proj(self):
        """
        Get all necessary projection information.
        """
        # get WKT string projection
        self.wkt = self.dgdal.GetProjectionRef()
        # get Spatial Reference (projection)
        self.crs = osr.SpatialReference()
        self.crs.ImportFromWkt(self.wkt)
        # get Geo Transform
        self.gt = self.dgdal.GetGeoTransform()
        # get proj4 string
        self.proj4 = self.crs.ExportToProj4()
        # get rasterio object 
        self.rasterio = rasterio.crs.CRS.from_proj4(self.proj4)
        # get pyproj element for tif file
        self.pyproj = pyproj.Proj(self.proj4)
        # projection short string
        self.projstr = self.rasterio.to_dict()['proj']
        # WPS projections from proj4 attribute +proj
        self.projwrf = {'lcc': 'lambert',
            'stere': 'polar',
            'merc': 'mercator',
            'latlong': 'regular_ll',
            'aea': 'albers_nad83',
            'stere': 'polar_wgs84'
        }.get(self.projstr,None)
        # If not found, resample to lambert projection
        if not self.projwrf:
            self.projwrf = 'lambert'
            self.reproj = True
        else:
            self.reproj = False

    def get_array(self,resample=None):
        """
        Get array information.
        """
        if resample:
            logging.warning('GeoTIFF.get_array(resample=coord) is not implemented yet')
            # resample and maybe reproject (if self.reproj, then reproject and set new projection variables)
        else:
            return self.dgdal.ReadAsArray()

    def geogrid_index(self):
        """
        Geolocation in a form suitable for geogrid index.
        :return: dictionary key:value 
        """
        # known points in the center of the array (mass-staggered mesh from 1 at origin)
        known_x = (self.nx+1)/2 
        known_y = (self.ny+1)/2 
        # known points in the center of the array (unstaggered mesh from 0 at origin)
        known_x_uns = known_x-.5
        known_y_uns = known_y-.5
        # get pyproj element for WGS84
        ref_proj = pyproj.Proj(proj='lonlat',ellps='WGS84',datum='WGS84',no_defs=True)
        # lat/lon known points
        posX, posY = gdal.ApplyGeoTransform(self.gt,known_x_uns,known_y_uns)
        known_lon,known_lat = pyproj.transform(self.pyproj,ref_proj,posX,posY) 
        # print ceter lon-lat coordinates to check with gdalinfo
        logging.info('GeoTIFF.geogrid_index - center lon/lat coordinates using pyproj.transform: %s' % coord2str(known_lon,known_lat)) 
        logging.info('GeoTIFF.geogrid_index - center lon/lat coordinates using gdal.Info: %s' % gdal.Info(self.path).split('Center')[1].split(') (')[1].split(')')[0])  
        # row order depending on sign of dy
        if self.gt[5] > 0:
            row_order = 'bottom_top'
        else:
            row_order = 'top_bottom'
        # write geogrid index
        return Dict({'projection' : self.projwrf,
            'dx' : self.gt[1],
            'dy' : abs(self.gt[5]),
            'truelat1' : self.crs.GetProjParm("standard_parallel_1"),
            'truelat2' : self.crs.GetProjParm("standard_parallel_2"),
            'stdlon' : self.crs.GetProjParm("longitude_of_center"),
            'known_x' : known_x,
            'known_y' : known_y,
            'known_lon' : known_lon,
            'known_lat' : known_lat,
            'row_order' : row_order
        })

    def to_geogrid(self,path,var):
        """
        Transform to geogrid files
        :param path: path to write the geogrid files
        :param var: variable name for WPS (available options in src/geo/var_wisdom.py)
        """
        logging.info('GeoTIFF.to_geogrid() - getting array')
        array = self.get_array()
        logging.info('GeoTIFF.to_geogrid() - creating index')
        index = self.geogrid_index()
        logging.info('GeoTIFF.to_geogrid() - writting geogrid')
        write_geogrid_var(path,var,array,index)

    def __str__(self):
        """
        Returns a string representation of the message as provided by gdal.
        
        :return: descriptive string
        """
        return str(self.dgdal)