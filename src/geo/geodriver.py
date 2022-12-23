# geodriver.py
# Angel Farguell, March 2020

from osgeo import gdal, osr
import pyproj
import logging
from utils import Dict
import numpy as np
from geo.write_geogrid import write_geogrid_var
from geo.geo_utils import coord2str

class GeoDriverError(Exception):
    """
    Raised when a GeoDriver failed.
    """
    pass

class GeoDriver(object):
    """
    Represents the content of one GeoDriver file.
    """

    def __init__(self, gdal_ds):
        """
        Initializes metadata information.

        :param gdal_ds: gdal dataset object
        """
        if not isinstance(gdal_ds,gdal.Dataset):
            raise GeoDriverError('GeoDriver: need a gdal.Dataset object instead of %s.\nTry using: GeoDriver.from_file or GeoDriver.from_elements' % type(gdal_ds))
        self.ds = gdal_ds
        # raster size and spacing
        self.nx = self.ds.RasterXSize
        self.ny = self.ds.RasterYSize
        self.bands = self.ds.RasterCount
        self.resampled = False
        # define projection objects
        self.init_proj()

    def close(self):
        """
        Close the current file.
        """
        self.ds = None

    def init_proj(self):
        """
        Get all necessary projection information.
        """
        # get WKT string projection
        self.wkt = self.ds.GetProjectionRef()
        # get Spatial Reference (projection)
        self.crs = osr.SpatialReference()
        self.crs.ImportFromWkt(self.wkt)
        # get Geo Transform
        self.gt = self.ds.GetGeoTransform()
        # get proj4 string
        self.proj4 = self.crs.ExportToProj4()
        # get pyproj element for tif file
        self.pyproj = pyproj.Proj(self.proj4)
        # projection short string
        proj_list = [p for p in self.proj4.split(' ') if 'proj' in p]
        if len(proj_list):
            self.projstr = proj_list[0].split('=')[-1]
        else:
            self.projstr = ''
        # WPS projections from proj4 attribute +proj
        self.projwrf = {'lcc': 'lambert',
            'stere': 'polar',
            'merc': 'mercator',
            'longlat': 'regular_ll',
            'aea': 'albers_nad83',
            'stere': 'polar_wgs84'
        }.get(self.projstr,None)
        # If not found, resample to lambert projection
        if not self.projwrf:
            self.projwrf = 'lambert'
            self.reproj = True
            # the GeoTIFF needs to be reprojected to lambert projection
            raise GeoDriverError('GeoTIFF.init_proj() with reprojection is not implemented yet')
        else:
            self.reproj = False

    def get_array(self,bbox=None):
        """
        Get array information.

        :param bbox: optional, WGS84 bounding box (min_lon,max_lon,min_lat,max_lat)
        """
        if isinstance(bbox,(tuple,list,np.ndarray)):
            return self.resample_bbox(bbox)
        else:
            return self.ds.ReadAsArray()

    def get_coord(self):
        """
        Get coordinate arrays for the whole GeoTIFF file using initial projection.
        """
        x0,dx,_,y0,_,dy = self.gt
        xx = np.linspace(x0,x0+dx*self.nx,self.nx)
        yy = np.linspace(y0,y0+dy*self.ny,self.ny)
        return np.meshgrid(xx,yy)

    def resample_bbox(self,bbox):
        """
        Resample using lon-lat bounding box.

        :param bbox: optional, WGS84 bounding box (min_lon,max_lon,min_lat,max_lat)
        """
        logging.info('GeoTIFF.resample_bbox - resampling GeoDriver into bounding box: %s' % list(bbox))  
        # geotransform
        x0,dx,_,y0,_,dy = self.gt
        xx = np.arange(x0,x0+dx*self.nx,dx)
        yy = np.arange(y0,y0+dy*self.ny,dy)
        # get pyproj element for WGS84
        ref_proj = pyproj.Proj(proj='longlat',ellps='WGS84',datum='WGS84',no_defs=True)
        # get bounding box in projection of the GeoTIFF
        ref_corners = ((bbox[0],bbox[2]),(bbox[0],bbox[3]),(bbox[1],bbox[3]),(bbox[1],bbox[2]))
        proj_corners = [pyproj.transform(ref_proj,self.pyproj,c[0],c[1]) for c in ref_corners]
        x_min = min([c[0] for c in proj_corners])
        x_max = max([c[0] for c in proj_corners])
        y_min = min([c[1] for c in proj_corners])
        y_max = max([c[1] for c in proj_corners])
        # find resample indexes
        i_mins = np.where(x_min <= xx)[0]
        i_maxs = np.where(xx <= x_max)[0]
        j_mins = np.where(y_min <= yy)[0]
        j_maxs = np.where(yy <= y_max)[0]
        if len(i_mins)==0 or len(i_maxs)==0 or len(j_mins)==0 or len(j_maxs)==0:
            logging.error('GeoTIFF.resample_bbox - bbox provided does not intersect the image')
            raise GeoDriverError('Incorrect bbox {}'.format(bbox))
        if dx > 0:
            i_min = i_mins.min()
            i_max = i_maxs.max()
        else:
            i_max = i_mins.max()
            i_min = i_maxs.min()
        if dy > 0:
            j_min = j_mins.min()
            j_max = j_maxs.max()
        else:
            j_max = j_mins.max()
            j_min = j_maxs.min()
        # save resample indexes
        self.resample_indxs = (i_min,i_max,j_min,j_max)
        # only working on Linux
        try:
            # get virtual array and resample
            if self.bands == 1:
                a_r = self.ds.GetVirtualMemArray()[j_min:j_max,i_min:i_max].copy()
            else:
                a_r = self.ds.GetVirtualMemArray()[:,j_min:j_max,i_min:i_max].copy()
        except:
            logging.warning('GeoTIFF.resample_bbox - reading the whole array and sampling after')
            # get array and resample
            if self.bands == 1:
                a_r = self.ds.ReadAsArray()[j_min:j_max,i_min:i_max]
            else:
                a_r = self.ds.ReadAsArray()[:,j_min:j_max,i_min:i_max]
        # update other elements
        self.gt = (xx[i_min],dx,0,yy[j_min],0,dy)
        if self.bands == 1:
            self.ny,self.nx = a_r.shape
        else:
            _,self.ny,self.nx = a_r.shape
        self.resampled = True
        self.bbox = bbox
        return a_r

    def geogrid_index(self):
        """
        Geolocation in a form suitable for geogrid index.

        :return: dictionary key:value 
        """
        # spacing
        dx,dy = self.gt[1],self.gt[5]
        # known points in the center of the array (mass-staggered mesh from 1 at origin)
        known_x = (self.nx+1)/2. 
        known_y = (self.ny+1)/2. 
        # known points in the center of the array (unstaggered mesh from 0 at origin)
        known_x_uns = known_x-.5
        known_y_uns = known_y-.5
        # get pyproj element for reference lat-long coordinates
        ref_proj = self.pyproj.to_latlong()
        # lat/lon known points
        posX, posY = gdal.ApplyGeoTransform(self.gt,known_x_uns,known_y_uns)
        known_lon,known_lat = pyproj.transform(self.pyproj,ref_proj,posX,posY) 
        # print center lon-lat coordinates to check with gdalinfo
        if not self.resampled:
            logging.info('GeoTIFF.geogrid_index - center lon/lat coordinates using pyproj.transform: %s' % coord2str(known_lon,known_lat)) 
            logging.info('GeoTIFF.geogrid_index - center lon/lat coordinates using gdal.Info: %s' % gdal.Info(self.ds).split('Center')[1].split(') (')[1].split(')')[0])  
        # row order depends on dy
        if dy > 0:
            row_order = 'bottom_top'
        else:
            row_order = 'top_bottom'
        # write geogrid index
        return Dict({'projection' : self.projwrf,
            'dx' : dx,
            'dy' : dy,
            'truelat1' : self.crs.GetProjParm("standard_parallel_1"),
            'truelat2' : self.crs.GetProjParm("standard_parallel_2"),
            'stdlon' : self.crs.GetProjParm("longitude_of_center",self.crs.GetProjParm("central_meridian",self.crs.GetProjParm("longitude_of_origin"))),
            'known_x' : known_x,
            'known_y' : known_y,
            'known_lon' : known_lon,
            'known_lat' : known_lat,
            'row_order' : row_order
        })

    def to_geogrid(self,path,var,bbox=None):
        """
        Transform to geogrid files

        :param path: path to write the geogrid files
        :param var: variable name for WPS (available options in src/geo/var_wisdom.py)
        :param bbox: optional, WGS84 bounding box (min_lon,max_lon,min_lat,max_lat)
        """
        logging.info('GeoTIFF.to_geogrid() - getting array')
        array = self.get_array(bbox)
        if len(array.shape) > 2:
            raise GeoDriverError('GeoTIFF.to_geogrid() - generation of 3d geogrid arrays is not supported')
        # WPS is generating regular lat-lon arrays upside down
        if self.projwrf == 'regular_ll':
            array = np.flipud(array)
        logging.info('GeoTIFF.to_geogrid() - creating index')
        index = self.geogrid_index()
        logging.info('GeoTIFF.to_geogrid() - writting geogrid')
        write_geogrid_var(path,var,array,index,coord=self.get_coord())

    def to_geotiff(self,path,**kwargs):
        """
        Transform to geotiff file

        :param path: path to write the geotiff file
        :param **kwargs: dictionary with optional arguments:
               - bbox: WGS84 bounding box (min_lon,max_lon,min_lat,max_lat)
               - unit: string with units information
               - desc: string with description information
               - ndv: number to define the non-data value
        """
        bbox = kwargs.get('bbox')
        logging.info('GeoTIFF.to_geotiff() - getting array')
        array = self.get_array(bbox)
        dims = array.shape
        if len(dims) == 2:
            array = array.reshape((1,)+dims)
        logging.info('GeoTIFF.to_geotiff() - writting geotiff')
        ds = gdal.GetDriverByName('GTiff').Create(path,self.nx,self.ny,self.bands,gdal.GDT_Float32)
        ds.SetGeoTransform(self.gt)
        ds.SetProjection(self.crs.ExportToWkt())
        for b in range(self.bands):
            band = ds.GetRasterBand(b+1)
            band.WriteArray(array[b,:,:])
            if 'unit' in kwargs.keys():
                if isinstance(kwargs['unit'],(tuple,list,np.ndarray)):
                    band.SetUnitType(kwargs['unit'][b])
                else:
                    band.SetUnitType(kwargs['unit'])
            if 'desc' in kwargs.keys():
                if isinstance(kwargs['desc'],(tuple,list,np.ndarray)):
                    band.SetDescription(kwargs['desc'][b])
                else:
                    band.SetDescription(kwargs['desc'])
            if 'ndv' in kwargs.keys():
                if isinstance(kwargs['ndv'],(tuple,list,np.ndarray)):
                    band.SetNoDataValue(kwargs['ndv'][b])
                else:
                    band.SetNoDataValue(kwargs['ndv'])
        ds.FlushCache()
        ds = None

    def __str__(self):
        """
        Returns a string representation of the message as provided by gdal.
        
        :return: descriptive string
        """
        return str(self.ds)

    @classmethod
    def from_file(cls, path):
        """
        Construct a GeoDriver from file.
        
        :param path: the path to the file
        """
        logging.info('Reading file %s' % path)
        # open GeoTIFF file
        ds = gdal.Open(path)
        if ds:
            gd = cls(ds)
            gd.path = path
        else:
            gd = None
        return gd

    @classmethod
    def from_elements(cls, data, crs, gt):
        """
        Construct a GeoDriver from file.
        
        :param data: data array 
        :oaran crs: spatial reference
        :param gt: geotransform tuple
        """
        dims = data.shape
        if len(dims) == 2:
            data = data.reshape((1,)+dims)
            bands,rows,cols = data.shape
        elif len(dims) == 3:
            bands,rows,cols = dims
        else:
            raise GeoDriverError('GeoDriver.from_elements - data provided must be a 2d or 3d array')
        ds = gdal.GetDriverByName('MEM').Create('',cols,rows,bands,gdal.GDT_Float32)
        ds.SetProjection(crs.ExportToWkt())
        ds.SetGeoTransform(gt)
        for b in range(bands):
            ds.GetRasterBand(b+1).WriteArray(data[b,:,:])
        gd = cls(ds)
        return gd
