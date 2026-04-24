#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
import matplotlib as mpl
from six.moves import range
from six.moves import zip
mpl.use('AGG')
import simplekml as kml
import numpy as np
import netCDF4 as nc4
from pyhdf.SD import SD, SDC
import re
import h5py
import sys
import os
import json
import time
import logging
import os.path as osp
import traceback
from datetime import timedelta
from subprocess import check_call
from utils import dump, traceargs, esmf_to_utc, utc_to_esmf
from vis.vis_utils import print_stats
from vis.rasterizer import make_colorbar, make_discrete_colorbar, basemap_raster_mercator, basemap_barbs_mercator, basemap_scatter_mercator
from vis.var_wisdom import convert_value, get_wisdom, is_windvec, is_fire_var, strip_end
from geo.geo_utils import ncwrfmeta
from geo.geodriver import GeoDriver

class PostprocError(Exception):
    pass

def scalar_field_to_raster(fa, lats, lons, wisdom):
    """
    Render a scalar variable into a geolocated raster and colorbar.

    :param fa: the field to render
    :param lats: the latitudes
    :param lons: the longitudes
    :param wisdom: a configuration dictionary controlling the visualization
    :return: a tuple with the raster as StringIO, its geolocation and the PNG colorbar StringIO (None if not requested)
    """
    # gather wisdom about the variable
    native_unit = wisdom['native_unit']
    cmap_name = wisdom['colormap']
    cmap = mpl.cm.get_cmap(cmap_name)

    if lats.shape != fa.shape:
        raise PostprocError("Variable size %s does not correspond to grid size %s." % (fa.shape, lats.shape))

    logging.info('scalar_field_to_raster: variable %s min %s max %s' % (wisdom['name'], np.nanmin(fa),np.nanmax(fa)))

    # mask 'transparent' color value
    if 'transparent_values' in wisdom:
        rng = wisdom['transparent_values']
        fa = np.ma.masked_array(fa, np.logical_and(fa >= rng[0], fa <= rng[1]))
        logging.info('scalar_field_to_raster: transparent from %s to %s, elements %s not masked %s' % (rng[0], rng[1], fa.size , fa.count()))
        logging.info('scalar_field_to_raster: masked variable %s min %s max %s' % (wisdom['name'], np.nanmin(fa),np.nanmax(fa)))
    else:
        fa = np.ma.masked_array(fa)

    # look at mins and maxes, transparent don't count
    if fa.count():
        fa_min,fa_max = np.nanmin(fa),np.nanmax(fa)
    else:
        fa_min, fa_max = 0.0, 0.0

    # determine if we will use the range in the variable or a fixed range
    scale = wisdom['scale']
    if scale != 'original':
        m = fa.mask.copy()
        fa_min, fa_max = scale[0], scale[1]
        fa[fa < fa_min] = fa_min
        fa[fa > fa_max] = fa_max
        fa.mask = m

    # define distribution of colormap if necessary 
    norm = None
    ticks = None
    ticklabels = None
    if 'norm_opt' in wisdom: 
        norm_opt = wisdom['norm_opt']
        if norm_opt == 'lognorm':
            linthresh = wisdom.get('linthresh',.01)
            linscale = wisdom.get('linscale',.001)
            norm = lambda xmin,xmax: mpl.colors.SymLogNorm(linthresh=linthresh,linscale=linscale,vmin=xmin,vmax=xmax,base=10) 
        elif norm_opt == 'boundary':
            bounds = wisdom.get('bounds',[0,1,2,4,6,8,12,16,20,25,30,40,60,100,200])
            if 'labels' in wisdom and wisdom['labels'] is not None:
                ticklabels = wisdom['labels']
            elif wisdom['colorbar'] is not None:
                cb_unit = wisdom['colorbar']
                ticklabels = [convert_value(native_unit, cb_unit, b) for b in bounds[1:]]
            ticks = bounds[1:]
            bounds = bounds + [1e10]
            colors = wisdom.get('colors',np.array([(255,255,255),(197,234,252),(148,210,240),
                                                    (107,170,213),(72,149,176),(74,167,113),
                                                    (114,190,75),(203,217,88),(249,201,80),
                                                    (245,137,56),(234,84,43),(217,45,43),
                                                    (188,28,32),(156,22,27),(147,32,205)])/255.)
            cmap = mpl.colors.LinearSegmentedColormap.from_list('custom',colors,N=len(colors))
            norm = lambda xmin,xmax: mpl.colors.BoundaryNorm(boundaries=bounds,ncolors=len(bounds))
        else:
            norm = lambda xmin,xmax: mpl.colors.Normalize(xmin,xmax) 

    # only create the colorbar if requested
    cb_png_data = None
    if wisdom['colorbar'] is not None:
        if scale == 'original':
            logging.warning('postprocessor: Colorbar %s %s specified with scaling original' % (wisdom['name'],wisdom['colorbar']))
            logging.warning('postprocessor: Colorbar is not updated by the online visualization system correctly, use only with explicit scaling')
        cb_unit = wisdom['colorbar']
        cbu_min,cbu_max = convert_value(native_unit, cb_unit, fa_min), convert_value(native_unit, cb_unit, fa_max)
        #  colorbar + add it to the KMZ as a screen overlay
        logging.info('scalar_field_to_raster: making colorbar from %s to %s' % (cbu_min, cbu_max))
        spacing = wisdom.get('spacing', 'proportional')
        cb_png_data,levels = make_colorbar(
            [cbu_min, cbu_max], 'vertical', 2, cmap, wisdom['name'] + ' ' + cb_unit,
            ticks=ticks, spacing=spacing, norm=norm, ticklabels=ticklabels
        )
    else:
        levels = None

    # replace masked values by nans just in case
    fa.data[fa.mask]=np.nan
    fa.fill_value=np.nan

    # create the raster & get coordinate bounds
    raster_png_data,corner_coords = basemap_raster_mercator(lons, lats, fa, fa_min, fa_max, cmap, norm=norm)

    return raster_png_data, corner_coords, cb_png_data, levels


def vector_field_to_raster(u, v, lats, lons, wisdom):
    """
    Postprocess a vector field into arrows (used for wind).

    :param u: the vector in the U direction
    :param v: the vector in the V direction
    :param lats: the latitudes
    :param lons: the longitudes
    :param wisdom: a configuration dictionary controlling the visualization
    :return: the raster png as a StringIO, and the coordinates
    """

    if u.shape != lats.shape:
        raise PostprocError("U field does not correspond to grid size: field %s grid %s." % (u.shape, lats.shape))

    if v.shape != lats.shape:
        raise PostprocError("V field does not correspond to grid size: field %s grid %s." % (v.shape, lats.shape))

    # look at mins and maxes
    fa_min,fa_max = min(np.nanmin(u), np.nanmin(v)),max(np.nanmax(u), np.nanmax(v))

    # determine if we will use the range in the variable or a fixed range
    scale = wisdom['scale']
    if scale != 'original':
        fa_min, fa_max = scale[0], scale[1]
        u[u < fa_min] = fa_min
        u[u > fa_max] = fa_max
        v[v < fa_min] = fa_min
        v[v > fa_max] = fa_max

    # create the raster & get coordinate bounds, HACK to get better quiver resolution
    s = wisdom.get('ref', 2)
    raster_png_data,corner_coords = basemap_barbs_mercator(u[::s,::s], v[::s,::s], lats[::s,::s], lons[::s,::s])

    return raster_png_data, corner_coords

def scatter_to_raster(fa, lats, lons, wisdom):
    """
    Render a scatter variable into a geolocated raster and colorbar.

    :param fa: the array of values
    :param lats: the latitudes
    :param lons: the longitudes
    :param wisdom: a configuration dictionary controlling the visualization
    :return: a tuple with the raster as StringIO, its geolocation and the PNG colorbar StringIO (None if not requested)
    """
    # gather wisdom about the variable
    scale = wisdom['scale']
    native_unit = wisdom['native_unit']
    cmap_name = wisdom['colormap']
    cmap = mpl.cm.get_cmap(cmap_name)
    
    if lats.shape != fa.shape:
        raise PostprocError("Variable size %s does not correspond to grid size %s." % (fa.shape, lats.shape))
    
    if len(fa):
        logging.info('scatter_to_raster: variable %s min %s max %s' % (wisdom['name'], np.nanmin(fa), np.nanmax(fa)))
        
        # mask 'transparent' color value
        if 'transparent_values' in wisdom:
            rng = wisdom['transparent_values']
            fa = np.ma.masked_array(fa, np.logical_and(fa >= rng[0], fa <= rng[1]))
            logging.info('scatter_to_raster: transparent from %s to %, elements %s not masked %s' % (rng[0], rng[1], fa.size , fa.count()))
            logging.info('scatter_to_raster: masked variable %s min %s max %s' % (wisdom['name'], np.nanmin(fa), np.nanmax(fa)))
        else:
            fa = np.ma.masked_array(fa)
    
        # look at mins and maxes, transparent don't count
        if fa.count():
            fa_min,fa_max = np.nanmin(fa),np.nanmax(fa)
        else:
            fa_min, fa_max = 0.0, 0.0
        
        # determine if we will use the range in the variable or a fixed range
        if scale != 'original':
            m = fa.mask.copy()
            fa_min, fa_max = scale[0], scale[1]
            fa[fa < fa_min] = fa_min
            fa[fa > fa_max] = fa_max
            fa.mask = m
    else:
        fa = np.ma.masked_array(fa)
        # determine if we will use the range in the variable or a fixed range
        if scale != 'original':
            fa_min, fa_max = scale[0], scale[1]
        else:
            fa_min, fa_max = 0.0, 0.0

    # define distribution of colormap if necessary 
    norm = None
    ticks = None
    ticklabels = None
    if 'norm_opt' in wisdom: 
        norm_opt = wisdom['norm_opt']
        if norm_opt == 'lognorm':
            linthresh = wisdom.get('linthresh', .01)
            linscale = wisdom.get('linscale', .001)
            norm = lambda xmin,xmax: mpl.colors.SymLogNorm(linthresh=linthresh,linscale=linscale,vmin=xmin,vmax=xmax,base=10) 
            logging.info(f'scatter_to_raster: norm_opt={norm_opt} linthresh={linthresh} linscale={linscale}')
        elif norm_opt == 'boundary':
            bounds = wisdom.get('bounds', [0,1,2,4,6,8,12,16,20,25,30,40,60,100,200])
            ticks = bounds[1:]
            if 'labels' in wisdom and wisdom['labels'] is not None:
                ticklabels = wisdom['labels']
            elif wisdom['colorbar'] is not None:
                cb_unit = wisdom['colorbar']
                ticklabels = [convert_value(native_unit, cb_unit, b) for b in ticks]
            colors = wisdom.get('colors',np.array([(255,255,255),(197,234,252),(148,210,240),
                                                    (107,170,213),(72,149,176),(74,167,113),
                                                    (114,190,75),(203,217,88),(249,201,80),
                                                    (245,137,56),(234,84,43),(217,45,43),
                                                    (188,28,32),(156,22,27),(147,32,205)])/255.)
            cmap = mpl.colors.LinearSegmentedColormap.from_list('custom',colors,N=len(colors))
            bounds = bounds + [1e10]
            norm = lambda xmin,xmax,b=bounds: mpl.colors.BoundaryNorm(boundaries=b, ncolors=len(b))
            logging.info(f'scatter_to_raster: norm_opt={norm_opt} bounds={bounds} ticks={ticks} labels={ticklabels}')
        else:
            norm = lambda xmin,xmax: mpl.colors.Normalize(xmin,xmax) 
            logging.info(f'scatter_to_raster: norm_opt={norm_opt}')

    # only create the colorbar if requested
    cb_png_data = None
    if wisdom['colorbar'] is not None:
        if scale == 'original':
            logging.warning('postprocessor: Colorbar %s %s specified with scaling original' % (wisdom['name'],wisdom['colorbar']))
            logging.warning('postprocessor: Colorbar is not updated by the online visualization system correctly, use only with explicit scaling')
        cb_unit = wisdom['colorbar']
        cbu_min,cbu_max = convert_value(native_unit, cb_unit, fa_min), convert_value(native_unit, cb_unit, fa_max)
        #  colorbar + add it to the KMZ as a screen overlay
        logging.info('scatter_to_raster: making colorbar from %s to %s' % (cbu_min, cbu_max))
        spacing = wisdom.get('spacing', 'proportional') 
        cb_png_data,levels = make_colorbar(
            [cbu_min, cbu_max], 'vertical', 2, cmap, wisdom['name'] + ' ' + cb_unit,
            ticks=ticks, spacing=spacing, norm=norm, ticklabels=ticklabels
        )
    else:
        levels = None

    # replace masked values by nans just in case
    fa.data[fa.mask]=np.nan
    fa.fill_value=np.nan
    
    # other parameters
    if len(lons):
        bounds = wisdom.get('bbox', (np.amin(lons), np.amax(lons), np.amin(lats), np.amax(lats)))
    else:
        bounds = wisdom.get('bbox', (0, 0, 0, 0))
    alpha = wisdom.get('alpha', 0.7)
    marker = wisdom.get('marker', 'o')
    text = wisdom.get('text', False)
    size = wisdom.get('size', 15)
    linewidth = wisdom.get('linewidth', 0.7)
    logging.info('scatter_to_raster: options: alpha={}, marker={}, text={}, size={}, linewidth={}'.format(alpha,marker,text,size,linewidth))

    # create the raster & get coordinate bounds
    raster_png_data,corner_coords = basemap_scatter_mercator(
        [fa], [lons], [lats], bounds, [alpha], fa_min, fa_max, cmap, 
        size=size, marker=marker, linewidths=linewidth, text=text, norm=norm
    )

    return raster_png_data, corner_coords, cb_png_data, levels


class Postprocessor(object):
    """
    Postprocessing of WRF data.
    """

    def __init__(self, output_path, prod_name, wisdom_update = {}):
        """
        Initialize postprocessor with output parameters.

        :param output_path: path where postprocessing files are stored
        :param prod_name: name of manifest json file and prefix of all output files
        :param wisdom_update: an optional dictionary that maps variables to
                              modifications requested in their visualization wisdom
        """
        logging.info("Postprocessor: output_path=%s prod_name=%s" % (output_path, prod_name))
        dump(wisdom_update,"Postprocessor: wisdom_update")
        self.output_path = output_path
        self.product_name = prod_name
        self.manifest = {}
        self.wisdom_update = wisdom_update

        # in case the manifest exists, load the existing version
        mf_path = os.path.join(output_path, prod_name + '.json')
        if osp.exists(mf_path):
            self.manifest = json.load(open(mf_path))
            logging.info('postprocessor: Loaded manifest at %s' % mf_path)
            # dump(self.manifest,"postprocessor: manifest")
        else:
            logging.info('postprocessor: manifest at %s does not exist yet' % mf_path)


    def _scalar2raster(self, d, var, tndx, **tif_args):
        """
        Convert a single variable into a raster and colorbar.

        :param d: the netcdf file
        :param var: the variable name
        :param tndx: the time index
        :return: two StringIO objects, first is the PNG raster, second is PNG colorbar
        """
        # gather wisdom about the variable
        wisdom = get_wisdom(var).copy()
        wisdom.update(self.wisdom_update.get(var, {}))
        native_unit = wisdom['native_unit']
        cmap_name = wisdom['colormap']
        cmap = mpl.cm.get_cmap(cmap_name)

        # extract variable
        fa = wisdom['retrieve_as'](d,tndx) # this calls a lambda defined to read the required 2d field
        lat, lon = wisdom['grid'](d)
        if is_fire_var(var):
            fm,fn = strip_end(d)
            fa = fa[:fm,:fn]
            lon = lon[:fm,:fn] 
            lat = lat[:fm,:fn]

        if lat.shape != fa.shape:
            raise PostprocError("Variable %s size does not correspond to grid size." % var)

        # check for 'transparent' color value and mask
        if 'transparent_values' in wisdom:
            rng = wisdom['transparent_values']
            logging.info('_scalar2raster: variable %s min %s max %s masking transparent from %s to %s'
                % (var, np.nanmin(fa),np.nanmax(fa), rng[0], rng[1]))
            fa = np.ma.masked_array(fa, np.logical_and(fa >= rng[0], fa <= rng[1]))
        else:
            fa=np.ma.masked_array(fa)

        # look at mins and maxes
        fa_min,fa_max = np.nanmin(fa),np.nanmax(fa)

        # determine if we will use the range in the variable or a fixed range
        scale = wisdom['scale']
        if scale != 'original':
            m = fa.mask.copy()             # save mask, resetting values below will destroy the mask
            fa_min, fa_max = scale[0], scale[1]
            fa[fa < fa_min] = fa_min
            fa[fa > fa_max] = fa_max
            fa.mask = m                    # restore the mask

        # define distribution of colormap if necessary 
        norm = None
        ticks = None
        ticklabels = None
        if 'norm_opt' in wisdom: 
            norm_opt = wisdom['norm_opt']
            if norm_opt == 'lognorm':
                linthresh = wisdom.get('linthresh',.01)
                linscale = wisdom.get('linscale',.001)
                norm = lambda xmin,xmax: mpl.colors.SymLogNorm(linthresh=linthresh,linscale=linscale,vmin=xmin,vmax=xmax,base=10) 
            elif norm_opt == 'boundary':
                bounds = wisdom.get('bounds',[0,1,2,4,6,8,12,16,20,25,30,40,60,100,200])
                if 'labels' in wisdom and wisdom['labels'] is not None:
                    ticklabels = wisdom['labels']
                elif wisdom['colorbar'] is not None:
                    cb_unit = wisdom['colorbar']
                    ticklabels = [convert_value(native_unit, cb_unit, b) for b in bounds[1:]]
                ticks = bounds[1:]
                bounds = bounds + [1e10]
                colors = wisdom.get('colors',np.array([(255,255,255),(197,234,252),(148,210,240),
                                                       (107,170,213),(72,149,176),(74,167,113),
                                                       (114,190,75),(203,217,88),(249,201,80),
                                                       (245,137,56),(234,84,43),(217,45,43),
                                                       (188,28,32),(156,22,27),(147,32,205)])/255.)
                cmap = mpl.colors.LinearSegmentedColormap.from_list('custom',colors,N=len(colors))
                norm = lambda xmin,xmax: mpl.colors.BoundaryNorm(boundaries=bounds,ncolors=len(bounds))
            else:
                norm = lambda xmin,xmax: mpl.colors.Normalize(xmin,xmax) 

        # only create the colorbar if requested
        cb_png_data = None
        if wisdom['colorbar'] is not None:
            cb_unit = wisdom['colorbar']
            cbu_min,cbu_max = convert_value(native_unit, cb_unit, fa_min), convert_value(native_unit, cb_unit, fa_max)
            #  colorbar + add it to the KMZ as a screen overlay
            legend = wisdom['name'] + ' ' + cb_unit
            logging.info('_scalar2raster: variable %s colorbar from %s to %s %s' % (var, cbu_min,cbu_max, legend))
            spacing = wisdom.get('spacing','proportional')
            cb_png_data,levels = make_colorbar([cbu_min, cbu_max],'vertical',2,cmap,legend,ticks=ticks,spacing=spacing,norm=norm,ticklabels=ticklabels)
        else:
            levels = None

        # replace masked values by nans just in case
        fa.data[fa.mask]=np.nan
        fa.fill_value=np.nan

        logging.info('_scalar2raster: variable %s elements %s count %s not masked %s min %s max %s'
            % (var, fa.size , fa.count(), np.count_nonzero(fa.mask == False), np.nanmin(fa),np.nanmax(fa) ))
        
        # create the raster & get coordinate bounds
        raster_png_data,corner_coords = basemap_raster_mercator(lon,lat,fa,fa_min,fa_max,cmap,norm=norm)

        # create GeoTIFF file
        if tif_args:
            if 'ndv' in tif_args.keys():
                ndv = tif_args['ndv']
            else:
                ndv = -9999.0
            fa[np.isnan(fa)] = ndv
            tif_path = tif_args.get('tif_path')
            logging.info('_scalar2raster: writting GeoTIFF file {}'.format(tif_path))
            crs = tif_args.get('crs')
            geot = tif_args.get('geot')
            GeoDriver.from_elements(fa, crs, geot).to_geotiff(tif_path, desc = wisdom['name'], unit = native_unit, ndv = ndv)

        return raster_png_data, corner_coords, cb_png_data, levels


    def _vector2raster(self, d, var, tndx, **tif_args):
        """
        Postprocess a vector field into barbs (used for wind) and contains hacks for WRF staggered grids.

        :param d: the open netCDF file
        :param var: the name of the variable in var_wisdom
        :param tndx: the time index
        :return: the raster png as a StringIO, and the coordinates
        """
        # gather wisdom about the variable
        wisdom = get_wisdom(var).copy()
        wisdom.update(self.wisdom_update.get(var, {}))
        native_unit = wisdom['native_unit']
        cmap_name = wisdom.get('colormap')
        u_name, v_name = wisdom['components']
        lat, lon = wisdom['grid'](d)

        # extract variable
        uw, vw = get_wisdom(u_name), get_wisdom(v_name)
        uw.update(self.wisdom_update.get(u_name, {}))
        vw.update(self.wisdom_update.get(v_name, {}))
        u = uw['retrieve_as'](d, tndx)
        v = vw['retrieve_as'](d, tndx)
        if is_fire_var(var):
            fm,fn = strip_end(d)
            u = u[:fm,:fn]
            v = v[:fm,:fn]
            lon = lon[:fm,:fn] 
            lat = lat[:fm,:fn]

        if u.shape != lat.shape:
            raise PostprocError("Variable %s size does not correspond to grid size: var %s grid %s." % (u_name, u.shape, lat.shape))

        if v.shape != lat.shape:
            raise PostprocError("Variable %s size does not correspond to grid size." % v_name)

        # look at mins and maxes
        fa_min,fa_max = min(np.nanmin(u), np.nanmin(v)),max(np.nanmax(u), np.nanmax(v))

        # determine if we will use the range in the variable or a fixed range
        scale = wisdom['scale']
        if scale != 'original':
            fa_min, fa_max = scale[0], scale[1]
            u[u < fa_min] = fa_min
            u[u > fa_max] = fa_max
            v[v < fa_min] = fa_min
            v[v > fa_max] = fa_max

        # HACK to get better quiver resolution
        s = 3
        cb_png_data = None
        levels = None
        if cmap_name is not None:
            # create speed array and calculat min and max
            speed_formula = wisdom.get('speed_formula',lambda u,v: np.sqrt(u**2+v**2))
            fa = speed_formula(u,v)
            fa_min,fa_max = np.nanmin(fa),np.nanmax(fa)
            if 'speed_scale' in wisdom: 
                speed_scale = wisdom['speed_scale']
                fa_min, fa_max = speed_scale[0], speed_scale[1]
                m = fa.mask.copy()             # save mask, resetting values below will destroy the mask
                fa[fa < fa_min] = fa_min
                fa[fa > fa_max] = fa_max
                fa.mask = m
            # replace masked values by nans just in case
            fa.data[fa.mask]=np.nan
            fa.fill_value=np.nan
            logging.info('_vector2scalar: variable %s elements %s count %s not masked %s min %s max %s scale min %s max %s'
                        % (var, fa.size , fa.count(), np.count_nonzero(fa.mask == False), np.nanmin(fa), np.nanmax(fa), fa_min, fa_max))
            # create cmap
            cmap = mpl.cm.get_cmap(cmap_name)
            # check for 'transparent' color value and mask
            if 'transparent_values' in wisdom:
                rng = wisdom['transparent_values']
                fa = np.ma.masked_array(fa, np.logical_and(fa >= rng[0], fa <= rng[1]))
            else:
                fa=np.ma.masked_array(fa)
            # define distribution of colormap if necessary
            norm = None
            ticks = None
            ticklabels = None
            if 'norm_opt' in wisdom:
                norm_opt = wisdom['norm_opt']
                if norm_opt == 'lognorm':
                    linthresh = wisdom.get('linthresh',.01)
                    linscale = wisdom.get('linscale',.001)
                    norm = lambda xmin,xmax: mpl.colors.SymLogNorm(linthresh=linthresh,linscale=linscale,vmin=xmin,vmax=xmax,base=10)
                elif norm_opt == 'boundary':
                    bounds = wisdom.get('bounds',[0,1,2,4,6,8,12,16,20,25,30,40,60,100,200])
                    if 'labels' in wisdom and wisdom['labels'] is not None:
                        ticklabels = wisdom['labels']
                    elif wisdom['colorbar'] is not None:
                        cb_unit = wisdom['colorbar']
                        ticklabels = [convert_value(native_unit, cb_unit, b) for b in bounds[1:]]
                    ticks = bounds[1:]
                    bounds = bounds + [1e10]
                    colors = wisdom.get('colors',np.array([(255,255,255),(197,234,252),(148,210,240),
                                                       (107,170,213),(72,149,176),(74,167,113),
                                                       (114,190,75),(203,217,88),(249,201,80),
                                                       (245,137,56),(234,84,43),(217,45,43),
                                                       (188,28,32),(156,22,27),(147,32,205)])/255.)
                    cmap = mpl.colors.LinearSegmentedColormap.from_list('custom',colors,N=len(colors))
                    norm = lambda xmin,xmax: mpl.colors.BoundaryNorm(boundaries=bounds,ncolors=len(bounds))
                else:
                    norm = lambda xmin,xmax: mpl.colors.Normalize(xmin,xmax)
            # only create the colorbar if requested
            if wisdom['colorbar'] is not None:
                cb_unit = wisdom['colorbar']
                cbu_min,cbu_max = convert_value(native_unit, cb_unit, fa_min), convert_value(native_unit, cb_unit, fa_max)
                #  colorbar + add it to the KMZ as a screen overlay
                legend = wisdom['name'] + ' ' + cb_unit
                logging.info('_vector2raster: variable %s colorbar from %s to %s %s' % (var, cbu_min, cbu_max, legend))
                spacing = wisdom.get('spacing','proportional')
                cb_png_data,levels = make_colorbar([cbu_min, cbu_max],'vertical',2,cmap,legend,ticks=ticks,spacing=spacing,norm=norm,ticklabels=ticklabels)
            raster_png_data,corner_coords = basemap_barbs_mercator(u[::s,::s],v[::s,::s],lat[::s,::s],lon[::s,::s],fa[::s,::s],fa_min,fa_max,cmap,norm=norm)
        else:
            # create the raster & get coordinate bounds
            raster_png_data,corner_coords = basemap_barbs_mercator(u[::s,::s],v[::s,::s],lat[::s,::s],lon[::s,::s])

        # create GeoTIFF file
        if tif_args:
            if 'ndv' in tif_args.keys():
                ndv = tif_args['ndv']
            else:
                ndv = -9999.0
            u[np.isnan(u)] = ndv
            v[np.isnan(v)] = ndv
            tif_path = tif_args.get('tif_path')
            logging.info('_vector2raster: writting GeoTIFF file {}'.format(tif_path))
            crs = tif_args.get('crs')
            gt = tif_args.get('geot')
            GeoDriver.from_elements(np.array([u, v]), crs, gt).to_geotiff(tif_path, desc = [uw['name'], vw['name']], unit = native_unit, ndv = ndv)

        return raster_png_data, corner_coords, cb_png_data, levels


    def _sat2raster(self, dgs, dfs, gts, sat, bounds):
        """
        Postprocess a satellite granule into raster.

        :param dgs: list of open geolocation files
        :param dfs: list of open fire files
        :param gts: list of granule times
        :param sat: the name of the satellite variable in var_wisdom
        :param bounds: bounds for the bounding box satellite data
        :return: the raster png as a StringIO, and the coordinates
        """
        # gather wisdom about the satellite variable
        wisdom = get_wisdom(sat).copy()
        wisdom.update(self.wisdom_update.get(sat, {}))
        cmap_name = wisdom['colormap']

        values = wisdom['options'].get('values',(3,5,7,8,9))
        alphas = wisdom['options'].get('alphas',(.5,.5,.6,.7,.8))
        labels = wisdom['options'].get('labels',('Water','Ground','Fire low','Fire nominal','Fire high'))
        colors = wisdom['options'].get('colors',((0,0,.5),(0,.5,0),(1,1,0),(1,.65,0),(.5,0,0)))

        # for each pair of files
        N = len(values)
        lons,lats,vals = [[]]*N,[[]]*N,[[]]*N
        for dg,df in zip(dgs,dfs):
            # extract variables
            lat, lon = wisdom['grid'](dg)
            fa = wisdom['retrieve_as'](df)

            # check variables
            if fa.shape != lat.shape:
                raise PostprocError("Variable %s size does not correspond to grid size." % sat)

            # process variables
            flon = np.array(lon).ravel()
            flat = np.array(lat).ravel()
            mask = np.array(fa).ravel()
            bbox = np.logical_and(np.logical_and(np.logical_and(flon>bounds[0],flon<bounds[1]),flat>bounds[2]),flat<bounds[3])
            categories = [np.array(mask[bbox] == value) for value in values]

            for i,cat in enumerate(categories):
                lons[i] = np.append(lons[i],flon[bbox][cat])
                lats[i] = np.append(lats[i],flat[bbox][cat])
                vals[i] = np.append(vals[i],np.ones(cat.sum())*i)
                logging.info('_sat2raster: variable %s, category %s, shape of data (%s,%s,%s)' % (sat,i,lons[i].shape,lats[i].shape,vals[i].shape))

        # create discrete colormap
        if cmap_name == 'discrete':
            cmap = mpl.colors.LinearSegmentedColormap.from_list('cb_'+cmap_name, colors, N=N)
        else:
            cmap = mpl.cm.get_cmap(cmap_name)

        # only create the colorbar if requested
        cb_png_data = None
        if wisdom['colorbar'] is not None:
            if cmap_name == 'discrete':
                cb_N = 5
                cb_labels = ('Water','Ground','Fire low','Fire nominal','Fire high')
                cb_colors = ((0,0,.5),(0,.5,0),(1,1,0),(1,.65,0),(.5,0,0))
                cb_cmap = mpl.colors.LinearSegmentedColormap.from_list('cb_'+cmap_name, cb_colors, N=cb_N)
                #  colorbar + add it to the KMZ as a screen overlay
                legend = ''
                logging.info('_sat2raster: variable %s colorbar from %s to %s %s' % (sat, 0, cb_N, legend))
                cb_png_data = make_discrete_colorbar(cb_labels, cb_colors, 'vertical', 2, cb_cmap, legend)

        # create the raster & get coordinate bounds
        cmin = -.5
        cmax = N-.5
        raster_png_data,corner_coords = basemap_scatter_mercator(vals,lons,lats,bounds,alphas,cmin,cmax,cmap)

        return raster_png_data, corner_coords, cb_png_data, np.array(values).tolist()


    def _sat2raster_empty(self, sat, bounds):
        """
        Postprocess an empty satellite granule into raster.

        :param sat: the name of the satellite variable in var_wisdom
        :param bounds: bounds for the bounding box satellite data
        :return: the raster png as a StringIO, and the coordinates
        """
        # gather wisdom about the satellite variable
        wisdom = get_wisdom(sat).copy()
        wisdom.update(self.wisdom_update.get(sat, {}))
        cmap_name = wisdom['colormap']

        values = wisdom['options'].get('values',(3,5,7,8,9))
        labels = wisdom['options'].get('labels',('Water','Ground','Fire low','Fire nominal','Fire high'))
        colors = wisdom['options'].get('colors',((0,0,.5),(0,.5,0),(1,1,0),(1,.65,0),(.5,0,0)))
        N = len(labels)

        # create discrete colormap
        if cmap_name == 'discrete':
            cmap = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=N)
        else:
            cmap = mpl.cm.get_cmap(cmap_name)

        # only create the colorbar if requested
        cb_png_data = None
        if wisdom['colorbar'] is not None:
            cb_N = 5
            cb_labels = ('Water','Ground','Fire low','Fire nominal','Fire high')
            cb_colors = ((0,0,.5),(0,.5,0),(1,1,0),(1,.65,0),(.5,0,0))
            cb_cmap = mpl.colors.LinearSegmentedColormap.from_list('cb_'+cmap_name, cb_colors, N=cb_N)
            #  colorbar + add it to the KMZ as a screen overlay
            legend = ''
            logging.info('_sat2raster: variable %s colorbar from %s to %s %s' % (sat, 0, cb_N, legend))
            cb_png_data = make_discrete_colorbar(cb_labels, cb_colors, 'vertical', 2, cb_cmap, legend)

        # create the raster & get coordinate bounds
        raster_png_data,corner_coords = basemap_scatter_mercator([],[],[],bounds,[],0,0,cmap)

        return raster_png_data, corner_coords, cb_png_data, np.array(values).tolist()


    def _scalar2png(self, d, var, tndx, out_path, **tif_args):
        """
        Postprocess a scalar field into a raster file and a colorbar file, both PNG.

        :param d: the open netCDF file
        :param var: the variable name
        :param tndx: the temporal index of the grid to store
        :param out_path: the path to the KMZ output
        :return: the path to the raster, to the colorbar and the bounding coordinates
        """
        # render the raster & colorbar
        raster_png_data, corner_coords, cb_png_data, levels = self._scalar2raster(d, var, tndx, **tif_args)

        # write raster file
        raster_path = out_path + "-raster.png"
        logging.info("writing file %s size %s" % (raster_path, sys.getsizeof(raster_png_data)))
        with open(raster_path, 'wb') as f:
            f.write(raster_png_data)

        # write colorbar file
        colorbar_path = None
        if cb_png_data is not None:
            colorbar_path = out_path + "-cb.png"
            with open(colorbar_path, "wb") as f:
                f.write(cb_png_data)

        return raster_path, colorbar_path, corner_coords, levels


    def _vector2png(self, d, var, tndx, out_path, **tif_args):
        """
        Postprocess a single vector variable ``var`` and stores result in a raster file.

        :param d: the open netCDF file
        :param var: the variable name
        :param tndx: the temporal index of the grid to store
        :param out_path: the path to the KMZ output
        :return: the path to the raster and the bounding coordinates
        """
        # render the raster & colorbar
        raster_png_data, corner_coords, cb_png_data, levels = self._vector2raster(d, var, tndx, **tif_args)

        raster_path = out_path + '-raster.png'
        with open(raster_path, 'wb') as f:
            f.write(raster_png_data)

        # write colorbar file
        colorbar_path = None
        if cb_png_data is not None:
            colorbar_path = out_path + "-cb.png"
            with open(colorbar_path, "wb") as f:
                f.write(cb_png_data)

        return raster_path, colorbar_path, corner_coords, levels


    def _sat2png(self, dgs, dfs, gts, sat, bounds, out_path):
        """
        Postprocess a single sat granule variable ``sat`` and stores result in a raster file.

        :param dgs: list of open geolocation files
        :param dfs: list of open fire files
        :param gts: list of granule times
        :param sat: the satellite variable name
        :param bounds: bounds for the bounding box satellite data
        :param out_path: the path to the KMZ output
        :return: the path to the raster and the bounding coordinates
        """
        # render the raster & colorbar
        raster_png_data, corner_coords, cb_png_data, levels = self._sat2raster(dgs, dfs, gts, sat, bounds)

        raster_path = out_path + '-raster.png'
        logging.info("writing file %s size %s" % (raster_path, sys.getsizeof(raster_png_data)))

        with open(raster_path, 'wb') as f:
                f.write(raster_png_data)

        # write colorbar file
        colorbar_path = None
        if cb_png_data is not None:
            colorbar_path = out_path + "-cb.png"
            with open(colorbar_path, "wb") as f:
                f.write(cb_png_data)

        return raster_path, colorbar_path, corner_coords, levels


    def _sat2png_empty(self, sat, bounds, out_path):
        """
        Postprocess an empty sat granule variable ``sat`` and stores result in an empty raster file.

        :param sat: the satellite variable name
        :param bounds: bounds for the bounding box satellite data
        :param out_path: the path to the KMZ output
        :return: the path to the raster and the bounding coordinates
        """
        # render the raster & colorbar
        raster_png_data, corner_coords, cb_png_data, levels = self._sat2raster_empty(sat, bounds)

        raster_path = out_path + '-raster.png'
        logging.info("writing file %s size %s" % (raster_path, sys.getsizeof(raster_png_data)))

        with open(raster_path, 'wb') as f:
                f.write(raster_png_data)

        # write colorbar file
        colorbar_path = None
        if cb_png_data is not None:
            colorbar_path = out_path + "-cb.png"
            with open(colorbar_path, "wb") as f:
                f.write(cb_png_data)

        return raster_path, colorbar_path, corner_coords, levels


    def _scalar2kmz(self, d, var, tndx, ts_esmf_begin, ts_esmf_end, out_path, cleanup = True, **tif_args):
        """
        Postprocess a single raster variable ``fa`` and store result in out_path.

        :param d: the open netCDF file
        :param var: the variable name
        :param tndx: the temporal index of the grid to store
        :param ts_esmf: time string yyyy-mm-dd_HH:MM:SS
        :param out_path: the path to the KMZ output
        :param cleanup: if true, delete png files
        :return: the path to the generated KMZ
        """
        # construct kml file

        name = ts_esmf_begin + ' ' + var
        file = kml.Kml(name = name)
        doc = file.newdocument(name = name)
        doc.timespan.begin=ts_esmf_begin.replace('_','T')+'Z'
        if ts_esmf_end is not None:
            doc.timespan.end=ts_esmf_end.replace('_','T')+'Z'

        # generate the png files
        raster_path, cb_path, corner_coords, levels = self._scalar2png(d, var, tndx, out_path, **tif_args)

        # add colorbar to KMZ
        if cb_path is not None:
            cbo = doc.newscreenoverlay(name='colorbar')
            cbo.overlayxy = kml.OverlayXY(x=0,y=1,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.screenxy = kml.ScreenXY(x=0.02,y=0.95,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.size = kml.Size(x=150,y=300,xunits=kml.Units.pixels,yunits=kml.Units.pixels)
            cbo.color = kml.Color.rgb(255,255,255,a=150)
            cbo.visibility = 1
            #doc.addfile(cb_path)
            cbo.icon.href=cb_path

        # add ground overlay
        ground = doc.newgroundoverlay(name=var,color='80ffffff')
        ground.gxlatlonquad.coords = corner_coords
        #doc.addfile(raster_path)
        ground.icon.href = raster_path

        # build output file
        kmz_path = out_path + ".kmz"
        file.savekmz(kmz_path)

        # cleanup
        if cleanup:
            os.remove(raster_path)
            if cb_path is not None:
                os.remove(cb_path)

        return kmz_path, raster_path, cb_path, corner_coords, levels


    def _vector2kmz(self, d, var, tndx, ts_esmf_begin, ts_esmf_end, out_path, cleanup = True, **tif_args):
        """
        Postprocess a single vector variable ``var`` and store result in out_path.

        :param d: the open netCDF file
        :param var: the variable name
        :param tndx: the temporal index of the grid to store
        :param ts_esmf: time string yyyy-mm-dd_HH:MM:SS
        :param out_path: the path to the KMZ output
        :param cleanup: if True, PNG files are deleted after KMZ is build
        :return: the path to the generated KMZ
        """
        # construct kml file

        name = ts_esmf_begin + ' ' + var
        file = kml.Kml(name = name)

        doc = file.newdocument(name = name)
        doc.timespan.begin=ts_esmf_begin.replace('_','T')+'Z'
        if ts_esmf_end is not None:
            doc.timespan.end=ts_esmf_end.replace('_','T')+'Z'

        # generate the png files
        raster_path, cb_path, corner_coords, levels = self._vector2png(d, var, tndx, out_path, **tif_args)

        # add colorbar to KMZ
        if cb_path is not None:
            cbo = doc.newscreenoverlay(name='colorbar')
            cbo.overlayxy = kml.OverlayXY(x=0,y=1,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.screenxy = kml.ScreenXY(x=0.02,y=0.95,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.size = kml.Size(x=150,y=300,xunits=kml.Units.pixels,yunits=kml.Units.pixels)
            cbo.color = kml.Color.rgb(255,255,255,a=150)
            cbo.visibility = 1
            #doc.addfile(cb_path)
            cbo.icon.href=cb_path

        # add ground overlay
        ground = doc.newgroundoverlay(name=var,color='80ffffff')
        ground.gxlatlonquad.coords = corner_coords
        #doc.addfile(raster_path)
        ground.icon.href = raster_path

        # build output file
        kmz_path = out_path + ".kmz"
        file.savekmz(kmz_path)

        # cleanup
        if cleanup:
            os.remove(raster_path)

        return kmz_path, raster_path, cb_path, corner_coords, levels


    def _sat2kmz(self, dgs, dfs, gts, sat, ts_esmf_begin, ts_esmf_end, bounds, out_path, cleanup = True):
        """
        Postprocess a single satellite variable ``sat`` and store result in out_path.

        :param dgs: list of open geolocation files
        :param dfs: list of open fire files
        :param gts: list of granule times
        :param sat: the sat name
        :param ts_esmf_begin: time string yyyy-mm-dd_HH:MM:SS
        :param ts_esmf_end: time string yyyy-mm-dd_HH:MM:SS
        :param bounds: bounds for the bounding box satellite data
        :param out_path: the path to the KMZ output
        :param cleanup: if True, PNG files are deleted after KMZ is build
        :return: the path to the generated KMZ
        """
        # construct kml file

        name = ts_esmf_begin + ' ' + sat
        file = kml.Kml(name = name)
        doc = file.newdocument(name = name)
        doc.timespan.begin=ts_esmf_begin.replace('_','T')+'Z'
        if ts_esmf_end is not None:
            doc.timespan.end=ts_esmf_end.replace('_','T')+'Z'

        # generate the png files
        raster_path, cb_path, corner_coords, levels = self._sat2png(dgs, dfs, gts, sat, bounds, out_path)

        # add colorbar to KMZ
        if cb_path is not None:
            cbo = doc.newscreenoverlay(name='colorbar')
            cbo.overlayxy = kml.OverlayXY(x=0,y=1,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.screenxy = kml.ScreenXY(x=0.02,y=0.95,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.size = kml.Size(x=150,y=300,xunits=kml.Units.pixels,yunits=kml.Units.pixels)
            cbo.color = kml.Color.rgb(255,255,255,a=150)
            cbo.visibility = 1
            cbo.icon.href=cb_path

        # add ground overlay
        ground = doc.newgroundoverlay(name=sat,color='80ffffff')
        ground.gxlatlonquad.coords = corner_coords
        ground.icon.href = raster_path

        # build output file
        kmz_path = out_path + ".kmz"
        file.savekmz(kmz_path)

        # cleanup
        if cleanup:
            os.remove(raster_path)
            if cb_path is not None:
                os.remove(cb_path)

        return kmz_path, raster_path, cb_path, corner_coords, levels


    def _sat2kmz_empty(self, sat, ts_esmf_begin, ts_esmf_end, bounds, out_path, cleanup = True):
        """
        Postprocess an empty satellite variable ``sat`` and store result in out_path.

        :param sat: the sat name
        :param ts_esmf_begin: time string yyyy-mm-dd_HH:MM:SS
        :param ts_esmf_end: time string yyyy-mm-dd_HH:MM:SS
        :param bounds: bounds for the bounding box satellite data
        :param out_path: the path to the KMZ output
        :param cleanup: if True, PNG files are deleted after KMZ is build
        :return: the path to the generated KMZ
        """
        # construct kml file

        name = ts_esmf_begin + ' ' + sat
        file = kml.Kml(name = name)
        doc = file.newdocument(name = name)
        doc.timespan.begin=ts_esmf_begin.replace('_','T')+'Z'
        if ts_esmf_end is not None:
            doc.timespan.end=ts_esmf_end.replace('_','T')+'Z'

        # generate the png files
        raster_path, cb_path, corner_coords, levels = self._sat2png_empty(sat, bounds, out_path)

        # add colorbar to KMZ
        if cb_path is not None:
            cbo = doc.newscreenoverlay(name='colorbar')
            cbo.overlayxy = kml.OverlayXY(x=0,y=1,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.screenxy = kml.ScreenXY(x=0.02,y=0.95,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
            cbo.size = kml.Size(x=150,y=300,xunits=kml.Units.pixels,yunits=kml.Units.pixels)
            cbo.color = kml.Color.rgb(255,255,255,a=150)
            cbo.visibility = 1
            cbo.icon.href=cb_path

        # add ground overlay
        ground = doc.newgroundoverlay(name=sat,color='80ffffff')
        ground.gxlatlonquad.coords = corner_coords
        ground.icon.href = raster_path

        # build output file
        kmz_path = out_path + ".kmz"
        file.savekmz(kmz_path)

        # cleanup
        if cleanup:
            os.remove(raster_path)
            if cb_path is not None:
                os.remove(cb_path)

        return kmz_path, raster_path, cb_path, corner_coords, levels


    def open_file(self, path_file):
        """
        Open file depending on its extension

        :param path_file: local path to the file
        """
        if not osp.exists(path_file):
            logging.error('open_file: file %s does not exist locally' % path_file)
            raise PostprocError('failed opening file {}'.format(path_file)) 
        ext = osp.splitext(path_file)[1]
        logging.info('open_file: open file %s with extension %s' % (path_file, ext))
        if ext == ".nc":
            try:
                d = nc4.Dataset(str(path_file),'r')
            except Exception as e:
                logging.error('open_file: can not open file %s with exception %s' % (path_file,e))
                raise PostprocError('failed opening file {}'.format(path_file)) 
        elif ext == ".hdf":
            try:
                d = SD(str(path_file),SDC.READ)
            except Exception as e:
                logging.error('open_file: can not open file %s with exception %s' % (path_file,e))
                raise PostprocError('failed opening file {}'.format(path_file)) 
        elif ext == ".h5":
            try:
                d = h5py.File(str(path_file),'r')
            except Exception as e:
                logging.error('open_file: can not open file %s with exception %s' % (path_file,e))
                raise PostprocError('failed opening file {}'.format(path_file)) 
        else:
            logging.error('open_file: unrecognized extension %s' % ext)
            raise PostprocError('failed opening file {}'.format(path_file)) 
        return d,ext


    def close_file(self, d, ext):
        """
        Close file depending on its extension

        :param d: open file to close
        :param ext: extension of the file
        """
        if ext == ".nc":
            d.close()
        elif ext == ".hdf":
            d.end()
        elif ext == ".h5":
            d.close()
        else:
            logging.error('close_file: unrecognized extension %s' % ext)


    def process_sats(self, jsat, dom_id, ts_esmf, sats):
        """
        Postprocess satellite data at a given simulation interval into PNG and KMZ files.

        :param jsat: dictionary of satellite products acquiered
        :param dom_id: the domain identifier
            :param ts_esmf: time stamp in ESMF format
        :param sats: list of sat variables to process
        """
        traceargs()

        dt_before = timedelta(minutes = jsat.sat_interval[str(dom_id)][0])
        dt_after = timedelta(minutes = jsat.sat_interval[str(dom_id)][1])
        dt = timedelta(minutes = jsat.dt[str(dom_id)])
        if dt_after + dt_before < dt:
            dt_after = dt
            dt_before = timedelta(minutes=0)
        ts_initial = (esmf_to_utc(ts_esmf) - dt_before).replace(second=0, microsecond=0)
        ts_final = (esmf_to_utc(ts_esmf) + dt_after).replace(second=0, microsecond=0)
        logging.info('process_sats: looking from time %s to time %s' % (ts_initial,ts_final))
        bounds = jsat.bounds[str(dom_id)]
        logging.info('process_sats: bounding box %s' % bounds)

        for sat in sats:
            logging.info('process_sats: postprocessing %s for time %s' % (sat, ts_esmf))
            sat_source = jsat['satprod_satsource'][sat]
            logging.info('process_sats: product %s from %s source' % (sat,sat_source))
            dgs,dfs,egs,efs,gts = [],[],[],[],[]
            for k,gran in jsat.granules[sat_source].items():
                gran_time = esmf_to_utc(gran['time_start_iso'])
                logging.debug('process_sats: evaluating product %s, granule %s, at time %s, and for time %s' % (sat, k, utc_to_esmf(gran_time), ts_esmf))
                if gran_time >= ts_initial and gran_time < ts_final:
                    try:
                        dg,eg = self.open_file(gran['geo_local_path'])
                        df,ef = self.open_file(gran['fire_local_path'])
                        dgs.append(dg)
                        dfs.append(df)
                        egs.append(eg)
                        efs.append(ef)
                        gts.append(gran_time)
                    except Exception as e:
                        logging.warning("Exception %s while evaluating granule %s from product %s for time %s" % (e, gran, sat, ts_esmf))
                        logging.warning(traceback.print_exc())
            if not dgs:
                logging.info('process_sats: any granule %s in output process interval %s - %s' % (sat, utc_to_esmf(ts_initial), utc_to_esmf(ts_final)))
                try:
                    outpath_base = osp.join(self.output_path, self.product_name + ("-%02d-" % dom_id) + "sat_empty")
                    kmz_path, raster_path, cb_path, coords, levels, mf_upd = None, None, None, None, None, {}
                    if osp.exists(outpath_base+".kmz"):
                        logging.info('process_sats: empty sat %s already processed' % (outpath_base+".kmz"))
                        kmz_path = outpath_base+".kmz"
                        raster_path = outpath_base+"-raster.png"
                        cb_path = outpath_base+"-cb.png"
                        numpy_bounds = [ (bounds[0],bounds[2]),
                                (bounds[1],bounds[2]),
                                (bounds[1],bounds[3]),
                                (bounds[0],bounds[3]) ]
                        float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
                        coords = float_bounds
                    else:
                        logging.info('process_sats: processing empty sat %s for the first time' % (outpath_base+".kmz"))
                        kmz_path,raster_path,cb_path,coords,levels = self._sat2kmz_empty(sat, ts_esmf, None, bounds, outpath_base, cleanup=False)
                    if levels is not None:
                        mf_upd['levels'] = levels
                    if cb_path is not None:
                        mf_upd['colorbar'] = osp.basename(cb_path)
                    mf_upd['kml'] = osp.basename(kmz_path)
                    mf_upd['raster'] = osp.basename(raster_path)
                    mf_upd['coords'] = coords
                    logging.info("updating manifest for variable %s at time %s with manifest %s" % (sat, ts_esmf, mf_upd))
                    self._update_manifest(dom_id, ts_esmf, sat, mf_upd)
                except Exception as e:
                    logging.warning("Exception %s while postprocessing %s for time %s" % (e, sat, ts_esmf))
                    logging.warning(traceback.print_exc())
            else:
                logging.info('process_sats: some granule %s is in output process interval %s - %s' % (sat, utc_to_esmf(ts_initial), utc_to_esmf(ts_final)))
                try:
                    outpath_base = osp.join(self.output_path, self.product_name + ("-%02d-" % dom_id) + ts_esmf + "-" + sat)
                    kmz_path, raster_path, cb_path, coords, mf_upd = None, None, None, None, {}
                    kmz_path,raster_path,cb_path,coords,levels = self._sat2kmz(dgs, dfs, gts, sat, ts_esmf, None, bounds, outpath_base, cleanup=False)
                    if levels is not None:
                        mf_upd['levels'] = levels
                    if cb_path is not None:
                        mf_upd['colorbar'] = osp.basename(cb_path)
                    mf_upd['kml'] = osp.basename(kmz_path)
                    mf_upd['raster'] = osp.basename(raster_path)
                    mf_upd['coords'] = coords
                    logging.info("updating manifest for variable %s at time %s with manifest %s" % (sat, ts_esmf, mf_upd))
                    self._update_manifest(dom_id, ts_esmf, sat, mf_upd)
                    for i in range(len(dgs)):
                        self.close_file(dgs[i],egs[i])
                        self.close_file(dfs[i],efs[i])
                except Exception as e:
                    logging.warning("Exception %s while postprocessing %s for time %s" % (e, sat, ts_esmf))
                    logging.warning(traceback.print_exc())


    def process_vars(self, wrfout_path, dom_id, ts_esmf, vars, tif_proc = False, tslist = None):
        """
        Postprocess a list of scalar or vector fields at a given simulation time into PNG and KMZ
        files.

        :param wrfout_path: WRF file to process
        :param dom_id: the domain identifier
        :param ts_esmf: time stamp in ESMF format
        :param vars: list of variables to process
        :param tif_proc: boolen if process tif files
        """
        traceargs()

        logging.info('process_vars: looking for time %s in %s' % (ts_esmf,wrfout_path))
        max_retries = 10
        for k in range(max_retries):
            try:
                # open the netCDF dataset
                with open('/dev/null','w') as f:
                    check_call(['ncdump','-h','%s'%wrfout_path],stdout=f,stderr=f)
                logging.info('process_vars: netCDF file checked, using python netCDF4')
                d = nc4.Dataset(wrfout_path,'r')
                # extract ESMF string times and identify timestamp of interest
                times = [''.join(x) for x in d.variables['Times'][:].astype(str)]
                logging.info('process_vars: time steps found %s' % str(times))
                # make sure time step is processed on file
                if ts_esmf in times:
                    logging.info('process_vars: time step %s found in wrfout %s in retry %s' % (ts_esmf,wrfout_path,str(k+1)))
                    tndx = times.index(ts_esmf)
                    break
                else:
                    if k == max_retries-1:
                        logging.error('process_vars: cannot find time %s in %s' % (ts_esmf,wrfout_path))
                        logging.info('process_vars: Available times: %s' % times)
                        raise PostprocError("process_vars: Time %s not in %s" % (ts_esmf,osp.basename(wrfout_path)))
                    else:
                        logging.warning('process_vars: cannot find time %s in %s in retry %s of %s' % (ts_esmf,wrfout_path,str(k+1),str(max_retries)))
                        logging.info('process_vars: Available times: %s' % times)
                        logging.info('process_vars: waiting for next retry...')
                        time.sleep(5)
            except Exception as e:
                logging.warning('Exception %s while reading wrfout file %s' % (e, wrfout_path))
                if k == max_retries-1:
                    logging.error('process_vars: cannot open file %s' % wrfout_path)
                    raise PostprocError("process_vars: Unable to open file %s" % wrfout_path)
                else:
                    logging.warning('process_vars: cannot open file %s in retry %s of %s' % (wrfout_path,str(k+1),str(max_retries)))
                    logging.info('process_vars: waiting for next retry...')
                    time.sleep(5)
                    continue

        if tif_proc:
            crs,gt_a,gt_f = ncwrfmeta(d)

        # build an output file per variable
        for var in vars:
            logging.info('process_vars: postprocessing %s for time %s' % (var, ts_esmf))
            try:
                outpath_base = os.path.join(self.output_path, self.product_name + ("-%02d-" % dom_id) + ts_esmf + "-" + var)
                if tif_proc:
                    if is_fire_var(var):
                        geot = gt_f
                    else:
                        geot = gt_a
                    tif_args = {'crs': crs, 'geot': geot, 'tif_path': outpath_base+'.tif'}
                else:
                    tif_args = {}
                kmz_path, raster_path, cb_path, coords, levels, mf_upd = None, None, None, None, None, {}
                if is_windvec(var):
                    kmz_path,raster_path,cb_path,coords,levels = self._vector2kmz(d, var, tndx, ts_esmf, None, outpath_base, cleanup=False, **tif_args)
                else:
                    kmz_path,raster_path,cb_path,coords,levels = self._scalar2kmz(d, var, tndx, ts_esmf, None, outpath_base, cleanup=False, **tif_args)
                if levels is not None:
                    mf_upd['levels'] = levels
                if cb_path is not None:
                    mf_upd['colorbar'] = osp.basename(cb_path)
                mf_upd['kml'] = osp.basename(kmz_path)
                mf_upd['raster'] = osp.basename(raster_path)
                mf_upd['coords'] = coords
                if tif_proc:
                    mf_upd['tif'] = osp.basename(tif_args.get('tif_path'))
                self._update_manifest(dom_id, ts_esmf, var, mf_upd)
            except Exception as e:
                logging.warning("Exception %s while postprocessing %s for time %s" % (e, var, ts_esmf))
                logging.warning(traceback.print_exc())

        if tslist is not None:
            try:
                logging.info('process_vars: postprocessing timeseries for time %s' % ts_esmf)
                ts_paths = tslist.write_timestep(d,dom_id,tndx,ts_esmf)
                self._update_empty_manifest(dom_id, ts_esmf)
            except Exception as e:
                logging.warning("Exception %s while postprocessing timeseries for time %s" % (e, ts_esmf))
                logging.warning(traceback.print_exc())

        d.close()


    def vars2kmz(self, wrfout_path, dom_id, ts_esmf, vars):
        """
        Postprocess a list of scalar fields at a given simulation time into KMZ files.

        :param wrfout_path: WRF file to process
        :param dom_id: the domain identifier
        :param ts_esmf: time stamp in ESMF format
        :param vars: list of variables to process
        """
        # open the netCDF dataset
        d = nc4.Dataset(wrfout_path,'r')

        # extract ESMF string times and identify timestamp of interest
        times = [''.join(x) for x in d.variables['Times'][:].astype(str)]
        if ts_esmf not in times:
            raise PostprocError("vars2kmz: Invalid timestamp %s" % ts_esmf)
        tndx = times.index(ts_esmf)

        # build one KMZ per variable
        for var in vars:
            try:
                outpath_base = os.path.join(self.output_path, self.product_name + ("-%02d-" % dom_id) + ts_esmf + "-" + var)
                kmz_path = None
                if is_windvec(var):
                    kmz_path,_,_,_,_ = self._vector2kmz(d, var, tndx, ts_esmf, outpath_base)
                else:
                    kmz_path,_,_,_,_ = self._scalar2kmz(d, var, tndx, ts_esmf, outpath_base)
                kmz_name = osp.basename(kmz_path)
                self._update_manifest(dom_id, ts_esmf, var, { 'kml' : kmz_name })


            except Exception as e:
                logging.warning("Exception %s while postprocessing %s for time %s into KMZ" % (e, var, ts_esmf))
                logging.warning(traceback.print_exc())


    def vars2png(self, wrfout_path, dom_id, ts_esmf, vars):
        """
        Postprocess a list of scalar fields into KMZ files.

        :param wrfout_path: WRF file to process
        :param dom_id: the domain identifier
        :param ts_esmf: time stamp in ESMF format
        :param vars: list of variables to process
        """
        # open the netCDF dataset
        d = nc4.Dataset(wrfout_path,'r')

        # extract ESMF string times and identify timestamp of interest
        times = [''.join(x) for x in d.variables['Times'][:].astype(str)]
        if ts_esmf not in times:
            raise PostprocError("vars2png: Invalid timestamp %s" % ts_esmf)
        tndx = times.index(ts_esmf)

        # build one KMZ per variable
        for var in vars:
            try:
                outpath_base = os.path.join(self.output_path, self.product_name + ("-%02d-" % dom_id) + ts_esmf + "-" + var)
                if is_windvec(var):
                    raster_path, cb_path, coords, levels = self._vector2png(d, var, tndx, outpath_base)
                else:
                    raster_path, cb_path, coords, levels = self._scalar2png(d, var, tndx, outpath_base)
                mf_upd = { 'raster' : osp.basename(raster_path), 'coords' : coords}
                if cb_path is not None:
                    mf_upd['colorbar'] = osp.basename(cb_path)
                if levels is not None:
                    mf_upd['levels'] = levels
                self._update_manifest(dom_id, ts_esmf, var, mf_upd)
            except Exception as e:
                logging.warning("Exception %s while postprocessing %s for time %s into PNG" % (e, var, ts_esmf))
                logging.warning(traceback.print_exc())

    def process_file(self, wrfout_path, var_list, skip=1):
        """
        Process an entire file, all timestamps and generate images for var_instr.keys().

        :param wrfout_path: the wrfout to process
        :param var_list: list of variables to process
        :param skip: only process every skip-th frame
        """
        traceargs()

        # get domain ID
        dom_id = re.match(r'.*wrfout_d(0[0-9])_[0-9_\-:]{19}', wrfout_path).groups()

        # open the netCDF dataset
        d = nc4.Dataset(wrfout_path,'r')

        # extract ESMF string times and identify timestamp of interest
        times = [''.join(x) for x in d.variables['Times'][:].astype(str)]

        # build one KMZ per variable
        fixed_colorbars = {}
        for tndx, ts_esmf in enumerate(times[::skip]):
            print(('Processing time %s ...' % ts_esmf))
            for var in var_list:
                try:
                    outpath_base = os.path.join(self.output_path, self.product_name + '-' + ts_esmf + '-' + var)
                    if is_windvec(var):
                        kmz_path,raster_path,cb_path,coords,levels = self._vector2kmz(d, var, tndx, ts_esmf, outpath_base, cleanup=False)
                    else:
                        kmz_path,raster_path,cb_path,coords,levels = self._scalar2kmz(d, var, tndx, ts_esmf, outpath_base, cleanup=False)
                    mf_upd = { 'raster' : osp.basename(raster_path), 'coords' : coords, 'kml' : osp.basename(kmz_path) }
                    if cb_path is not None:
                        # optimization for display when we know colorbar has fixed scale
                        # we memoize the colorbar computed at time 0 and use it for all frames
                        # although other colorbars are generated, they are deleted
                        scale = self.wisdom_update.get('scale', get_wisdom(var)['scale'])
                        if type(scale) == list and np.isfinite(scale[0]) and np.isfinite(scale[1]):
                            logging.info("Using fixed colorbar strategy for variable " + var)
                            if tndx == 0:
                                fixed_colorbars[var] = osp.basename(cb_path)
                            else:
                                os.remove(cb_path)
                            mf_upd['colorbar'] = fixed_colorbars[var]
                        else:
                            mf_upd['colorbar'] = osp.basename(cb_path)
                    if levels is not None:
                        mf_upd['levels'] = levels
                    self._update_manifest(dom_id, ts_esmf, var, mf_upd)
                except Exception as e:
                    logging.warning("Exception %s while postprocessing %s for time %s" % (e, var, ts_esmf))
                    logging.warning(traceback.print_exc())


    def _update_empty_manifest(self,dom_id,ts_esmf):
        """
        Adds empty structure to the manifest if not generated yet.

        :param dom_id: the domain id (1, 2, 3, ...)
        :param ts_esmf: ESMF time string
        """
        # update the manifest with the domain/ts_esmf/var info
        dom = self.manifest.get(str(dom_id), {})
        self.manifest[str(dom_id)] = dom

        # extract timestamp
        td = dom.get(ts_esmf, {})
        dom[ts_esmf] = td

        # synchronize the file
        mf_path = os.path.join(self.output_path, self.product_name + '.json')
        json.dump(self.manifest, open(mf_path, 'w'),indent=1, separators=(',',':'))


    def _update_manifest(self,dom_id,ts_esmf,var,kv):
        """
        Adds a key-value set to the dictionary storing metadata for time ts_esmf and variable var.

        :param dom_id: the domain id (1, 2, 3, ...)
        :param ts_esmf: ESMF time string
        :param var: variable name
        :param kv: key-value dictionary to merge
        """
        # update the manifest with the domain/ts_esmf/var info
        dom = self.manifest.get(str(dom_id), {})
        self.manifest[str(dom_id)] = dom

        # extract timestamp
        td = dom.get(ts_esmf, {})
        dom[ts_esmf] = td

        # store variable in timestamp
        vd = td.get(var, {})
        td[var] = vd
        vd.update(kv)

        # synchronize the file
        mf_path = os.path.join(self.output_path, self.product_name + '.json')
        json.dump(self.manifest, open(mf_path, 'w'),indent=1, separators=(',',':'))



if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print(('usage: %s <wrfout_path> <var_instr> <prefix> [skip]' % sys.argv[0]))
        sys.exit(1)

    wrf_path = sys.argv[1]
    var_instr = None
    if sys.argv[2][0] == '@':
        var_instr = json.load(open(sys.argv[2][1:]))
    else:
        var_instr = {x:{} for x in sys.argv[2].split(',')}

    prefix = sys.argv[3]
    skip = 1
    if len(sys.argv) == 5:
        skip = int(sys.argv[4])

    p = Postprocessor(os.path.dirname(prefix), os.path.basename(prefix), var_instr)
    p.process_file(wrf_path, list(var_instr.keys()), skip)


