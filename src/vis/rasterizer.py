# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


from __future__ import absolute_import
import matplotlib as mpl
import matplotlib.pyplot as plt
import simplekml as kml
from mpl_toolkits.basemap import Basemap, interp
import numpy as np
import netCDF4 as nc4
import sys
import os
import logging
from six.moves import range
try:
    # python 2
    from StringIO import StringIO
except ImportError:
    # python 3
    from io import BytesIO as StringIO

def make_colorbar(rng,orientation,size_in,cmap,cb_label,dpi=200,spacing='proportional',ticks=None,norm=None,ticklabels=None):
    """
    Create a colorbar to accompany a raseter image.

    See: http://matplotlib.org/examples/api/colorbar_only.html for more examples

    :param rng: the minimum/maximum to display on the colorbar
    :param orientation: 'vertical' or 'horizontal'
    :param size_in: larger dimension in inches
    :param cmap: the colormap in use
    :param cb_label: the colorbar label
    :param dpi: optional, dots per inch
    :param spacing: optional, spacing between ticks
    :param ticks: optional, list of ticks to show
    :param norm: optional, scale to plot the data
    :return: a StringIO object with the colorbar as PNG data
    """
    if norm:
        norm = norm(rng[0],rng[1])
    else:
        norm = mpl.colors.Normalize(rng[0],rng[1])

    kwargs = { 'norm': norm,
               'orientation': orientation,
               'spacing': spacing,
               'cmap': cmap,
               'ticks': ticks}

    # build figure according to requested orientation
    hsize, wsize = (size_in,size_in*0.5) if orientation == 'vertical' else (size_in*0.5,size_in)
    fig = plt.figure(figsize=(wsize,hsize))

    # proportions that work with axis title (horizontal not tested)
    ax = fig.add_axes([.5,.03,.12,.8]) if orientation=='vertical' else fig.add_axes([0.03,.4,.8,.12])

    # construct the colorbar and modify properties
    cb = mpl.colorbar.ColorbarBase(ax,**kwargs)
    cb.set_label(cb_label,color='0',fontsize=8,labelpad=-40)
    if ticklabels:
        cb.set_ticklabels(ticklabels)
        levels = ticklabels        
    else:
        levels = cb.get_ticks().tolist()

    # move ticks to left side
    ax.yaxis.set_ticks_position('left')
    for tick_lbl in ax.get_yticklabels():
        tick_lbl.set_color('0')
        tick_lbl.set_fontsize(8)

    # save png to a StringIO
    str_io = StringIO()
    fig.savefig(str_io,dpi=dpi,format='png',transparent=True)
    plt.close()

    return str_io.getvalue(),levels


def make_discrete_colorbar(labels,colors,orientation,size_in,cmap,cb_label,dpi=200):
    """
    Create a discrete colorbar to accompany a scatter raseter image.

    See: http://matplotlib.org/examples/api/colorbar_only.html for more examples

    :param labels: list of labels for each colorbar tick
    :param colors: list of colors for each colorbar tick
    :param orientation: 'vertical' or 'horizontal'
    :param size_in: larger dimension in inches
    :param cmap: the colormap in use
    :param cb_label: the colorbar label
    :param dpi: dots per inch
    :return: a StringIO object with the colorbar as PNG data
    """
    N = len(labels)

    kwargs = { 
        'norm': mpl.colors.Normalize(-.5, N-.5),
        'orientation': orientation,
        'spacing': 'proportional',
        'ticks': list(range(0,N)),
        'cmap': cmap
    }

    # build figure according to requested orientation
    hsize, wsize = (size_in,size_in*0.5) if orientation == 'vertical' else (size_in*0.5,size_in)
    fig = plt.figure(figsize=(wsize,hsize))

    # proportions that work with axis title (horizontal not tested)
    ax = fig.add_axes([.5,.03,.12,.8]) if orientation=='vertical' else fig.add_axes([0.03,.4,.8,.12])

    # construct the colorbar and modify properties
    cb = mpl.colorbar.ColorbarBase(ax,**kwargs)
    cb.ax.tick_params(length=0)
    cb.ax.set_yticklabels(labels)
    cb.set_label(cb_label,color='0',fontsize=8,labelpad=-40)

    # move ticks to left side
    ax.yaxis.set_ticks_position('left')
    for tick_lbl in ax.get_yticklabels():
        tick_lbl.set_color('0')
        tick_lbl.set_fontsize(5)

    # save png to a StringIO
    str_io = StringIO()
    fig.savefig(str_io,dpi=dpi,format='png',transparent=True)
    plt.close()

    return str_io.getvalue()


def basemap_raster_mercator(lon, lat, grid, cmin, cmax, cmap_name, norm=None):
    if norm:
        norm = norm(cmin,cmax)
    
    # longitude/latitude extent
    lons = (np.amin(lon), np.amax(lon))
    lats = (np.amin(lat), np.amax(lat))

    logging.info('basemap_raster_mercator: bounding box %s %s %s %s' % (lons + lats))

    # construct spherical mercator projection for region of interest
    m = Basemap(projection='merc',llcrnrlat=lats[0], urcrnrlat=lats[1],
                llcrnrlon=lons[0],urcrnrlon=lons[1])

    masked_grid = np.ma.array(grid,mask=np.isnan(grid))
    fig = plt.figure(frameon=False,figsize=(12,8),dpi=72)
    plt.axis('off')
    cmap = mpl.cm.get_cmap(cmap_name)
    m.pcolormesh(lon,lat,masked_grid,latlon=True,norm=norm,cmap=cmap,vmin=cmin,vmax=cmax)

    str_io = StringIO()
    plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
    plt.close()

    numpy_bounds = [ (lons[0],lats[0]),(lons[1],lats[0]),(lons[1],lats[1]),(lons[0],lats[1]) ]
    float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
    return str_io.getvalue(), float_bounds


def basemap_barbs_mercator(u,v,lat,lon,grid=None,cmin=0,cmax=0,cmap_name=None,norm=None):

    # lon/lat extents
    lons = (np.amin(lon), np.amax(lon))
    lats = (np.amin(lat), np.amax(lat))

    logging.info('basemap_barbs_mercator: bounding box %s %s %s %s' % (lons + lats))
    
    # construct spherical mercator projection for region of interest
    m = Basemap(projection='merc',llcrnrlat=lats[0], urcrnrlat=lats[1],
                llcrnrlon=lons[0],urcrnrlon=lons[1])

    fig = plt.figure(frameon=False,figsize=(12,8),dpi=72*4)
    plt.axis('off')
    # if cmap_name is set, create speed coloring using all the other parameters
    if cmap_name is not None:
        masked_grid = np.ma.array(grid,mask=np.isnan(grid))
        if norm:
            norm = norm(cmin,cmax)
        cmap = mpl.cm.get_cmap(cmap_name)
        m.quiver(
            lon, lat, u, v, masked_grid, latlon=True, norm=norm,
            cmap=cmap, clim=(cmin,cmax), edgecolor='k', linewidth=.2,
            units='width'
        )
    else:
        m.quiver(lon, lat, u, v, latlon=True, units='width')

    str_io = StringIO()
    plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
    plt.close()

    numpy_bounds = [ (lons[0],lats[0]),(lons[1],lats[0]),(lons[1],lats[1]),(lons[0],lats[1]) ]
    float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
    return str_io.getvalue(), float_bounds


def basemap_scatter_mercator(val, lon, lat, bounds, alphas, cmin, cmax, cmap, size = 2, marker = 's', linewidths = 0, text = False, norm=None):
    if norm:
        norm = norm(cmin,cmax)
    
    # number of scatter elements
    N = len(val)	
    border = .05
    bounds = (bounds[0]-border, bounds[1]+border, bounds[2]-border, bounds[3]+border)
   
    logging.info('basemap_scatter_mercator: bounding box {}'.format(bounds))

    # construct spherical mercator projection for region of interest
    m = Basemap(projection='merc',llcrnrlat=bounds[2], urcrnrlat=bounds[3],
                		llcrnrlon=bounds[0],urcrnrlon=bounds[1])
    
    fig = plt.figure(frameon=False,figsize=(12,8),dpi=72*4)
    plt.axis('off')
    for i in range(N):
        m.scatter(
            lon[i], lat[i], size, c=val[i], latlon=True, marker=marker, norm=norm,
            cmap=cmap, vmin=cmin, vmax=cmax, alpha=alphas[i], linewidths=linewidths, edgecolors='k'
        )
    if text:
        for i in range(N):
            for x1,x2,x3 in zip(lon[i],lat[i],val[i]):
                x,y = m(x1+.05,x2+.05)
                s = '{:.2f}'.format(x3)
                plt.text(x,y,s,fontsize=size/7,fontweight='bold')

    # save png to a StringIO
    str_io = StringIO()
    plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
    plt.close()

    numpy_bounds = [ (bounds[0],bounds[2]),(bounds[1],bounds[2]),(bounds[1],bounds[3]),(bounds[0],bounds[3]) ]
    float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
    return str_io.getvalue(), float_bounds

