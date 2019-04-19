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


import matplotlib as mpl
import matplotlib.pyplot as plt
import simplekml as kml
from mpl_toolkits.basemap import Basemap, interp
import numpy as np
import netCDF4 as nc4
import sys
import os
import StringIO
import logging


def make_colorbar(rng,orientation,size_in,cmap,cb_label,dpi=200):
    """
    Create a colorbar to accompany a raseter image.

    See: http://matplotlib.org/examples/api/colorbar_only.html for more examples

    :param rng: the minimum/maximum to display on the colorbar
    :param orientation: 'vertical' or 'horizontal'
    :param size_in: larger dimension in inches
    :param cmap: the colormap in use
    :param cb_label: the colorbar label
    :param dpi: dots per inch
    :return: a StringIO object with the colorbar as PNG data
    """
    kwargs = { 'norm': mpl.colors.Normalize(rng[0],rng[1]),
               'orientation':orientation,
               'spacing':'proportional',
               'cmap':cmap}

    # build figure according to requested orientation
    hsize, wsize = (size_in,size_in*0.5) if orientation == 'vertical' else (size_in*0.5,size_in)
    fig = plt.figure(figsize=(wsize,hsize))

    # proportions that work with axis title (horizontal not tested)
    ax = fig.add_axes([.5,.03,.12,.8]) if orientation=='vertical' else fig.add_axes([0.03,.4,.8,.12])

    # construct the colorbar and modify properties
    cb = mpl.colorbar.ColorbarBase(ax,**kwargs)
    cb.set_label(cb_label,color='1',fontsize=8,labelpad=-40)

    # move ticks to left side
    ax.yaxis.set_ticks_position('left')
    for tick_lbl in ax.get_yticklabels():
        tick_lbl.set_color('1')
        tick_lbl.set_fontsize(8)

    # save png to a StringIO
    str_io = StringIO.StringIO()
    fig.savefig(str_io,dpi=dpi,format='png',transparent=True)
    plt.close()

    return str_io.getvalue()


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

    kwargs = { 'norm': mpl.colors.Normalize(-.5,N-.5),
	       'orientation': orientation,
               'spacing': 'proportional',
	       'ticks': range(0,N),
               'cmap': cmap}

    # build figure according to requested orientation
    hsize, wsize = (size_in,size_in*0.5) if orientation == 'vertical' else (size_in*0.5,size_in)
    fig = plt.figure(figsize=(wsize,hsize))

    # proportions that work with axis title (horizontal not tested)
    ax = fig.add_axes([.5,.03,.12,.8]) if orientation=='vertical' else fig.add_axes([0.03,.4,.8,.12])

    # construct the colorbar and modify properties
    cb = mpl.colorbar.ColorbarBase(ax,**kwargs)
    cb.ax.tick_params(length=0)
    cb.ax.set_yticklabels(labels)
    cb.set_label(cb_label,color='1',fontsize=8,labelpad=-40)

    # move ticks to left side
    ax.yaxis.set_ticks_position('left')
    for tick_lbl in ax.get_yticklabels():
        tick_lbl.set_color('1')
        tick_lbl.set_fontsize(5)

    # save png to a StringIO
    str_io = StringIO.StringIO()
    fig.savefig(str_io,dpi=dpi,format='png',transparent=True)
    plt.close()

    return str_io.getvalue()


def basemap_raster_mercator(lon, lat, grid, cmin, cmax, cmap_name):

    # longitude/latitude extent
    lons = (np.amin(lon), np.amax(lon))
    lats = (np.amin(lat), np.amax(lat))

    # construct spherical mercator projection for region of interest
    m = Basemap(projection='merc',llcrnrlat=lats[0], urcrnrlat=lats[1],
                llcrnrlon=lons[0],urcrnrlon=lons[1])

    #vmin,vmax = np.nanmin(grid),np.nanmax(grid)
    masked_grid = np.ma.array(grid,mask=np.isnan(grid))
    # logging.info('basemap_raster_mercator: not masked %s %s' % (grid.count(),masked_grid.count()))
    fig = plt.figure(frameon=False,figsize=(12,8),dpi=72)
    plt.axis('off')
    cmap = mpl.cm.get_cmap(cmap_name)
    m.pcolormesh(lon,lat,masked_grid,latlon=True,cmap=cmap,vmin=cmin,vmax=cmax)

    str_io = StringIO.StringIO()
    plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
    plt.close()

    numpy_bounds = [ (lons[0],lats[0]),(lons[1],lats[0]),(lons[1],lats[1]),(lons[0],lats[1]) ]
    float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
    return str_io.getvalue(), float_bounds


def basemap_barbs_mercator(u,v,lat,lon):

    # lon/lat extents
    lons = (np.amin(lon), np.amax(lon))
    lats = (np.amin(lat), np.amax(lat))

    # construct spherical mercator projection for region of interest
    m = Basemap(projection='merc',llcrnrlat=lats[0], urcrnrlat=lats[1],
                llcrnrlon=lons[0],urcrnrlon=lons[1])

    #vmin,vmax = np.nanmin(grid),np.nanmax(grid)
    fig = plt.figure(frameon=False,figsize=(12,8),dpi=72*4)
    plt.axis('off')
    m.quiver(lon,lat,u,v,latlon=True)

    str_io = StringIO.StringIO()
    plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
    plt.close()

    numpy_bounds = [ (lons[0],lats[0]),(lons[1],lats[0]),(lons[1],lats[1]),(lons[0],lats[1]) ]
    float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
    return str_io.getvalue(), float_bounds


def basemap_scatter_mercator(val, lon, lat, bounds, alphas, cmin, cmax, cmap):
    # number of scatter elements
    N = len(val)	
   
    # construct spherical mercator projection for region of interest
    m = Basemap(projection='merc',llcrnrlat=bounds[2], urcrnrlat=bounds[3],
                		llcrnrlon=bounds[0],urcrnrlon=bounds[1])
    
    fig = plt.figure(frameon=False,figsize=(12,8),dpi=72*2)
    plt.axis('off')
    for i in range(N):
    	m.scatter(lon[i],lat[i],50,c=val[i],latlon=True,marker='.',cmap=cmap,vmin=cmin,vmax=cmax,alpha=alphas[i],linewidths=0)

    # save png to a StringIO
    str_io = StringIO.StringIO()
    plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
    plt.close()

    numpy_bounds = [ (bounds[0],bounds[2]),(bounds[1],bounds[2]),(bounds[1],bounds[3]),(bounds[0],bounds[3]) ]
    float_bounds = [ (float(x), float(y)) for x,y in numpy_bounds ]
    return str_io.getvalue(), float_bounds

