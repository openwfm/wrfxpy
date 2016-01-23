import matplotlib as mpl
import matplotlib.pyplot as plt
import simplekml as kml
from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc4
import sys
import os
import StringIO


def make_colorbar(rng,orientation,size_in,cmap,cb_label,cb_title,dpi=200):
    """
    rng - [min max]
    orientation - 'vertical' or 'horizontal'
    size_in - larger dimension (height for vertical orientation, width for horizontal)
    cmap - the colormap in use
    units - the colorbar label (on left side vertically)
    dpi - dots per inch

    See: http://matplotlib.org/examples/api/colorbar_only.html for more examples
    """
    kwargs = { 'norm':mpl.colors.Normalize(rng[0],rng[1]),
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
  fig = plt.figure(frameon=False,figsize=(12,8),dpi=72)
  plt.axis('off')
  cmap = mpl.cm.get_cmap(cmap_name)
  m.pcolormesh(lon,lat,masked_grid,latlon=True,cmap=cmap,vmin=cmin,vmax=cmax)

  str_io = StringIO.StringIO()
  plt.savefig(str_io,bbox_inches='tight',format='png',pad_inches=0,transparent=True)
  bounds = [ (lons[0],lats[0]),(lons[1],lats[0]),(lons[1],lats[1]),(lons[0],lats[1]) ]

  return str_io.getvalue(), bounds


