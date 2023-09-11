.. wrfx documentation master file, created by
   sphinx-quickstart on Fri Jul 30 16:57:06 2021.

wrfxpy
======

*wrfxpy* is a set of software modules that provide functionality
related to running `WPS and WRF <http://www.openwfm.org/>`_.

In particular, the modules herein can:

* manipulate wps/input/fire namelists
* place and setup domains dynamically 
* download GRIB files from various GRIB sources
* download L2 Active Fires satellite data from different sources
* execute geogrid, ungrib, metgrid, real, WRF
* monitor WRF execution
* perform fuel moisture data assimilation using RAWS observations from the Mesowest network
* postprocess `netCDF files <http://www.unidata.ucar.edu/software/netcdf/>`_ to generate raster images, KMZ raster files, and GeoTIFF files
* postprocess satellite data files to generate raster images, KMZ raster files, and GeoTIFF files
* assemble simulation outputs into coherent packages for visualization and synchronize them with a remote *wrfxweb* server


Basic topics
------------

.. toctree::
   :maxdepth: 1

   installation
   quickstart
   forecasting
   standalones

Advanced topics
---------------

.. toctree::
   :maxdepth: 1

   catalog


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

