.. wrfxpy documentation master file, created by
   sphinx-quickstart on Wed Jan 13 17:40:44 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

wrfxpy
======

*wrfxpy* is a set of software modules that provide functionality
related to running `WPS and WRF <http://www.openwfm.org/>`_.

In particular, the modules herein can:

* manipulate wps/input/fire namelists
* place and setup domains dynamically 
* download GRIB files from various GRIB sources
* execute geogrid, ungrib, metgrid, real, WRF
* monitor WRF execution
* perform fuel moisture data assimilation using RAWS observations from the Mesowest network
* postprocess `netCDF files <http://www.unidata.ucar.edu/software/netcdf/>`_ to generate raster images or KMZ raster files
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

