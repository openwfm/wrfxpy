Introduction
************

*wrfxpy* is a set of software modules that provide functionality
related to running `WPS and WRF <http://www.openwfm.org/>`_.

In particular, the modules herein can:

* manipulate wps/input/fire namelists
* place and setup domains dynamically 
* download GRIB files from various GRIB sources
* execute geogrid, ungrib, metgrid, real, WRF
* monitor WRF execution
* postprocess `netCDF files <http://www.unidata.ucar.edu/software/netcdf/>`_ to generate raster images or KMZ raster files
* update catalogs with products

Prerequisites
=============

The *wrfxpy* system needs the following software to be installed:

* `Python 2.7+ <https://www.python.org/download/releases/2.7/>`_
* `Basemap <http://matplotlib.org/basemap/>`_  to render the rasters
* `simplekml <https://simplekml.readthedocs.org/en/latest/>`_ to build KMZ files
* `f90nml <https://pypi.python.org/pypi/f90nml>`_ to manipulate Fortran namelists
* `pyproj <https://pypi.python.org/pypi/pyproj>`_ to place domains dynamically in LCC projection

The simplest way to satisfy these requirements is to install `Anaconda Python <https://www.continuum.io/downloads>`_ and then run the following commands:

::

  conda install basemap netcdf4 pyproj paramiko
  conda install --channel https://conda.anaconda.org/IOOS simplekml
  pip install f90nml

This should install all prerequisites.

*wrfxpy* is installed by cloning a GitHub repository

::

  git clone https://github.com/vejmelkam/wrfxpy.git
