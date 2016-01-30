Introduction
************

*wrfxpy* is a set of software modules that provide functionality
related to running `WPS and WRF <http://www.openwfm.org/>`_.

In particular, the modules herein can:

* manipulate wps/input/fire namelists
* download GRIB files from various GRIB sources
* execute geogrid, ungrib, metgrid, real, WRF
* monitor WRF execution
* postprocess `netCDF files <http://www.unidata.ucar.edu/software/netcdf/>` to generate raster images or KML files

Prerequisites
=============

The *wrfxpy* system needs the following software to be installed:

* `Python 2.7+ <https://www.python.org/download/releases/2.7/>`_
* `Basemap <http://matplotlib.org/basemap/>`_
* `simplekml <https://simplekml.readthedocs.org/en/latest/>`_

The simplest way to satisfy these requirements is to install `Anaconda Python <https://www.continuum.io/downloads>`_ and then run the following commands:

  conda install basemap netcdf4 
  conda install --channel https://conda.anaconda.org/IOOS simplekml
  pip install f90nml

This should install all prerequisites.

