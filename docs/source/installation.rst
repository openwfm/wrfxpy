Installation
************

Recommended method
==================

Download and install the Python 2 `Anaconda Python <https://www.continuum.io/downloads>`_ distribution for your platform.  We recommend an installation into the users home directory.

Install pre-requisites: 

::

  conda install basemap netcdf4 pyproj paramiko
  conda install --channel https://conda.anaconda.org/IOOS simplekml
  pip install f90nml
  pip install MesoPy

Note that ``conda`` and ``pip`` are package managers available in the Anaconda Python distribution.

Next, clone the *wrfxpy* code:

::
  
  git clone https://github.com/vejmelkam/wrfxpy.git

And finally, configure the system:

::
  
  cd wrfxpy/etc
  cp conf.json.initial conf.json
  <your-favorite-editor-here> conf.json

And tell *wrfxpy* where your WPS and WRFV3 system is located by editing the directories and where it is located.  If you need to move the workspace directory somewhere else, change the ``workspace_dir`` key.


Custom method
=============

A different python distribution and installation can be used if needed.  Below is a list of libraries the system requires:

* `Python 2.7+ <https://www.python.org/download/releases/2.7/>`_
* `Basemap <http://matplotlib.org/basemap/>`_  to render the rasters
* `simplekml <https://simplekml.readthedocs.org/en/latest/>`_ to build KMZ files
* `f90nml <https://pypi.python.org/pypi/f90nml>`_ to manipulate Fortran namelists
* `pyproj <https://pypi.python.org/pypi/pyproj>`_ to place domains dynamically in LCC projection
* `paramiko <https://pypi.python.org/pypi/paramiko>`_ to communicate over SSH with remote hosts
* `netCDF4 <https://pypi.python.org/pypi/netCDF4>`_ to manipulate WPS and WRF files
* `MesoPy <https://pypi.python.org/pypi/MesoPy>`_ to retrieve fuel moisture observations from Mesowest

The simplest way to satisfy these requirements is to install `Anaconda Python <https://www.continuum.io/downloads>`_ and then run the following commands:

::

  conda install basemap netcdf4 pyproj paramiko
  conda install --channel https://conda.anaconda.org/IOOS simplekml
  pip install f90nml

This should install all prerequisites.

*wrfxpy* is installed by cloning a GitHub repository

::

  git clone https://github.com/vejmelkam/wrfxpy.git
