Standalone scripts
******************

Although *wrfxpy* is meant to be an integrated system, some functionalities
are exposed through separate scripts.  These are detailed in this section.

Domain setup
============

The script ``domain_setup.sh`` accepts a domain configuration description and
injects the domain configuration into a WPS nmelist file and into an input
namelist file.  Please refer to the domain configuration description in :doc:`forecasting`.

Example::
  ./domain_setup.sh my_domains.json namelist.wps namelist.input

Assuming that ``my_domains.json`` contains the following::

  {
    "1" : {
      "cell_size" : [1000, 1000],
      "domain_size" : [91, 91],
      "center_latlon" : [39.1, -105.9],
      "truelats" : [38.5, 39.6],
      "stand_lon" : -105.9,
      "time_step" : 5,
      "history_interval" : 15,
      "geog_res" : "0.3s",
      "subgrid_ratio" : [50, 50]
    }
  } 

Then both namelists will be setup for a single-domain configuration (1km grid
cell sie, 91 x 91 domain size, 20m fire grid).

.. note::
  The namelist files are *overwritten*.

Grib retrieval
==============

The script ``grib_retr.sh`` accepts fourth arguments, the grib source identifier,
the UTC start, the end time of a simulation in ESMF format and the ingest directory.

Example::
  ./grib_retr.sh HRRR 2016-03-26_14:00:00 2016-03-26_19:00:00 ingest

This will find out which GRIB2 files are required to perform this simulation and
will download them into subdirectories of the ``ingest`` directory.

.. tip::
  Using the *wrfxpy* ingest directory (or the same directory) consistently will make
  best use of the transparent local caching functionality.  Any files that have already
  been downloaded are not re-downloaded.


Fuel moisture DA
================

The script ``apply_fmda.sh`` accepts a single wrfinput path argument and
performs a data assimilation step using background covariance.

Example::
  ./apply_fmda.sh wrfinput_d01

The script will read in the timestamp from the wrfinput file, determine it's
physical extent (lat/lon) and download all observations of 10-hr fuel moisture
valid at that time available in the region.  Then the equilibrium fuel moisture
content is computed and adjusted with respect to the observations using the
background covariance.  The updated values are written back into the fuel moisture
file.


