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


Postprocessing
==============

The script ``postprocess.sh`` accepts four arguments, the wrfout file to process,
the variables to postprocess (or an instruction file, see below), the prefix on which
to base the filenames and the skip (the script will process every skip-th frame).
The script always generates PNG files and KMZ files for each variable and timestamp.

Example::

  ./postprocess.sh /path/to/wrfout T2,PSFC my_directory/file_prefix 1

Alternatively, instead of listing the variables, a more detailed configuration controlling
the colormaps, ranges and other parameters can be specified::

  ./postprocess.sh /path/to/wrfout @var_instructions my_directory/file_prefix 1

Where the file ``var_instructions`` contains::

  {
    "FGRNHFX" : {
        "name" : "Grnd Heat flux",
        "colorbar" : "W/m^2",
        "colormap" : "jet",
        "transparent_values" : [0, 1],
        "scale" : [0, 6]
    }
  }

Will show the colorbar in ``W/m^2`` units and change the displayed variable name to
``Grnd Heat flux``, set the colormap to ``jet``, ensure that values between 0 and 1
are not shown and fix the scale from 0 to 6.

.. tip::
  For the default and more information on values that can be set, examine ``src/vis/var_wisdom.py``.


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


SSH Shuttle
===========

The script ``ssh_shuttle.sh`` accepts a local directory a remote directory name and an identifier
and uploads the entire local directory with simulation results to the remote host configured in ``conf.json`` and registers the simulation in the ``catalog.json`` file on the remote server.

Examples::

  ./ssh_shuttle.sh wksp/mu-simulation/products test_fire_april test_fire_april

The script scans all the files in ``wksp/mu-simulation/products`` and uses SFTP to put them onto the remote hots.  The remote directory must be either an absolute path or (recommended) should be relative to the remote host root setup in ``conf.json``.  The identifier will be used as the description and also as the key under which the simulation is stored in ``catalog.json`` on the remote host.

