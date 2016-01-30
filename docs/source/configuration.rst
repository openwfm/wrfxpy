Configuration keys
******************

The script `forecast.py` requires a JSON configuration file to control its execution,
for an example refer to the :doc:`tutorial`.  The configuration file is JSON dictionary
with the keys described in the following sections.  Not all keys are required.

Directories
===========

All of the following keys are required.

* ``workspace_dir : [string]`` the path in which jobs are executed
* ``wps_install_dir : [string]`` the path to a working installation of WPS
* ``wrf_install_dir : [string]`` the path to a working installation of WRF
* ``sys_install_dir : [string]`` the system installation directory

WRF-SFIRE inputs
================

All of the following keys except ``precomputed`` are required.

* ``grid_code : [string]`` the grid code is part of the job id and semantically should identify the configured grid
* ``grib_source : [string]`` must be HRRR
* ``geogrid_path : [string]`` the path to ``WPS GEOG`` data
* ``start_utc : [esmf_time]`` the start time of the simulation in ESMF format
* ``end_utc : [esmf_time]`` the end time of the simulation in ESMF format

The keys in the remainder of this section are optional.  Without the ignitions key, the fire section of
the template ``namelist.input`` will remain unaffected.  If the ``ignitions``
key is included, the fire model is switched on.

* ``ignitions : [dict]`` (optional) is a dictionary of domains (string identifier, e.g. "1") to a list of ignitions that should be added to the domain, each being a dictionary with keys as shown in example.  Including this option causes the fire model to be switched on in each domain with ignitions.  Limitation: each ignition can only be setup in one domain.
* ``precomputed : [dict]`` (optional) the system has support for precomputed geogrids, which can be linked in instead of running ``geogrid.exe``.  The value must be a dictionary mapping ``geo_em.dYY.nc`` files to their actual location.


Namelist templates
==================

All of the following keys are required.

* ``wps_namelist_path : [string]`` the WPS namelist template
* ``wrf_namelist_path : [string]`` the WRF namelist template
* ``fire_namelist_path : [string]`` the fire namelist template


Parallel job configuration
==========================

The following keys are compulstory.

* ``num_nodes : [int]`` the number of parallel nodes to use for WRF execution
* ``ppn : [int]`` the number of processors per node to request
* ``wall_time_hrs : [int]`` the wall time to request from the schedule in hours
* ``qman : [string]`` the queue manager to use, must be ``sge``


Postprocessing
==============

The key ``postproc``, when present contains a dictionary keyed by domain id (string),
which identifies the variables to postprocess for each domain.

  
