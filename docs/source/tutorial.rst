Tutorial
********

To perform a fire forecast, the script ``forecast.py`` has to be executed with
a JSON configuration file as an argument, for example:

::
  ./forecast.py colorado.json

An example configuration script is ``colorado.json`` in the working directory,
also listed here:

::
  {
    "grid_code": "colo2dv1",
    "workspace_dir": "wksp",
    "wps_install_dir": "/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/WPS",
    "wrf_install_dir": "/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/WRFV3",
    "grib_source": "HRRR",
    "wps_namelist_path": "etc/nlists/colorado-3k.wps",
    "wrf_namelist_path": "etc/nlists/colorado-3k.input",
    "fire_namelist_path": "etc/nlists/colorado-3k.fire",
    "geogrid_path": "/share_home/mvejmelka/Packages/WPS-GEOG",
    "sys_install_dir": "/share_home/mvejmelka/Projects/wrfxpy",
    "num_nodes": 6,
    "ppn": 12,
    "wall_time_hrs": 3,
    "qman": "sge",
    "start_utc": "2016-01-23_08:00:00",
    "end_utc": "2016-01-23_11:00:00",
    "precomputed" : { "geo_em.d01.nc" : "precomputed/colorado/geo_em.d01.nc", "geo_em.d02.nc" : "precomputed/colorado/geo_em.d02.nc" }
  }


Configuration keys
==================

* ``grid_code : [string]`` the grid code is part of the job id and semantically should identify the configured grid
* ``workspace_dir : [string]`` the path in which jobs are executed
* ``wps_install_dir : [string]`` the path to a working installation of WPS
* ``wrf_install_dir : [string]`` the path to a working installation of WRF
* ``grib_source : [string]`` must be HRRR
* ``wps_namelist_path : [string]`` the WPS namelist template
* ``wrf_namelist_path : [string]`` the WRF namelist template
* ``fire_namelist_path : [string]`` the fire namelist template
* ``geogrid_path : [string]`` the path to ``WPS GEOG`` data
* ``sys_install_dir : [string]`` the system installation directory
* ``num_nodes : [int]`` the number of parallel nodes to use for WRF execution
* ``ppn : [int]`` the number of processors per node to request
* ``wall_time_hrs : [int]`` the wall time to request from the schedule in hours
* ``qman : [string]`` the queue manager to use, must be ``sge``
* ``start_utc : [esmf_time]`` the start time of the simulation in ESMF format
* ``end_utc : [esmf_time]`` the end time of the simulation in ESMF format
* ``precomputed : [dict]`` (optional) the system has support for precomputed geogrids



  
