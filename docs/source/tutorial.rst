Tutorial
********
.. warning::

  It is imperative that a working WRF-SFIRE installation is available.
  Installation and configuration of WRF-SFIRE is outside the scope of this documentation,
  please refer to the `OpenWFM <http://www.openwfm.org>`_ website for instructions.


To perform a fire forecast, the script ``forecast.py`` has to be executed with
a JSON configuration file as an argument, for example:

::

  ./forecast.py examples/colorado.json

An example configuration script is ``examples/colorado.json`` in the working directory,
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
    "sys_install_dir": "/share_home/mvejmelka/Projects/wrfxpy-dev",
    "num_nodes": 10,
    "ppn": 12,
    "wall_time_hrs": 3,
    "qman": "sge",
    "start_utc": "2016-01-23_08:00:00",
    "end_utc": "2016-01-23_12:00:00",
    "precomputed" : {
      "geo_em.d01.nc" : "precomputed/colorado/geo_em.d01.nc",
      "geo_em.d02.nc" : "precomputed/colorado/geo_em.d02.nc"
    },
    "ignitions" : {
      "1" : [],
      "2" : [ {
        "start_delay_s" : 600,
        "duration_s" : 240,
        "lat" : 39.894264,
        "long" : -103.903222
        } ]
    },
    "postproc" : {
        "1" : [ "T2" ],
        "2" : [ "T2" ]
    }
  }


This example configuration runs a fire simulation with the following settings:

  - WPS and WRF in the directory ``/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/``
  - use `HRRR <http://ruc.noaa.gov/hrrr/>`_ as the source for initial and boundary conditions
  - start the simulation at 8 AM UTC 23rd Jan 2016 and simulate for 4 hours
  - do not run ``geogrid.exe`` but rather use the precomputed domains (obtained by running ``geogrid.exe`` previously)
  - use 10 nodes, 12 CPU cores per node, allow a wall time of 3 hrs, the queue manager is ``SGE`` (Sun Grid Engine)
  - use the WPS/WRF/fire namelist templates as indicated (a template is simply a valid namelist that will be manipulated by the system)
  - ignite the fire in domain 2 600s after the start of the simulation  (8:10 AM UTC) at the given location and deactivate the ignition after 4 minutes.
  - in domain 1, don't ignite any fires but switch on the fire model (computes spread rates, etc.)
  - generate surface temperature maps (T2) for domain 1 and for domain 2
   
For a detailed overview of the configuration keys, refer to :doc:`configuration`.


