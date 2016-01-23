Tutorial
********

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
    "sys_install_dir": "/share_home/mvejmelka/Projects/wrfxpy",
    "num_nodes": 6,
    "ppn": 12,
    "wall_time_hrs": 3,
    "qman": "sge",
    "start_utc": "2016-01-23_08:00:00",
    "end_utc": "2016-01-23_11:00:00",
    "precomputed" : { "geo_em.d01.nc" : "precomputed/colorado/geo_em.d01.nc", "geo_em.d02.nc" : "precomputed/colorado/geo_em.d02.nc" },
    "ignitions" : {
     "2" : [ {
      "start_delay_s" : 600,
      "duration_s" : 240,
      "lat" : 39.894264,
      "long" : -103.903222
       } ]
  }

For a detailed overview of the configuration keys, refer to :doc:`configuration`.


