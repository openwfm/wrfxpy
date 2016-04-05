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
    "wps_install_dir": "/replace/with/abs/path/to/WPS",
    "wrf_install_dir": "/replace/with/abs/path/to/WRFV3",
    "grib_source": "HRRR",
    "wps_namelist_path": "etc/nlists/colorado-3k.wps",
    "wrf_namelist_path": "etc/nlists/colorado-3k.input",
    "fire_namelist_path": "etc/nlists/colorado-3k.fire",
    "emissions_namelist_path": "etc/nlists/colorado-3k.fire_emissions",
    "geogrid_path": "/replace/with/path/to/WPS-GEOG",
    "sys_install_dir": "/replace/with/path/to/installation/directory/of/wrfxpy",
    "num_nodes": 10,
    "ppn": 12,
    "wall_time_hrs": 3,
    "qman": "sge",
    "start_utc": "T-30",
    "end_utc": "T+180",
    "domains" : {
      "1" : {
        "cell_size" : [1000, 1000],
        "domain_size" : [101, 101],
        "center_latlon" : [39.894264, -103.903222],
        "truelats" : [39.2, 40.5],
        "stand_lon" : -103.903222,
        "time_step" : 5,
        "history_interval" : 15,
        "geog_res" : ".3s",
        "subgrid_ratio" : [40, 40]
      }
    },
    "ignitions" : {
      "1" : [ {
        "start_delay_s" : 600,
        "duration_s" : 240,
        "lat" : 39.894264,
        "long" : -103.903222
        } ]
    },
    "email_notifications" : {
      "to" : "your@email.here",
      "events" : ["start", "wrf_submit", "complete" ]
    },
    "postproc" : {
        "1" : [ "T2", "WINDSPD", "WINDVEC", "FIRE_AREA", "SMOKE_INT", "FLINEINT", "FGRNHFX" ]
    }
  }


This example configuration runs a fire simulation with the following settings:

  - WPS and WRF in the given directories
  - use `HRRR <http://ruc.noaa.gov/hrrr/>`_ as the source for initial and boundary conditions
  - start the simulation at time now minus 30 mins and run for 180 mins
  - place a single domain centered on the (single) ignition point
  - use 10 nodes, 12 CPU cores per node, allow a wall time of 3 hrs, the queue manager is ``SGE`` (Sun Grid Engine)
  - use the WPS/WRF/fire namelist templates as indicated (a template is simply a valid namelist)
  - ignite the fire 600s after the start of the simulation  (8:10 AM UTC) at the given location and deactivate the ignition after 4 minutes.
  - generate surface temperature maps, wind information and fire fields in the domain
   
For a detailed overview of the configuration keys, refer to :doc:`configuration`.


