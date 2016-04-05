Quickstart
**********

.. warning::

  It is imperative that a working WRF-SFIRE installation is available.
  Installation and configuration of WRF-SFIRE is outside the scope of this documentation,
  please refer to the `OpenWFM <http://www.openwfm.org>`_ website for instructions.


System setup
============
  First, the system must be configured.  Copy the file ``etc/conf.json.initial`` into
  ``etc/conf.json`` and set the paths therein to point to the three correct directories:
  ``WPS`` directory, ``WRF`` directory and the directory ``wrfxpy`` is installed in.


First fire forecast
===================

To perform a fire forecast, the script ``src/forecast.py`` has to be executed with
a JSON configuration file as an argument, for example:

::

  python src/forecast.py examples/colorado.json

An example configuration script is ``examples/colorado.json`` in the working directory,
also listed here:

**NOTE** The ``geogrid_path`` below must be replaced with the path to your WPS-GEOG data

::

  {
    "grid_code": "colo2dv1",
    "grib_source": "NAM",
    "wps_namelist_path": "etc/nlists/colorado-3k.wps",
    "wrf_namelist_path": "etc/nlists/colorado-3k.input",
    "fire_namelist_path": "etc/nlists/colorado-3k.fire",
    "emissions_namelist_path": "etc/nlists/colorado-3k.fire_emissions",
    "geogrid_path": "/path/to/WPS geo data",
    "num_nodes": 10,
    "ppn": 12,
    "wall_time_hrs": 3,
    "qman": "sge",
    "start_utc": "T-30",
    "end_utc": "T+300",
    "domains" : {
      "1" : {
        "time_step" : 60,
        "cell_size" : [12000,12000],
        "domain_size" : [101, 101],
        "center_latlon" : [39, -105.663],
        "truelats" : [37, 40.5],
        "stand_lon" : -105.5,
        "history_interval" : 120,
        "subgrid_ratio" :[1,1],
        "geog_res" : "30s"
      },
      "2" : {
        "parent_id" : 1,
        "parent_cell_size_ratio" : 4,
        "parent_time_step_ratio" : 4,
        "subgrid_ratio" : [1, 1],
        "parent_start" : [30, 30],
        "parent_end" : [71, 71],
        "history_interval" : 120,
        "geog_res" : "30s"
      },
      "3" : {
        "parent_id" : 2,
        "parent_cell_size_ratio" : 3,
        "parent_time_step_ratio" : 3,
        "subgrid_ratio" : [40, 40],
        "bounding_box" : [-104.5, 39.1, -103.8, 40.5],
        "history_interval" : 30,
        "geog_res" : ".3s"
      }
    },
    "ignitions" : {
       "3" : [ {
        "start_delay_s" : 600,
        "duration_s" : 240,
        "lat" : 39.894264,
        "long" : -103.903222
         } ]
    },
    "postproc" : {
        "3" : [ "T2", "WINDSPD", "WINDVEC", "FIRE_AREA", "SMOKE_INT", "FGRNHFX", "FLINEINT" ]
    }
  }

This example configuration runs a fire simulation with the following settings:

  - a three domain configuration, where the first two are statically placed and the third is computed
    from a required bounding box as a child of the second
  - use `NAM <http://www.nco.ncep.noaa.gov/pmb/products/nam/>`_ as the source for initial and boundary conditions
  - start the simulation at time now minus 30 mins and execute a 5 hour simulation
  - use 10 nodes, 12 CPU cores per node, allow a wall time of 3 hrs, the queue manager is ``SGE`` (Sun Grid Engine)
  - use the WPS/WRF/fire namelist templates as indicated (a template is simply a valid namelist)
  - ignite the fire 600s after the start of the simulation  (8:10 AM UTC) at the given location and deactivate the ignition after 4 minutes.
  - generate surface temperature maps, wind information and fire fields in domain 3, the where the fire is burning
   
For a detailed overview of the configuration keys, refer to :doc:`configuration`.


