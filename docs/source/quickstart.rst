Quickstart
**********

.. important::

  It is imperative that a working WRF-SFIRE installation is available.
  Please first follow the installation instructions :doc:`installation`.
  

First fire forecast
===================

To perform a fire forecast, the script ``forecast.sh`` has to be executed with
a JSON configuration file as an argument, for example:

::

  ./forecast.sh <json-configuration-file>

An example configuration script is ``examples/simple_fire.json``, also listed here for
convenience.  The script has most of the values filled out but there are some placeholders.

Please set the following values:

  * ``geogrid_path`` should point to the directory with your WPS-GEOG data
  * ``num_nodes`` are the number of nodes to use for the parallel job
  * ``ppn`` the number of processors per node to use
  * ``wall_time_hrs`` number of hours of wall time to reserve for the job
  * ``qsys`` the queueing subsystem id which point into ``etc/clusters.json``

::

  {
    "grid_code": "test",
    "grib_source": "NAM",
    "wps_namelist_path": "etc/nlists/default.wps",
    "wrf_namelist_path": "etc/nlists/default.input",
    "fire_namelist_path": "etc/nlists/default.fire",
    "emissions_namelist_path": "etc/nlists/default.fire_emissions",
    "geogrid_path": "/path/to/your/WPS-GEOG",
    "num_nodes": 10,
    "ppn": 12,
    "wall_time_hrs": 3,
    "qsys": "sge",
    "start_utc": "T-30",
    "end_utc": "T+300",
    "domains" : {
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
    },
    "ignitions" : {
       "1" : [ {
        "start_delay_s" : 600,
        "duration_s" : 240,
        "lat" : 39.894264,
        "long" : -103.903222
       } ]
    },
   "postproc" : {
      "1" : [ "T2", "PSFC", "WINDSPD", "WINDVEC", "FIRE_AREA", "FGRNHFX", "FLINEINT", "SMOKE_INT" ]
    }
  }


This example configuration runs a fire simulation with the following settings:

  - a single domain configuration with a domain placed approximately around an ignition point 
  - use `NAM <http://www.nco.ncep.noaa.gov/pmb/products/nam/>`_ as the source for initial and boundary conditions
  - start at time now minus 30 mins and run a 5 hour simulation
  - use 10 nodes, 12 CPU cores per node, allow a wall time of 3 hrs, the queue manager is ``SGE`` (Sun Grid Engine)
  - use the default WPS/WRF/fire/emissions namelists as base
  - ignite the fire 600s after the start of the simulation and deactivate the ignition after 4 minutes.
  - generate surface temperature maps, wind information and fire fields for domain 1, the where the fire is burning

.. tip::
  To learn how to configure jobs in more detail, refer to :doc:`forecasting`.


