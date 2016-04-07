Forecasting 
******************

The script ``forecast.sh`` serves to run weather forecasts, fire danger forecasts and fire
simulations depending on its settings.

The script requires a JSON configuration file to control its execution, for an example refer
to the :doc:`quickstart`.  The configuration file is JSON dictionary with the keys described
in the following sections.  Not all keys are required.


Domains
=======

Domains are shared by WPS and by WRF and as they are very important, they have their own section.
All domains are given in the key ``domains : [dict]``.

The first domain is always the top-level domain, and subsequent domains are always child domains
such that their parent is always defined before they are.  Each domain can be precomputed or
dynamically placed.  In the following, we detail configuration of each type.

Note: in the input namelist, wrfxpy will set the top level domain as ``specified``, while
all other domains will be ``nested``.

All domains must have the following keys:

* ``history_interval : [integer]`` the history interval in minutes, default is 60
* ``parent_id : [integer]`` the id of the parent domain (can omit for parent domain with id=1)

A top-level precomputed domain must have the following keys:

Example:

::

  "domains" : {
    "1" : {
      "time_step" : 5,
      "precomputed" : "path/to/precomputed/geo_em.dYY.nc",
      "history_interval" : 15
    }


This will cause the use of the domain as precomputed in ``colorado_domain_1.nc`` with a time step of 5 seconds and
a history interval of 15 mins.  All spatial information will be retrieved from the netCDF file.

Top-level dynamically placed
----------------------------

A top-level dynamically placed domain must have the following keys:


* ``time_step : [integer]`` time step in seconds
* ``subgrid_ratio : [integer, integer]`` the refinement ratio for x and y direction for the fire grid, default is 1, 1
* ``cell_size : [integer, integer]`` the size of one grid cell in meters in DX, DY order (placed)
* ``domain_size : [integer, integer]`` the number of grid points in longitudinal and latitudinal direction (placed)
* ``geog_res : [string]`` the resolution of geographical/fuel data to use for the domain (placed)
* ``center_latlon : [float, float]`` the latitude and longitude of grid center (placed)
* ``truelats : [float, float]`` the true lattitudes of the LCC projection
* ``stand_lon : [float]`` the standard longitude of the LCC projection

Example:

::

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
  }
 
This will set up a top level domain domain with cell size 1km and 101x101 grid points, centered on the location [39.894264, -103.903222],
with true latitudes and standard longitude as provided.  The time step will be 5 seconds and history interval will be 15 mins.
High resolution ``.3s`` geographical/fuel data will be used to construct the domain.  The domain will have a fire grid 25m x 25m.

If the domain is precomputed, the following key must be set:

* ``precomputed : [string]`` the precomputed ``geo_em.dYY.nc`` file path relative to wrfxpy installation

Then, all values are loaded from the ``geo_em.dYY.nc`` file except for ``time_step`` and ``history_interval``.

Child domain statically placed
------------------------------

A child domain requires the following keys:

* ``parent_time_step_ratio : [int]`` ratio of child time step to parent time step
* ``parent_cell_size_ratio : [int]`` ratio of the size of the child cell to the parent cell
* ``parent_start : [int, int]`` the x, y coordinates of the parent grid where this grid starts
* ``parent_end : [int, int]`` the x, y coordinates of the parent grid where this grid ends

Example

::

  "domains" : {
    "1" : {
      "time_step" : 5,
      "history_interval" : 30,
      "precomputed" : "precomputed/my_grids/colorado_domain_1.nc"
    },
    "2" : {
      "parent_id" : 1,
      "parent_time_step_ratio" : 4,
      "history_interval" : 30,
      "precomputed" : "precomputed/my_grids/colorado_domain_2.nc"
    }
  }

If the child domain is precomputed, again all these values are read in from the ``geo_em`` file automatically except
timing information: ``parent_time_step_ratio`` must still be set.

* ``precomputed : [string]`` the precomputed ``geo_em.dYY.nc`` file path relative to wrfxpy installation


Child domain placed by bounding box
-----------------------------------

* ``parent_cell_size_ratio : [int]`` the ratio of cell size to parent cell size
* ``parent_time_step_ratio : [int]`` ratio of child time step to parent time step
* ``bounding_box : [float, float, float, float]`` the bounding box the domain should enclose as [min_lon, min_lat, max_lon, max_lat]

Examples

::

  "domains" : {
    "1" : {
      "time_step" : 50,
      "history_interval" : 30,
      "precomputed" : "precomputed/my_grids/colorado_domain_1.nc"
    },
    "2" : {
      "parent_cell_size_ratio" : 3,
      "parent_time_step_ratio" : 3,
      "bounding_box" : [-105, 39, -105.5, 39.5],
      "history_interval" : 15,
      "geog_res" : ".3s",
      "subgrid_ratio" : [50, 50]
      "parent_time_step_ratio: [int]`` 
    }
  }


The value must be a dictionary mapping ``geo_em.dYY.nc`` files to their actual location.


WRF-SFIRE inputs
================

All of the following keys are required.

* ``grid_code : [string]`` the grid code is part of the job id and semantically should identify the configured grid
* ``grib_source : [string]`` must be one of NAM, NARR or HRRR
* ``geogrid_path : [string]`` the path to ``WPS GEOG`` data directory 
* ``start_utc : [esmf_time]`` the start time of the simulation in ESMF format
* ``end_utc : [esmf_time]`` the end time of the simulation in ESMF format

The keys in the remainder of this section are optional.

* ``ignitions : [dict]`` (optional) is a dictionary of domains (string identifier, e.g. "1") to a list of ignitions that should be added to the domain, each being a dictionary with keys as shown in example.  Including this option causes the fire model to be switched on in each domain listed.  A total of five ignitions is allowed (combined for all domains).  
  
.. tip::  
  If a domain is listed without any ignitions, the fire model is switched on and computes quantities
  related to fire danger, such as fire spread rates, fuel moisture values, etc.


Namelist templates
==================

All of the following keys are required.

* ``wps_namelist_path : [string]`` the WPS namelist template
* ``wrf_namelist_path : [string]`` the WRF namelist template
* ``fire_namelist_path : [string]`` the fire namelist template
* ``emissions_namelist_path : [string]`` the file_emissions namelist template


Parallel job configuration
==========================

The following keys are required.

* ``num_nodes : [int]`` the number of parallel nodes to use for WRF execution
* ``ppn : [int]`` the number of processors per node to request
* ``wall_time_hrs : [int]`` the wall time to request from the schedule in hours
* ``qman : [string]`` the queue manager to use, must be ``sge``

Fuel moisture data assimilation
===============================

The key ``fuel_moisture_da`` is optional.  If given, it needs to contain two keys:

* ``domains : [list(int)]`` a list of domains for which to run data assimilation

.. important::
  In addition to this, the file ``etc/tokens.json`` must contain the key ``mesowest_token : [string]``,
  which will be used to access the Mesowest API (you must obtain one here `Mesowest <http://mesowest.org/>`_).

The data assimilation code will download 10-hr fuel moisture observations from stations in the domain area and assimilate them into the equilibrium.

Example::

  "fuel_moisture_da" : {
    "domains" : [ 1 ]
  }


Postprocessing
==============

The key ``postproc``, when present contains a dictionary keyed by domain id (string), which identifies the variables to postprocess for each domain.
For each listed variable, a PNG and a KMZ file is created and if required, a colorbar (configured in ``var_wisdom``).

Example::

  "postproc" : {
    "1" : ["T2", "PSFC", "WINDSPD" ]
  }


  
