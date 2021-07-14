Catalog and manifest
********************

This document describes the catalog format and the manifest format.
The catalog file collects computed simulations for the visualizaion server
and points to the manifest file for each simulation.
The manifest contains postprocessed rasters pertaining to a single simulation.

Catalog
=======

The file ``catalog.json`` in the root visualization directory of the *wrfxweb*
visualization system is a JSON file which stores the following information about each simulation:

* ``manifest_path : [string]`` the path to the manifest string
* ``description : [string]`` a description that is shown in the selection menu to a user
* ``from_utc : [esmf_time]`` the start time of the simulation in ESMF format
* ``to_utc : [esmf_time]`` the end time of the simulation in ESMF format

Example::

  { 
    "patch_fire": {
      "manifest_path": "patch3/patch.json",
      "description": "Patch Springs Fire [UT]",
      "to_utc": "2013-08-19_09:00:00",
      "from_utc": "2013-08-11_00:00:00"
    },
    "test_fire_3": {
      "manifest_path": "test_fire_3/wfc-two-domain-fire.json",
      "description": "2-domain test fire, viscosity=0",
      "to_utc": "2016-04-08_23:00:00",
      "from_utc": "2016-04-08_18:00:00"
    },
    .
    .
    .
  } 


Manifest
========

The manifest file is a JSON file that collects information on which domains, timestamps
and variables are generated from a simulation.  The top-level object is a dictionary keyed
by string domain identifier ("1", "2", ...) and contains an object keyed by ESMF time
("2016-03-30_00:00:00", ...) which in turn contains a dictionary keyed by variable names
(e.g. "T2", "WINDVEC", ...).  The postprocessing results for each variable are represented
by dictionary keys as follows:

* ``colobar : [string]`` (optional) path to the colorbar, if any
* ``raster : [string]`` path to the display raster (PNG file)
* ``kml : [string]`` (optional) path to the KMZ file for possible download
* ``coords : [string]`` the corner coordinates of the PNG file (geolocation)

An partial example of one variable in a file is below:

::

  {
    "1" : {
      .
      .
      "2016-03-30_16:15:00" : {
      .
      .
        "FMC_G": {
          "colorbar": "patch-2013-08-14_12:00:00-FMC_G-cb.png",
          "raster": "patch-2013-08-14_12:30:00-FMC_G-raster.png",
          "kml": "patch-2013-08-14_12:30:00-FMC_G.kmz",
          "coords": [
              [ -113.15496826171875, 39.978614807128906 ],
              [ -112.26193237304688, 39.978614807128906 ],
              [ -112.26193237304688, 40.65946960449219 ],
              [ -113.15496826171875, 40.65946960449219 ]
          ]
        }
      }
    }
  }


