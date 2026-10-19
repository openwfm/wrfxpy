wrfxpy
******

|docs|

*wrfxpy* is a set of software modules that provide functionality
related to running `WPS and WRF <http://www.openwfm.org/>`_.

In particular, the code enables the user to:

* manipulate wps/input/fire namelists
* place and setup domains dynamically 
* download GRIB files from various GRIB sources
* execute geogrid, ungrib, metgrid, real, WRF
* monitor WRF execution
* postprocess `netCDF files <http://www.unidata.ucar.edu/software/netcdf/>`_ to generate raster images or KMZ raster files
* arrange simulation results into catalogs for display `wrfxweb <https://github.com/vejmelkam/wrfxweb>`_

FMDA with pretrained RNN models
==============================

Follow the `RNN installation guide <docs/rnn_integration.md>`_ for package
installation, pretrained model assets, and test commands. The
`FMDA setup guide <https://wiki.openwfm.org/wiki/Running_FMDA_in_the_WRFx_system>`_
describes observation access and operational configuration.

Copy the supplied JSON configuration files over the blank configuration files
listed below, then edit the destination copies for your installation. Paths
in this table are relative to the wrfxpy repository root. Keep the supplied
originals unchanged; when updating an existing installation, retain its
configured destination files.

.. list-table:: Configuration files to copy
   :header-rows: 1

   * - Supplied file
     - Destination
   * - ``etc/conf.json.initial``
     - ``etc/conf.json``
   * - ``etc/rnn_cycler.json.initial``
     - ``etc/rnn_cycler.json``
   * - ``etc/tokens.json.initial``
     - ``etc/tokens.json``
   * - Either ``etc/fmda_cycler_conus.json`` or ``etc/fmda_cycler_all.json``
     - ``etc/fmda_cycler.json``

Choose one FMDA preset. The CONUS preset supplies one region and a lookback
length of 0; the all-regions preset supplies CONUS plus regional domains and
a lookback length of 4. Use the supplied settings as the starting point,
including regions, lookback length, workspace paths, downloads, and retries.

Edit the destination copies as follows:

* In ``etc/conf.json``, replace the installation-path placeholders needed by
  FMDA, including ``sys_install_path`` and ``wget``. Configure the ``shuttle_*``
  settings for remote publication, or remove those settings for local output.
* In ``etc/rnn_cycler.json``, set ``rnn_model_path`` and ``transfer_path`` to
  the directories containing the supplied pretrained assets.
* In ``etc/tokens.json``, set ``mesowest`` to the observation-access token.
  The ``appkey`` entry is for satellite access.
* In ``etc/fmda_cycler.json``, retain the chosen preset's defaults unless the
  intended run requires different values.

From the wrfxpy repository root, validate the edited JSON files before running
tests or the cycler::

   python -m json.tool etc/conf.json >/dev/null
   python -m json.tool etc/rnn_cycler.json >/dev/null
   python -m json.tool etc/tokens.json >/dev/null
   python -m json.tool etc/fmda_cycler.json >/dev/null

These commands check JSON syntax. The configured paths and credentials must
also be valid for the installation.

.. |docs| image:: https://readthedocs.org/projects/wrfx/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: https://wrfx.readthedocs.org/en/latest/



