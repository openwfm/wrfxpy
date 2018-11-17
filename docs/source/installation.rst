Installation
************

Recommended method
==================
We recommend using the `Anaconda Python <https://www.continuum.io/downloads>`_ distribution.
Alternative installation instructions are given at the bottom.

WPS/WRF-SFIRE
-------------
Please follow the instructions on `OpenWFM <http://www.openwfm.org>`_ to run WPS/WRF with real data.
Ensure that you have working WPS/WRF installation is available and compatible with the MPI on your cluster.
Note that *wrfxpy* will not modify the WPS/WRF installation, instead for each job, it will clone their directories
into its workspace.

.. attention::

  Past this point, you should be able to run a fire simulation yourself,
  that is have a working ``WPS/WRF-SFIRE`` installation including ``WPS-GEOG``
  and fuels/topography downloaded.  You should be able to submit a parallel
  job into your cluster/supercomputer to run ``wrf.exe``


Python and packages
-------------------
Download and install the Python 2 `Anaconda Python <https://www.continuum.io/downloads>`_ distribution for your platform.  We recommend an installation into the users home directory.

::
  wget https://repo.continuum.io/archive/Anaconda2-5.3.0-Linux-x86_64.sh
  chmod +x Anaconda2-5.3.0-Linux-x86_64.sh
  ./Anaconda2-5.3.0-Linux-x86_64.sh

Install pre-requisites: 

::

  conda install basemap netcdf4 pyproj paramiko dill
  conda install -c conda-forge simplekml pygrib f90nml
  conda install -c anaconda numpy
  pip install MesoPy
  pip install python-cmr


Add your environment (~/.custom.csh) the following line
  setenv PROJ_LIB "$HOME/anaconda2/share/proj"
or if you are using bash to ~/.profile or ~/.bashrc
  setenv PROJ_LIB "$HOME/anaconda2/share/proj"

Note that ``conda`` and ``pip`` are package managers available in the Anaconda Python distribution.

wrfxpy
------

Next, clone the *wrfxpy* code:

::
  
  git clone https://github.com/openwfm/wrfxpy.git

configuration
-------------

And finally, a ``etc/conf.json`` file must be created with the keys discussed below.  A template file ``etc/conf.json.initial`` is provided as a starting point.

Configure the system directories, ``WPS/WRF-SFIRE`` locations and workspace locations by editing the following keys:

::

  "workspace_path": "wksp"
  "wps_install_path": "path/to/WPS"
  "wrf_install_path": "path/to/WRF"
  "sys_install_path": "/path/to/wrfxpy"
  "wps_geog_path" : "/path/to/wps-geogrid"

Optionally, the *wrfxpy* installation can be connected to a visualization server `wrfxweb <https://github.com/vejmelkam/wrfxweb>`_.  The following keys are all optional (and only used if the postprocessed results of simulations are uploaded).

::

  "shuttle_ssh_key": "path/to/your/priv_ssh_key/to/remote/host",
  "shuttle_remote_user" : "remote_username",
  "shuttle_remote_host" : "remote_hostname",
  "shuttle_remote_root" : "remote directory for output storage"

This concludes the ``etc/conf.json`` file.

Next, *wrfxpy* needs to know how jobs are submitted on your cluster.  Create an entry for your cluster, here we use ``speedy`` as an example::

  {
    "speedy" : {
      "qsub_cmd" : "qsub",
      "qsub_script" : "etc/qsub/speedy.sub"
    }
  }

And then the file ``etc/qsub/speedy.sub`` should contain a submission script template, that makes use of the following variables supplied by *wrfxpy* based on job configuration:

* ``%(nodes)d`` the number of nodes requested
* ``%(ppn)d`` the number of processors per node requested
* ``%(wall_time_hrs)d`` the number of hours requested
* ``%(exec_path)d`` the path to the wrf.exe that should be executed
* ``%(cwd)d`` the job working directory
* ``%(task_id)d`` a task id that can be used to identify the job
* ``%(np)d`` the total number of processes requested, equals ``nodes`` x ``ppn``

Note that not all keys need to be used, as shown in the ``speedy`` example::

  #$ -S /bin/bash
  #$ -N %(task_id)s
  #$ -wd %(cwd)s
  #$ -l h_rt=%(wall_time_hrs)d:00:00
  #$ -pe mpich %(np)d
  mpirun_rsh -np %(np)d -hostfile $TMPDIR/machines %(exec_path)s

The script template should be derived from a working submission script.


.. attention::
  You are now ready for your first fire simulation, continue with :doc:`quickstart`.


  


Custom installation
===================

If Anaconda python is not practical, a different python distribution can be used.  Below is a list of packages the system requires:

* `Python 2.7+ <https://www.python.org/download/releases/2.7/>`_
* `Basemap <http://matplotlib.org/basemap/>`_  to render the rasters
* `simplekml <https://simplekml.readthedocs.org/en/latest/>`_ to build KMZ files
* `f90nml <https://pypi.python.org/pypi/f90nml>`_ to manipulate Fortran namelists
* `pyproj <https://pypi.python.org/pypi/pyproj>`_ to place domains dynamically in LCC projection
* `paramiko <https://pypi.python.org/pypi/paramiko>`_ to communicate over SSH with remote hosts
* `netCDF4 <https://pypi.python.org/pypi/netCDF4>`_ to manipulate WPS and WRF files
* `MesoPy <https://pypi.python.org/pypi/MesoPy>`_ to retrieve fuel moisture observations from Mesowest

*wrfxpy* is installed by cloning a GitHub repository

::

  git clone https://github.com/openwfm/wrfxpy.git

Configure *wrfxpy* by editing ``etc/conf.json`` as above and then continue with :doc:`quickstart`.

