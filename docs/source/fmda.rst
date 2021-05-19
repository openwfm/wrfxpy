fmda module
***********

The ``fmda`` module contains code that performs fuel moisture data assimilation.

In particular, the module can:

* open a wrfinput/wrfoutput file, extract meteorological data,
* compute the equilibrium moisture (or run fuel moisture model for a few steps)
* retrieve fuel moisture observations from the `Mesowest <http://mesowest.org>`_ API service
* perform a Kalman filter update on the computed moisture
* store the analysis back in the processed file

.. automodule:: fmda.fuel_moisture_da
   :members:
