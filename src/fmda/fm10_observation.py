# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


class FM10Observation:
  """
  An observation of an FM-10 value by a RAWS station.
  """

  def __init__(self,tm,lat,lon,elev,obs,ngp):
    """
    Constructs an fm-10 observation from given information.

    :param tm: timestamp
    :param lat: latitude
    :param lon: longitude
    :param elev: elevation (meters)
    :param obs: the observation value
    :param ngp: the nearest grid point (depends on WRF grid)
    """
    self.tm = tm
    self.lon = lon
    self.lat = lat
    self.elevation = elev
    self.obs_val = obs
    self.ngp = ngp


  def get_value(self):
    """
    Return the observation value.
    """
    return self.obs_val


  def get_variance(self):
    """
    Estimate and compute the variance of the measurement of the given field.
    The variance estimate are taken from the Campbell Scientific CS506 manual as a representative
    sample.  We don't generally know from which type of sensor a value is obtained.

    Also, we don't consider the age of the fuel stick (even harder to find out) on which the variance
    depends.

    :return: the estimated variance
    """
    val = self.obs_val
    if val < 0.1:
        return (0.0125/4)**2
    elif val < 0.2:
        return (0.02/4)**2
    elif val < 0.3:
        return (0.034/4)**2
    else:
        return (0.041/4)**2


  def get_elevation(self):
    """
    Return the elevation, where the observation was taken.
    """
    return self.elevation


  def get_location(self):
    """
    Longitude and lattitude of the originating station (shortcut).
    """
    return (self.lat, self.lon)


  def get_nearest_grid_point(self):
    """
    Return the indices that identify the nearest grid point.
    """
    return self.ngp

  
  def get_time(self):
    """
    Return the GMT time of the observations.
    """
    return self.tm
  

  def __str__(self):
    """
    Returns a string representation.
    """
    return "%s loc: [%g,%g] val: %g var: %g" % (str(self.tm),self.lat,self.lon,self.obs_val,self.get_variance())

