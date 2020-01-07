from __future__ import absolute_import
import numpy as np
import logging

from vis.vis_utils import interpolate2height, height8p, height8p_terrain, \
      u8p, v8p, cloud_to_level_hPa, smoke_to_height_terrain, density, print_stats, \
      smoke_concentration
from six.moves import range

smoke_threshold_int = 300
smoke_threshold = 100
smoke_integrated_unit = 'g/m^2'
smoke_integrated_transparent = 1e-2
smoke_concentration_scale = 400
smoke_concentration_transparent=10

def smoke_to_height_terrain_u(var,d,t,h):
      v=convert_value('ug/m^2', smoke_integrated_unit,smoke_to_height_terrain(d,t,h))
      print_stats(var,v,smoke_integrated_unit)
      return v

def plume_center(d,t):
      """
      Compute plume center of mass
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      tr = d.variables['tr17_1'][t,:,:,:]
      smoke_int = np.sum(tr, axis = 0)
      z = height8p(d,t)
      h = np.sum(z * tr, axis = 0)
      h[smoke_int <= smoke_threshold_int] = 0
      smoke_int[smoke_int <= smoke_threshold_int] = 1
      return h/smoke_int

def plume_height(d,t):
      """
      Compute plume height
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      z =  height8p(d,t)
      tr = d.variables['tr17_1'][t,:,:,:]
      h = np.zeros(tr.shape[1:])
      for i in range(0, tr.shape[2]):
          for j in range(0, tr.shape[1]):
              for k in range(tr.shape[0]-1, -1, -1):
                   if tr[k,j,i] > smoke_threshold:
                        h[j,i] = z[k,j,i]
                        break
      return h

def smoke_at_height_terrain_ft(varname,d,t,level_ft):
      return smoke_at_height_terrain(varname,d,t,convert_value('ft','m',level_ft))

def interpolate2height_terrain(d,t,var,level):
      return interpolate2height(var,height8p_terrain(d,t),level)

def smoke_at_height_terrain(varname,d,t,level):
      s = interpolate2height_terrain(d,t,smoke_concentration(d,t),level)
      print_stats(varname,s,'ug/m^3')
      return s

def u8p_m(d,t,level):
       return interpolate2height(u8p(d,t),height8p_terrain(d,t),level)

def v8p_m(d,t,level):
       return interpolate2height(v8p(d,t),height8p_terrain(d,t),level)

def u8p_ft(d,t,level_ft):
       return u8p_m(d,t,convert_value('ft','m',level_ft))

def v8p_ft(d,t,level_ft):
       return v8p_m(d,t,convert_value('ft','m',level_ft))

def is_windvec(name):
       return name in ['WINDVEC1000FT','WINDVEC4000FT','WINDVEC6000FT','WINDVEC']

def is_fire_var(name):
       return name in ['FGRNHFX','FIRE_AREA','FLINEINT','FIRE_HFX','F_ROS','F_INT','NFUEL_CAT','ZSF','FMC_G']

_var_wisdom = {
     'CLOUDTO700HPA' : {
        'name' : 'Cloud up to 700hPa',
        'native_unit' : 'kg/m^2',
        'colorbar' : 'kg/m^2',
        'colorbar' : None,
        'colormap' : 'gray_r',
        'transparent_values' : [-np.inf,0.00001],
        'scale' : 'original',
        'retrieve_as' : lambda d,t: cloud_to_level_hPa(d,t,700),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'CLOUD700TO400HPA' : {
        'name' : 'Cloud 700hPa to 400hPa',
        'native_unit' : 'kg/m^2',
        'colorbar' : 'kg/m^2',
        'colorbar' : None,
        'colormap' : 'gray_r',
        'transparent_values' : [-np.inf,0.00001],
        'scale' : 'original',
        'retrieve_as' : lambda d,t: cloud_to_level_hPa(d,t,400) - cloud_to_level_hPa(d,t,700),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'CLOUDABOVE400HPA' : {
        'name' : 'Cloud above 400hPa',
        'native_unit' : 'kg/m^2',
        'colorbar' : 'kg/m^2',
        'colorbar' : None,
        'colormap' : 'gray_r',
        'transparent_values' : [-np.inf,0.00001],
        'scale' : 'original',
        'retrieve_as' : lambda d,t: cloud_to_level_hPa(d,t,0) - cloud_to_level_hPa(d,t,400),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'SMOKE1000FT' : {
        'name' : 'Smoke 1000ft above terrain',
        'native_unit' : 'ug/m^3',
        'colorbar' : 'ug/m^3',
        'colormap' : 'rainbow',
        'transparent_values' : [-np.inf,10],
        'scale' : [0, 500],
        'retrieve_as' : lambda d,t: smoke_at_height_terrain_ft('SMOKE1000FT',d,t,1000),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'SMOKE4000FT' : {
        'name' : 'Smoke 4000ft above terrain',
        'native_unit' : 'ug/m^3',
        'colorbar' : 'ug/m^3',
        'colormap' : 'rainbow',
        'transparent_values' : [-np.inf,10],
        'scale' : [0, 500],
        'retrieve_as' : lambda d,t: smoke_at_height_terrain_ft('SMOKE4000FT',d,t,4000),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'SMOKE6000FT' : {
        'name' : 'Smoke 6000ft above terrain',
        'native_unit' : 'ug/m^3',
        'colorbar' : 'ug/m^3',
        'colormap' : 'rainbow',
        'transparent_values' : [-np.inf,10],
        'scale' : [0, 500],
        'retrieve_as' : lambda d,t: smoke_at_height_terrain_ft('SMOKE6000FT',d,t,6000),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'WINDSPD1000FT' : {
        'name' : 'wind speed at 1000ft',
        'native_unit' : 'm/s',
        'colorbar' : 'm/s',
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d, t: np.sqrt(u8p_ft(d,t,1000)**2.0 + v8p_ft(d,t,1000)**2.0),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
      },
    'WINDVEC1000FT' : {
        'name' : 'wind speed at 1000ft',
        'components' : [ 'U1000FT', 'V1000FT' ],
        'native_unit' : 'm/s',
        'scale' : 'original',
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    'U1000FT' : {
        'name' : 'longitudinal wind component 1000ft',
        'retrieve_as' : lambda d, t: u8p_ft(d,t,1000),
    },
    'V1000FT' : {
        'name' : 'latitudinal wind component 1000ft',
        'retrieve_as' : lambda d, t: v8p_ft(d,t,1000),
    },
     'WINDSPD4000FT' : {
        'name' : 'wind speed at 4000ft',
        'native_unit' : 'm/s',
        'colorbar' : 'm/s',
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d, t: np.sqrt(u8p_ft(d,t,4000)**2.0 + v8p_ft(d,t,4000)**2.0),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
      },
    'WINDVEC4000FT' : {
        'name' : 'wind speed at 4000ft',
        'components' : [ 'U4000FT', 'V4000FT' ],
        'native_unit' : 'm/s',
        'scale' : 'original',
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    'U4000FT' : {
        'name' : 'longitudinal wind component 4000ft',
        'retrieve_as' : lambda d, t: u8p_ft(d,t,4000),
    },
    'V4000FT' : {
        'name' : 'latitudinal wind component 4000ft',
        'retrieve_as' : lambda d, t: v8p_ft(d,t,4000),
    },
     'WINDSPD6000FT' : {
        'name' : 'wind speed at 6000ft',
        'native_unit' : 'm/s',
        'colorbar' : 'm/s',
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d, t: np.sqrt(u8p_ft(d,t,6000)**2.0 + v8p_ft(d,t,6000)**2.0),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
      },
    'WINDVEC6000FT' : {
        'name' : 'wind speed at 6000ft',
        'components' : [ 'U6000FT', 'V6000FT' ],
        'native_unit' : 'm/s',
        'scale' : 'original',
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    'U6000FT' : {
        'name' : 'longitudinal wind component 6000ft',
        'retrieve_as' : lambda d, t: u8p_ft(d,t,6000),
    },
    'V6000FT' : {
        'name' : 'latitudinal wind component 6000ft',
        'retrieve_as' : lambda d, t: v8p_ft(d,t,6000),
    },
     'PLUME_HEIGHT' : {
        'name' : 'plume height',
        'native_unit' : 'm',
        'colorbar' : 'm',
        'colormap' : 'jet',
        'transparent_values' : [-np.inf, 50],
        'scale' : [0, 8000],
        'retrieve_as' : lambda d,t: plume_height(d,t),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'PLUME_CENTER' : {
        'name' : 'plume height',
        'native_unit' : 'm',
        'colorbar' : 'm',
        'colormap' : 'jet',
        'transparent_values' : [-np.inf, 50],
        'scale' : [0, 8000],
        'retrieve_as' : lambda d,t: plume_center(d,t),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'FGRNHFX' : {
        'name' : 'Ground level heat flux [log]',
        'native_unit' : 'W/m^2',
        'colorbar' : 'W/m^2',
        'colormap' : 'jet',
        'transparent_values' : [-np.inf, 1],
        'scale' : [0, 6],
        'retrieve_as' : lambda d,t: np.ma.filled(np.ma.log10(np.ma.masked_less_equal(d.variables['FGRNHFX'][t,:,:], 0)), 0),
        'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:]),
      },
     'SMOKE_INT' : {
        'name' : 'vertically integrated smoke',
        'native_unit' : smoke_integrated_unit,
        'colorbar' : None,
        'colormap' : 'gray_r',
        'transparent_values' : [-np.inf, smoke_integrated_transparent],
        'scale' : 'original',
        'retrieve_as' : lambda d,t: smoke_to_height_terrain_u('SMOKE_INT',d,t,100000),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'SMOKETO10M' : {
        'name' : 'vertically integrated smoke to 10m',
        'native_unit' : smoke_integrated_unit,
        'colorbar' : None,
        'colormap' : 'gray_r',
        'transparent_values' : [-np.inf, smoke_integrated_transparent],
        'scale' : 'original',
        'retrieve_as' : lambda d,t: smoke_to_height_terrain_u('SMOKETO10M',d,t,10),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
     },
     'PM25_INT' : {
        'name' : 'vert integrated PM2.5',
        'native_unit' : smoke_integrated_unit,
        'colorbar' : smoke_integrated_unit,
        'colormap' : 'rainbow',
        'transparent_values' : [-np.inf, smoke_integrated_transparent],
        'scale' : [0, 5],
        'retrieve_as' : lambda d,t: smoke_to_height_terrain_u('PM25_INT',d,t,100000),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
     'PM25_SFC_CHEM' : {
        'name' : 'surface PM2.5 from WRF-CHEM',
        'native_unit' : smoke_integrated_unit,
        'colorbar' : 'ug/m^2',
        'colormap' : 'rainbow',
        'transparent_values' : [-np.inf, 10],
        'scale' : [20, 90],
        'retrieve_as' : lambda d,t: d.variables['PM2_5_DRY'][t,0,:,:],
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
    'T2' : {
        'name' : 'temperature at 2m',
        'native_unit' : 'K',
        'colorbar' : 'C',
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d,t: d.variables['T2'][t,:,:],
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:]),
      },
    'FIRE_AREA' : {
        'name' : 'fire area',
        'native_unit' : '-',
        'colorbar' : None,
        'colormap' : 'hot_r',
        'transparent_values' : [-np.inf, 0],
        'scale' : [0.0, 1.0],
        'retrieve_as' : lambda d,t: d.variables['FIRE_AREA'][t,:,:],
        'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
      },
    'FLINEINT' : {
        'name' : 'fireline intensity',
        'native_unit' : '-',
        'colorbar' : '-',
        'colormap' : 'jet',
        'transparent_values' : [-np.inf, 1],
        'scale' : 'original',
        'retrieve_as' : lambda d,t: np.ma.filled(np.ma.log10(np.ma.masked_less_equal(d.variables['FLINEINT'][t,:,:], 0)), 0),
        'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
      },
     'RH_FIRE' : {
        'name' : 'relative humidity',
        'native_unit' : '-',
        'colorbar' : '-',
        'colormap' : 'viridis_r',
        'scale' : [0.0, 1.0],
        'retrieve_as' : lambda d,t: d.variables['RH_FIRE'][t,:,:],
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
      },
     'FIRE_HFX' : {
        'name' : 'fire heat flux',
        'native_unit' : 'W/m^2',
        'colorbar' : 'W/m^2',
        'colormap' : 'jet',
        'scale' : 'original',
        'transparent_values' : [-np.inf, 0],
        'retrieve_as' : lambda d,t: d.variables['FIRE_HFX'][t,:,:],
        'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
      },
    'F_ROS' : {
        'name' : 'fire spread rate',
        'native_unit' : 'm/s',
        'colorbar' : 'm/s',
        'colormap' : 'jet',
        'scale' : [0.0, 2.0],
        'retrieve_as' : lambda d,t: d.variables['F_ROS'][t,:,:],
        'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
      },
    'PSFC' : {
        'name' : 'surface pressure',
        'native_unit' : 'Pa',
        'colorbar' : 'Pa',
        'colormap' : 'rainbow',
        'scale' : 'original',
        'retrieve_as' : lambda d, t: d.variables['PSFC'][t,:,:],
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
      },
     'WINDSPD' : {
        'name' : 'wind speed',
        'native_unit' : 'm/s',
        'colorbar' : 'm/s',
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d, t: np.sqrt(d.variables['U10'][t,:,:]**2.0 + d.variables['V10'][t,:,:]**2.0),
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
      },
    'WINDVEC' : {
        'name' : 'wind speed',
        'components' : [ 'U10', 'V10' ],
        'native_unit' : 'm/s',
        'scale' : 'original',
        'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    'U10' : {
        'name' : 'longitudinal wind component',
        'retrieve_as' : lambda d, t: d.variables['U10'][t,:,:],
    },
    'V10' : {
        'name' : 'latitudinal wind component',
        'retrieve_as' : lambda d, t: d.variables['V10'][t,:,:],
    },
    'F_INT' : {
        'name' : 'fireline intensity',
        'native_unit' : 'J/m/s^2',
        'colorbarct' : 'J/m/s^2',
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d,t: d.variables['F_INT'][t,:,:],
        'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
      },
    'NFUEL_CAT' : {
       'name' : 'fuel categories',
       'native_unit' : 'fuel_type',
       'colorbar' : 'fuel_type',
       'colormap' : 'Dark2',
       'scale' : 'original',
       'retrieve_as' : lambda d,t: d.variables['NFUEL_CAT'][t,:,:],
       'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
    },
    'ZSF' : {
       'name' : 'terrain height',
       'native_unit' : 'm',
       'colorbar' : 'm',
       'colormap' : 'terrain',
       'scale' : 'original',
       'retrieve_as' : lambda d,t: d.variables['ZSF'][t,:,:],
       'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
    },
    'FMC_G' : {
       'name' : 'fuel moisture',
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'gist_earth_r',
       'scale' : [0.0, 0.5],
       'retrieve_as' : lambda d,t: d.variables['FMC_G'][t,:,:],
       'grid' : lambda d: (d.variables['FXLAT'][0,:,:], d.variables['FXLONG'][0,:,:])
    },
    '1HR_FM' : {
       'name' : '1-HR fuel moisture',
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'jet_r',
       'scale' : [0.0, 0.5],
       'retrieve_as' : lambda d,t: d.variables['FMC_GC'][t,0,:,:],
       'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    '10HR_FM' : {
       'name' : '10-HR fuel moisture',
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'jet_r',
       'scale' : [0.0, 0.5],
       'retrieve_as' : lambda d,t: d.variables['FMC_GC'][t,1,:,:],
       'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    '100HR_FM' : {
       'name' : '100-HR fuel moisture',
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'jet_r',
       'scale' : [0.0, 0.5],
       'retrieve_as' : lambda d,t: d.variables['FMC_GC'][t,2,:,:],
       'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    'FMC_EQUI' : {
       'name' : 'fuel moisture equilibrium',
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'jet_r',
       'scale' : [0.0, 1.0],
       'retrieve_as' : lambda d,t: d.variables['FMC_EQUI'][t,0,:,:],
       'grid' : lambda d: (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
    },
    'TERRA_AF' : {
       'name' : 'MODIS Terra Active Fires satellite data',
       'source' : 'Terra',
       'native_unit' : '-',
       'colorbar' : '-',
       'colormap' : 'discrete',
       'scale' : 'discrete',
       'retrieve_as' : lambda d : d.select('fire mask').get(),
       'grid' : lambda d: (d.select('Latitude').get(), d.select('Longitude').get())
    },
    'AQUA_AF' : {
        'name' : 'MODIS Aqua Active Fires satellite data',
       	'source' : 'Aqua',
        'native_unit' : '-',
        'colorbar' : '-',
        'colormap' : 'discrete',
        'scale' : 'discrete',
        'retrieve_as' : lambda d : d.select('fire mask').get(),
        'grid' : lambda d : (d.select('Latitude').get(), d.select('Longitude').get())
    },
    'SNPP_AF' : {
        'name' : 'VIIRS S-NPP Active Fires satellite data',
       	'source' : 'SNPP',
        'native_unit' : '-',
        'colorbar' : '-',
        'colormap' : 'discrete',
        'scale' : 'discrete',
        'retrieve_as' : lambda d : d.variables['fire mask'][:],
        'grid' : lambda d : (d['HDFEOS']['SWATHS']['VNP_750M_GEOLOCATION']['Geolocation Fields']['Latitude'], d['HDFEOS']['SWATHS']['VNP_750M_GEOLOCATION']['Geolocation Fields']['Longitude'])
    }
}

# contains functions to transform values from one unit to another in a simple format.
# it's a dictionary with keys in the form (from_unit, to_unit) and the value is a lambda
# that maps the value from <from_unit> to <to_unit>.
_units_wisdom = {
    ('K',   'C') : lambda x: x - 273.15,
    ('K',   'F') : lambda x: 9.0 / 5.0 * (x - 273.15) + 32,
    ('m/s', 'ft/s') : lambda x: 3.2808399 * x,
    ('m',   'ft') : lambda x: 3.2808399 * x,
    ('ft/s','m/s') : lambda x: x / 3.2808399,
    ('ft',  'm') : lambda x: x / 3.2808399,
    ('ug/m^2', 'g/m^2') : lambda x: 1e-6 * x
}



def get_wisdom(var_name):
    """Return rendering wisdom for the variable <var_name>."""
    return _var_wisdom[var_name]

def get_wisdom_variables():
    """Return the variables for which wisdom is available."""
    return list(_var_wisdom.keys())

def convert_value(unit_from, unit_to, value):
    # handle the simple case
    if unit_from == unit_to:
        return value

    func = _units_wisdom.get((unit_from, unit_to))
    if func is None:
        return value
    else:
        return func(value)


