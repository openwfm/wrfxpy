
_var_wisdom = {
    'T2' : {
        'name' : 'temperature at 2m',
        'native_unit' : 'K',
        'colorbar_units' : ['C', 'F'],
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d,t: d.variables['T2'][t,:,:]
      },
    'RH_FIRE' : {
        'name' : 'relative humidity',
        'native_unit' : '-',
        'colorbar_units' : ['-'],
        'colormap' : 'jet_r',
        'scale' : [0.0, 1.0],
        'retrieve_as' : lambda d,t: d.variables['RH_FIRE'][t,:,:]
      },
     'FIRE_HFX' : {
        'name' : 'fire heat flux',
        'native_unit' : 'W/m^2',
        'colorbar_units' : ['W/m^2'],
        'colormap' : 'jet',
        'scale' : 'original',
        'retrieve_as' : lambda d,t: d.variables['FIRE_HFX'][t,:,:]
      },
    'F_ROS' : {
       'name' : 'fire spread rate',
       'native_unit' : 'm/s',
       'colorbar_units' : ['m/s', 'ft/s'],
       'colormap' : 'jet',
       'scale' : [0.0, 2.0],
       'retrieve_as' : lambda d,t: d.variables['F_ROS'][t,:,:]
      },
    'F_INT' : {
       'name' : 'fireline intensity',
       'native_unit' : 'J/m/s^2',
       'colorbar_units' : ['J/m/s^2'],
       'colormap' : 'jet',
       'scale' : 'original',
       'retrieve_as' : lambda d,t: d.variables['F_INT'][t,:,:]
    },
    'NFUEL_CAT' : {
       'name' : 'fuel categories',
       'native_unit' : 'fuel_type',
       'colorbar_units' : ['fuel_type'],
       'colormap' : 'jet',
       'scale' : 'original',
       'retrieve_as' : lambda d,t: d.variables['NFUEL_CAT'][t,:,:]
    },
    'ZSF' : {
       'name' : 'terrain height',
       'native_unit' : 'm',
       'colorbar_units' : ['m','ft'],
       'colormap' : 'jet',
       'scale' : 'original',
       'retrieve_as' : lambda d,t: d.variables['ZSF'][t,:,:]
    },
    'FMC_G' : {
       'name' : 'fuel moisture',
       'native_unit' : '-',
       'colorbar_units' : ['-'],
       'colormap' : 'jet_r',
       'scale' : [0.0, 0.5],
       'retrieve_as' : lambda d,t: d.variables['FMC_G'][t,:,:]
    },
    '1HR_FM' : {
       'name' : '1-HR fuel moisture',
       'native_unit' : '-',
       'colorbar_units' : ['-'],
       'colormap' : 'jet_r',
       'scale' : [0.0, 1.0],
       'retrieve_as' : lambda d,t: d.variables['FMC_GC'][t,0,:,:]
    },
    '10HR_FM' : {
       'name' : '10-HR fuel moisture',
       'native_unit' : '-',
       'colorbar_units' : ['-'],
       'colormap' : 'jet_r',
       'scale' : [0.0, 1.0],
       'retrieve_as' : lambda d,t: d.variables['FMC_GC'][t,1,:,:]
    },
    '100HR_FM' : {
       'name' : '100-HR fuel moisture',
       'native_unit' : '-',
       'colorbar_units' : ['-'],
       'colormap' : 'jet_r',
       'scale' : [0.0, 1.0],
       'retrieve_as' : lambda d,t: d.variables['FMC_GC'][t,2,:,:]
    },
    'FMC_EQUI' : {
       'name' : 'fuel moisture equilibrium',
       'native_unit' : '-',
       'colorbar_units' : ['-'],
       'colormap' : 'jet_r',
       'scale' : [0.0, 1.0],
       'retrieve_as' : lambda d,t: d.variables['FMC_EQUI'][t,0,:,:]
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
    ('ft',  'm') : lambda x: x / 3.2808399
}



def get_wisdom(var_name):
    """Return rendering wisdom for the variable <var_name>."""
    return _var_wisdom[var_name]

def get_wisdom_variables():
    """Return the variables for which wisdom is available."""
    return _var_wisdom.keys()

def convert_value(unit_from, unit_to, value):
    # handle the simple case
    if unit_from == unit_to:
        return value

    func = _units_wisdom.get((unit_from, unit_to))
    if func is None:
        return None
    else:
        return func(value)


