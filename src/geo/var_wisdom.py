from utils import Dict

_var_wisdom = {'NFUEL_CAT' : {
        'units': "fuel category",
        'description': "Anderson 13 fire behavior categories",
        'type': "categorical", 
        'category_range': [0,14],
        'missing_value': 99,       
        'fill' : Dict({
                    range(90,100): 14,
                    100: 'nearest'
                }),
        'scale': 1., # scale the array to be integer (default: depending on bits, not really good option for int array)
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:nearest_neighbor+average_16pt+search'             
    },
    'ZSF' : {
        'units': "meters",
        'description': "National Elevation Dataset 1/3 arcsecond resolution",
        'type': "continuous",
        'missing_value': 0.0,
        'scale': 1.,
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option': 'smth-desmth_special; smooth_passes=1'
    },
    'FMC_GC' : {
        'units': "1",
        'description': "1h, 10h, 100h fuel moisture",
        'type': "continuous",
        'tile_bdr': 0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt'
    },
    'FMEP' : {
        'units': "1",
        'description': "fuel moisture drying/wetting and rain equilibrium adjustments",
        'type': "continuous",
        'tile_bdr': 0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt'
    },
    'XLONG' : {
        'units': "degrees",
        'description': "Longitude in the domain",
        'type': "continuous",
        'missing_value': 0.0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option': 'smth-desmth_special; smooth_passes=1'
    },
    'XLAT' : {
        'units': "degrees",
        'description': "Latitude in the domain",
        'type': "continuous",
        'missing_value': 0.0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option': 'smth-desmth_special; smooth_passes=1'
    },
    'XI' : {
        'units': "1",
        'description': "Indexes x-direction for testing WPS",
        'type': "continuous",
        'missing_value': 0.0,
        'scale': 1.,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt'
    },
    'YI' : {
        'units' : "1",
        'description': "Indexes y-direction for testing WPS",
        'type': "continuous",
        'missing_value': 0.0,
        'scale': 1.,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt'
    }
}

def get_wisdom(var_name):
    """Return rendering wisdom for the variable <var_name>."""
    return _var_wisdom[var_name]

def get_wisdom_variables():
    """Return the variables for which wisdom is available."""
    return list(_var_wisdom.keys())
