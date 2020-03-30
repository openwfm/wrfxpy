
_var_wisdom = {'NFUEL_CAT' : {
        'units' : "fuel category",
        'description' : "Anderson 13 fire behavior categories",
        'type' : "categorical", 
        'category_range' : [0,99],
        'missing_value' : 99, 
        # missing and generated from nearest neighbours
        'missing' : [99],       
        # no value filled with from fill     
        'novalue' : range(90,99),  
        # value to fill the no value data       
        'fill_missing' : 14,
        'scale': 1., # scale the array to be integer (default: depending on bits, not really good option for int array)
        'signed': 'yes',
        'bits': 16,
        'interp_option' : 'default:nearest_neighbor+average_16pt+search'             
    },
    'ZSF' : {
        'units' : "meters",
        'description' : "National Elevation Dataset 1/3 arcsecond resolution",
        'type' : "continuous",
        'missing_value' : 0.0,
        'scale': 1.,
        'signed':'yes',
        'bits': 16,
        'interp_option' : 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option' : 'smth-desmth_special; smooth_passes=1'
    },
    'FMC_GC' : {
        'units' : "1",
        'description' : "1h, 10h, 100h fuel moisture",
        'type' : "continuous",
        'tile_bdr' : 0,
        'signed':'yes',
        'interp_option' : 'default:average_gcell(4.0)+four_pt+average_4pt'
    },
    'FMEP' : {
        'units' : "1",
        'description' : "fuel moisture drying/wetting and rain equilibrium adjustments",
        'type' : "continuous",
        'tile_bdr' : 0,
        'signed': 'yes',
        'interp_option' : 'default:average_gcell(4.0)+four_pt+average_4pt'
    }
}

def get_wisdom(var_name):
    """Return rendering wisdom for the variable <var_name>."""
    return _var_wisdom[var_name]

def get_wisdom_variables():
    """Return the variables for which wisdom is available."""
    return list(_var_wisdom.keys())