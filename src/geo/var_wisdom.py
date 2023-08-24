from utils import Dict
import os.path as osp
import json

_var_wisdom = {
    'NFUEL_CAT_13' : {
        'var_name': 'NFUEL_CAT',
        'units': "fuel category",
        'description': "Anderson 13 fire behavior categories",
        'type': "categorical", 
        'category_range': [0,14],
        'fill_missing': 14,       
        'fill' : Dict({
            range(15,91): 14, 92: 14, range(94,100): 14,
            -9999: 14, (0,91,93): 'nearest'
        }),
        'scale': 1., # scale the array to be integer (default: depending on bits, not really good option for int array)
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:nearest_neighbor+average_16pt+search',
        'subgrid': 'yes',
        'add_opts': {'dominant_only': 'NFUEL_CAT',
                    'z_dim_name': 'fuel_cat',
                    'halt_on_missing': 'no'}             
    },
    'NFUEL_CAT_40' : {
        'var_name': 'NFUEL_CAT',
        'units': 'fuel category',
        'description': '40 Scott and Burgan fire behavior categories',
        'type': 'categorical', 
        'category_range': [0,41],
        'fill_missing': 41,       
        'fill' : Dict({
            101: 1, 102: 2, 103: 3, 104: 4, 105: 5, 106: 6, 107: 7, 108: 8,
            109: 9, 121: 10, 122: 11, 123: 12, 124: 13, 141: 14, 142: 15,
            143: 16, 144: 17, 145: 18, 146: 19, 147: 20, 148: 21, 149: 22,
            161: 23, 162: 24, 163: 25, 164: 26, 165: 27, 181: 28, 182: 29,
            183: 30, 184: 31, 185: 32, 186: 33, 187: 34, 188: 35, 189: 36,
            201: 37, 202: 38, 203: 39, 204: 40, range(42,91): 41, 92: 41,
            range(94,100): 41, -9999: 41, (0,91,93): 'nearest'
        }),
        'scale': 1., # scale the array to be integer (default: depending on bits, not really good option for int array)
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:nearest_neighbor+average_16pt+search',
        'subgrid': 'yes',
        'add_opts': {'dominant_only': 'NFUEL_CAT',
                    'z_dim_name': 'fuel_cat',
                    'halt_on_missing': 'no'}             
    },
    'NTREE_CAT' : {
        'var_name': 'NTREE_CAT',
        'units': 'tree type',
        'description': 'Existing 200 vegetation type categories',
        'type': 'categorical',
        'category_range': [0,6],
        'fill_missing': 0,
        'scale': 1., # scale the array to be integer (default: depending on bits, not really good option for int array)
        'fill' : 'etc/vtables/fill_ntree_cat.csv',
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:nearest_neighbor+average_16pt+search',
        'subgrid': 'yes',
        'add_opts': {'dominant_only': 'NTREE_CAT',
                    'z_dim_name': 'tree_cat',
                    'halt_on_missing': 'no'}
    },
    'ZSF' : {
        'var_name': 'ZSF',
        'units': 'meters',
        'description': 'National Elevation Dataset 1/3 arcsecond resolution',
        'type': 'continuous',
        'fill_missing': 0,
        'scale': 1.,
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option': 'smth-desmth_special; smooth_passes=1',
        'subgrid': 'yes',
        'add_opts': {'df_dx': 'DZDXF',
                    'df_dy': 'DZDYF',
                    'halt_on_missing': 'no'}
    },
    'CAN_TOP' : {
        'var_name': 'CAN_TOP',
        'units': 'meters',
        'description': 'Forest Canopy Height from LANDFIRE',
        'type': 'continuous',
        'signed': 'yes',
        'unit_scale': 0.1, # scale the array to be in meters (original data is in decimeters)
        'bits': 16,
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'subgrid': 'yes',
    },
    'CAN_BOT' : {
        'var_name': 'CAN_BOTTOM',
        'units': 'meters',
        'description': 'Forest Canopy Base Height from LANDFIRE',
        'type': 'continuous',
        'signed': 'yes',
        'unit_scale': 0.1, # scale the array to be in meters (original data is in decimeters)
        'bits': 16,
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'subgrid': 'yes',
    },
    'CAN_BD' : {
        'var_name': 'CAN_BULK_DENSITY',
        'units': 'kg m-3',
        'description': 'Forest Canopy Bulk Density from LANDFIRE',
        'type': 'continuous',
        'signed': 'yes',
        'unit_scale': 0.01, # scale the array to be in kg/m3 (original data is multiplied by 100)
        'bits': 16,
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'subgrid': 'yes',
    },
    'CAN_COV' : {
        'var_name': 'CAN_COVER',
        'units': '%',
        'description': 'Forest Canopy Cover from LANDFIRE (%)',
        'type': 'continuous',
        'signed': 'yes',
        'bits': 16,
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'subgrid': 'yes',
    },
    'FMC_GC' : {
        'var_name': 'FMC_GC',
        'units': 'g/g',
        'description': '1h, 10h, 100h fuel moisture',
        'type': 'continuous',
        'tile_bdr': 0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt'
    },
    'FMEP' : {
        'var_name': 'FMEP',
        'units': '1',
        'description': 'fuel moisture drying/wetting and rain equilibrium adjustments',
        'type': 'continuous',
        'tile_bdr': 0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt'
    },
    'XLONG' : {
        'var_name': 'XLONG',
        'units': 'degrees',
        'description': 'Longitude in the domain',
        'type': 'continuous',
        'fill_missing': 0.0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option': 'smth-desmth_special; smooth_passes=1',
        'subgrid': 'no'
    },
    'XLAT' : {
        'var_name': 'XLAT',
        'units': 'degrees',
        'description': 'Latitude in the domain',
        'type': 'continuous',
        'fill_missing': 0.0,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'smooth_option': 'smth-desmth_special; smooth_passes=1',
        'subgrid': 'no'
    },
    'XI' : {
        'var_name': 'XI',
        'units': '1',
        'description': 'Indexes x-direction for testing WPS',
        'type': 'continuous',
        'fill_missing': 0.0,
        'scale': 1.,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'subgrid': 'yes'
    },
    'YI' : {
        'var_name': 'YI',
        'units' : '1',
        'description': 'Indexes y-direction for testing WPS',
        'type': 'continuous',
        'fill_missing': 0.0,
        'scale': 1.,
        'signed': 'yes',
        'interp_option': 'default:average_gcell(4.0)+four_pt+average_4pt',
        'subgrid': 'yes'
    }
}

class VarWisdom(dict):
    """
    A dictionary with the variable wisdom.
    """
    def __init__(self, d={}):
        """
        Updates itself with d.
        """
        self = _var_wisdom.copy()
        self.update(d)
        custom_path = 'etc/vtables/geo_vars_custom.json'
        if osp.exists(custom_path):
            self.update(json.load(open(custom_path)))

def get_wisdom(var_name):
    """Return rendering wisdom for the variable <var_name>."""
    return VarWisdom()[var_name]

def get_wisdom_variables():
    """Return the variables for which wisdom is available."""
    return list(VarWisdom().keys())
