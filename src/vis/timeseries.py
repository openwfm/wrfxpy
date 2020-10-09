import logging
import subprocess
import pandas as pd
import numpy as np
import os.path as osp
from vis.var_wisdom import get_wisdom, is_fire_var, strip_end
from clamp2mesh import nearest_idx

class Timeseries(object):
    """
    Timeseries of WRF data.
    """

    def __init__(self, output_path, prod_name, tslist, num_doms):
        """
        Initialize timeseries with output parameters.

        :param output_path: path where timeseries files are stored
        :param prod_name: name of manifest json file and prefix of all output files
        :param tslist: dictionary with time series information
        :param num_doms: number of domains
        """
        logging.info("Timeseries: output_path=%s prod_name=%s" % (output_path, prod_name))
        self.output_path = output_path
        self.product_name = prod_name
        self.stations = tslist['stations'].copy()
        self.variables = {var: get_wisdom(var).copy() for var in tslist['vars']}
        logging.info("Timeseries: stations=%s" % [st['name'] for st in self.stations.values()])
        logging.info("Timeseries: variables=%s" % list(self.variables.keys()))
        self.num_doms = num_doms
        # initialize the CSV files for each station
        self.initialize_stations()

    def initialize_stations(self):
        static_cols = ['station_name','station_lon','station_lat',
                       'datetime', 'domain', 'grid_i', 'grid_j', 
                       'grid_lon', 'grid_lat', 'grid_fire_i',
                       'grid_fire_j', 'grid_lon_fire', 'grid_lat_fire']
        var_cols = list(self.variables.keys())
        for st in self.stations.keys():
            self.stations[st].update({'local_path': {}})
            for dom_id in range(1,self.num_doms+1):
                st_path = osp.join(self.output_path, self.product_name + '-%02d-' % dom_id + st + '.csv')
                self.stations[st]['local_path'].update({str(dom_id): st_path})
                cols = static_cols + var_cols 
                df = pd.DataFrame({c: [] for c in cols}) 
                df.to_csv(st_path,index = False)

    def write_timestep(self,d,dom_id,tndx,ts_esmf):
        ts_paths = []
        logging.info('write_timestep: time series at time %s and domain %d' % (ts_esmf,dom_id))
        lats,lons = (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
        if 'FXLONG' in d.variables:
            lats_fire,lons_fire = (d.variables['XLAT'][0,:,:], d.variables['XLONG'][0,:,:])
            fm,fn = strip_end(d)
            lats_fire,lons_fire = (lats_fire[:fm,:fn], lons_fire[:fm,:fn]) 
        for k_st,station in self.stations.items():
            idx = nearest_idx(lons,lats,station['lon'],station['lat'])
            timestep = {
                'station_name': station['name'],
                'station_lon': station['lon'],
                'station_lat': station['lat'],
                'datetime': ts_esmf,
                'domain': dom_id,
                'grid_i': idx[0], 
                'grid_j': idx[1], 
                'grid_lon': lons[idx], 
                'grid_lat': lats[idx] 
            }
            if 'FXLONG' in d.variables:
                idx_fire = nearest_idx(lons_fire,lats_fire,station['lon'],station['lat'])
                timestep.update({
                    'grid_fire_i': idx_fire[0], 
                    'grid_fire_j': idx_fire[1], 
                    'grid_fire_lon': lons_fire[idx_fire], 
                    'grid_fire_lat': lats_fire[idx_fire] 
                })
            for k_v,var in self.variables.items():
                array = var['retrieve_as'](d,tndx)
                if is_fire_var(var): 
                    array = array[:fm,:fn]
                    val = array[idx_fire]
                else:
                    val = array[idx]
                timestep.update({k_v: val})
            df = pd.read_csv(station['local_path'][str(dom_id)])
            df = df.append(timestep,ignore_index=True)
            df.to_csv(station['local_path'][str(dom_id)],index=False) 
            ts_paths.append(osp.basename(station['local_path'][str(dom_id)]))
        return ts_paths

