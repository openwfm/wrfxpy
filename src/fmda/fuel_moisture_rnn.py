# Helper module for Operational RNN class
# The core RNN functionality is in a distributable package ml_fmda
#     Building model from weights path, cyclical prediction with states,  
#     weight warping for transfer, scaling 3d arrays, ...
# This module is for functionality specific to wrfxpy
# Utility functions plus some hard coded tables and objects
# Jonathon Hirschi, 2026

import copy
import warnings
import numpy as np
import logging
import os.path as osp
from utils import inq, ensure_dir
from wrf.wps_format import WPSFormat
from geo.write_geogrid import write_geogrid_var
from geo.var_wisdom import get_wisdom
import json

from ml_fmda.moisture_rnn import TimeWarpedFuelClassPredictors


# Namelist for converting feature names from ml_fmda project to wrfxpy 
# variable names from data dict build in cycle
RNN_FEATURE_TO_WRFXPY = {
    'solar': 'swdown',
    'wind': 'ws',
    'elev': 'hgt',
    'lat': 'lats',
    'lon': 'lons',
    'solar': 'swdown',
    'temp': 't2',
    'rh': 'rh',
    'pres': 'psfc'
} 

# Core Model Class
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class RNNMoistureModel(TimeWarpedFuelClassPredictors):

    """
    wrfxpy wrapper for the time-warped fuel-moisture predictors.

    The parent class constructs pretrained FM1, FM10, FM100, and FM1000
    OperationalRNNPredictor models. This wrapper adds wrfxpy-specific
    grid and output functionality.
    """

    def __init__(self, params, weights_path, warps):
        super().__init__(params, weights_path, warps)    
    
    # predict() and predict_cycle() inherited from parent

    def states_to_grid(self, grid_shape):
        """
        Convert stored recurrent states to a regular grid.

        :param grid_shape: (ny, nx) grid dimensions.
        :return: array with shape (ny, nx, 4, 2, nunits).
        """
        states = self.get_states()
        states_grid = []

        for fuel_class in self.FUEL_CLASSES:
            state = np.stack(states[fuel_class], axis=1)
            if state.shape[0] != np.prod(grid_shape):
                raise ValueError(
                    f"State batch size ({state.shape[0]}) does not match "
                    f"grid size ({np.prod(grid_shape)})."
                )
            states_grid.append(state.reshape(*grid_shape, *state.shape[1:]))

        return np.stack(states_grid, axis=2)        

    def to_netcdf(self, path, preds, grid_shape, data_vars):
        """
        Store the model in a netCDF file. Modified from FuelMoistureModel: accepts model outputs (FMC_GC) as input, since the RNN output is not stored internally, as only the recurrent state is needed to evolve the system forward.
        
        Expects a single time slice. RNN is compatible with extended time dimension, but not existing cycler

        :param path: the path where to store the model
        :param data_vars: dictionary of additional variables to store for visualization
        :param preds: ndarray of RNN model outputs
        """
        import netCDF4

        if preds.ndim != 3:
            raise ValueError(
                f"Expected preds with 3 dimensions (ny, nx, k), got shape {preds.shape}"
            )
   
        d0, d1, k = preds.shape
        if (d0, d1) != tuple(grid_shape):
            raise ValueError(
                f"Expected preds spatial shape {grid_shape}, got {(d0, d1)}"
            )


        rnn_states = self.states_to_grid((d0, d1))

        # Create file and grid
        ds = netCDF4.Dataset(path, 'w', format='NETCDF4')
        ds.createDimension('fuel_moisture_classes_stag', k)
        ds.createDimension('south_north', d0)
        ds.createDimension('west_east', d1)
        
        # Store model params as global attributes
        ds.warps = json.dumps(self.warps)
        ds.params = json.dumps(self.params)
        ds.weights_path = self.weights_path
        ds.FUEL_CLASSES = json.dumps(self.FUEL_CLASSES)

        # Store FMC predictions and recurrent state
        ncfmc = ds.createVariable('FMC_GC', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag'))
        ncfmc[:,:,:] = preds
        logging.info('rnn_moisture_model.to_netcdf: writing predictions as FMC_GC %s recurrent state as RNN_STATES %s' % (inq(preds),inq(rnn_states)))
        
        ds.createDimension('rnn_state_vars', 2)
        ds.createDimension('rnn_state_units', rnn_states.shape[-1])

        ncrnn_state = ds.createVariable(
            'RNN_STATES',
            'f4',
            ('south_north', 'west_east', 'fuel_moisture_classes_stag',
                'rnn_state_vars', 'rnn_state_units')
        )

        ncrnn_state[:, :, :, :, :] = rnn_states 
        #for v in data_vars:
        #    d.createVariable(v, 'f4', ('south_north', 'west_east'))[:,:]=data_vars[v]

        ds.close()

    def to_geogrid(self, preds, path, index, lats=[], lons=[]):
        """
        Store model to geogrid files
        """
        # TODO: check (xy) vs (yx) order. I think it's flipped in to_geogrid in existing fmda
        not_coord = len(lats) == 0 or len(lons) == 0
        test_latslons=True
        if not_coord:
            test_latslons=False

        logging.info("fmda.rnn_moisture_model.to_geogrid path=%s lats %s lons %s" % (path, inq(lats), inq(lons)))
        logging.info("fmda.rnn_moisture_model.to_geogrid: geogrid_index="+str(index))
        ensure_dir(path)
        
        xsize, ysize, n = preds.shape

        if not not_coord:
            x=int(xsize*0.5)
            y=int(ysize*0.5)
            index.update({'known_x':float(y),'known_y':float(x),'known_lat':float(lats[x-1,y-1]),'known_lon':float(lons[x-1,y-1])})
            logging.info("fmda.rnn_moisture_model.to_geogrid: geogrid updated="+str(index))

        FMC_GC = np.zeros((xsize, ysize, n+2))
        FMC_GC[:,:,:] = preds
        if test_latslons:
            logging.info("fmda.rnn_moisture_model.to_geogrid: storing lons lats to FMC_GC(:,:,-2:) to test in WRF against XLONG and XLAT")
            FMC_GC[:,:,-2] = lons
            FMC_GC[:,:,-1] = lats

        # Write recurrent states, flatten to (ny, nx, k*n_rnn_vars*n_units)
        RNN_STATES = self.states_to_grid((xsize, ysize))
        ny, nx, k, n_rnn_vars, n_units = RNN_STATES.shape
        assert (ny, nx) == (xsize, ysize)
        RNN_STATES = RNN_STATES.reshape(xsize, ysize, k*n_rnn_vars*n_units)


        write_geogrid_var(path,'FMC_GC',FMC_GC,index,bits=32)
        write_geogrid_var(path,'RNN_STATES',RNN_STATES,index,bits=32)


    def _make_rnn_state_labels():

        """
        Given internal recurrent state, create list of field names for use within params in to_wps_format

        Returns: list with naming convention, (N units)

        [FM1_H1, ..., FM1_HN, FM1_C1, ..., FM1_CN,
         FM10_H1, ..., FM10_HN, FM10_C1, ..., FM10_CN,
         FM100_H1, ..., FM100_HN, FM100_C1, ..., FM100_CN,
         FM1000_H1, ..., FM1000_HN, FM1000_C1, ..., FM1000_CN]

        TODO: decide on count from 0 or from 1, fortran vs python
        """

        states = self.get_states()
        
        # Get number of units, double check it makes sense
        n_units = states.size // (len(self.FUEL_CLASSES) * 2)
        lstm_idx = self.params["hidden_layers"].index("lstm")
        lstm_units = self.params["hidden_units"][lstm_idx]
        if n_units != lstm_units:
            raise ValueError(f"State units ({n_units}) do not match LSTM units ({lstm_units}).")
        
        labels = [
            f"{fuel_class.upper()}_{state}{unit}"
            for fuel_class in self.FUEL_CLASSES
            for state in ("H", "C")
            for unit in range(1, n_units + 1)
        ]
        return labels        

    def _make_rnn_state_descriptions():

        """
        Given internal recurrent state, create list of descriptive names for use
        within params in to_wps_format.

        Returns descriptions ordered by fuel class, state type, then unit.

        """

        states = self.get_states()
        n_units = states.size // (len(self.FUEL_CLASSES) * 2)

        lstm_idx = self.params["hidden_layers"].index("lstm")
        lstm_units = self.params["hidden_units"][lstm_idx]
        if n_units != lstm_units:
            raise ValueError(
                f"State units ({n_units}) do not match LSTM units ({lstm_units})."
            )

        return [
            f"{fuel_class.upper()} Fuel Moisture {state_name} State, Unit {unit}"
            for fuel_class in self.FUEL_CLASSES
            for state_name in ("Hidden", "Cell")
            for unit in range(1, n_units + 1)
        ]

    def to_wps_format(self, preds, path, index, lats, lons, time_tag):
        """
        Store model to wps format files
        """
        test_latslons=True

        logging.info("fmda.rnn_moisture_model.to_wps_format path=%s lats %s lons %s" % (path, inq(lats), inq(lons)))
        logging.info("fmda.rnn_moisture_model.to_wps_format: index="+str(index))
        ny, nx, n = preds.shape

        rnn_states = states_to_grid(grid_shape=(ny, nx)) 

        m = n + 2 + len(states)
        var = np.zeros((ny, nx, m))
        var[:,:,:n] = preds
        var[:,:,n] = lons
        var[:,:,n+1] = lats
        #var[:,:,-len(states)] = states
        #var[:,:,n:] = self.m_ext[:,:,n-2:]
        arrs = [var[:,:,k].swapaxes(0,1) for k in range(m)]
        startloc = "SWCORNER"
        startlat = lats[0,0]
        startlon = lons[0,0]
        fm_fields = ["{}H".format(10**i) for i in range(n-2)]
        fm_descrs = ["{}h Fuel Moisture Content".format(10**i) for i in range(n-2)]
        params = {
            'ifv': 5, 'hdate': "{}:00:00".format(time_tag), 'xfcst': 0.,
            'map_source': "WRF-SFIRE Wildland Fire Information and Forecasting System",
            'field': fm_fields + ["FMXLON", "FMXLAT", "FMEP0", "FMEP1"],
            'units': ["1" for _ in range(n-2)] + ["degrees", "degrees", "1", "1"],
            'desc': fm_descrs + ["Longitude for testing", "Latitude for testing",
                    "Drying/Wetting Equilibrium Adjustment", "Rain Equilibrium Adjustment"],
            'xlvl': 200100., 'nx': nx, 'ny': ny, 'iproj': 3, 'startloc': startloc,
            'startlat': startlat, 'startlon': startlon, 'dx': index['dx'], 'dy': index['dy'],
            'xlonc': index['stdlon'], 'truelat1': index['truelat1'], 'truelat2': index['truelat2'],
            'earth_radius': index['radius'], 'is_wind_earth_rel': 0, 'slab': arrs
        }
        WPSFormat.from_params(**params).to_file(osp.join(path,'FMDA:{}'.format(time_tag)))


    @classmethod
    def from_netcdf(cls, netcdf_path, weights_path):
        """
        Construct RNN from stored netCDF file, which has previous recurrent state, params for model architecture, and weights path. Weights stored externally are not saved in netCDF.

        :param netcdf_path: the path to the netCDF4 file
        :param weights_path: the path to the h5 file of RNN weights that should match netCDF params architecture
        """
        import netCDF4

        # Read and simple checks
        logging.info("reading from netCDF file " + netcdf_path+" and model weights "+weights_path)
        ds = netCDF4.Dataset(netcdf_path)
        params = json.loads(ds.params)
        warps = json.loads(ds.warps)
        FUEL_CLASSES = json.loads(ds.FUEL_CLASSES)
        stored_weights_path = ds.weights_path
        if stored_weights_path != weights_path:
            logging.warning(
                f"Supplied weights path differs from the path stored in the netCDF: "
                f"{weights_path} != {stored_weights_path}"
            )

        # Extract needed saved info
        rnn_states_grid = ds.variables["RNN_STATES"][:]
        d0, d1, k, n_rnn_vars, n_rnn_units = rnn_states_grid.shape
        rnn_states = {}
        for i, fuel_class in enumerate(FUEL_CLASSES):
            states_i = rnn_states_grid[:, :, i, :, :].reshape(
                (d0 * d1, n_rnn_vars, n_rnn_units)
            )
            states_i_list = [states_i[:, j, :] for j in range(n_rnn_vars)]
            rnn_states[fuel_class] = states_i_list


        #Construct RNN objects with stored states
        model = cls(params, weights_path, warps)
        model.set_states(rnn_states)
        
        return model
        


# Helper functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def convert_feature_names(features, name_map=RNN_FEATURE_TO_WRFXPY):
    """Convert feature names using name_map; leave unmapped names unchanged."""
    return [name_map.get(feature, feature) for feature in features]


def to_netcdf(path, arr, valid_times, fuel_vars=("FM1", "FM10")):
    """
    Standalone version of the model.to_netcdf from the FMDA object class. 
    The RNN class doesn't lend itself as cleanly to saving predictions as the model object,
    so making a standalone
    """
    import netCDF4
    d = netCDF4.Dataset(path, "w", format="NETCDF4")

    arr = np.asarray(arr)
    ny, nx, nt, nfuel = arr.shape
    if len(fuel_vars) != nfuel: raise ValueError( f"len(fuel_vars)={len(fuel_vars)} must match last dim nfuel={nfuel}" )


    d.createDimension("south_north", ny)
    d.createDimension("west_east", nx)
    d.createDimension("time", nt)

    time = d.createVariable("time", "f8", ("time",))
    time.units = "hours since 1970-01-01 00:00:00 UTC"
    time.calendar = "standard"
    time[:] = netCDF4.date2num(
        valid_times.astype("datetime64[ms]").astype(object),
        units=time.units,
        calendar=time.calendar,
    )
    for i, var_name in enumerate(fuel_vars):
        v = d.createVariable( var_name, "f4", ("south_north", "west_east", "time"), zlib=True, complevel=4, )
        v.long_name = f"{var_name} fuel moisture prediction"
        v.coordinates = "time"
        v[:] = arr[..., i]

    d.Conventions = "CF-1.8"

    d.close()
