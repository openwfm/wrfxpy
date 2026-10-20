import sys
sys.path.insert(0, "src")
import os.path as osp
import os
from ml_fmda.utils import read_yml, Dict
from ml_fmda.moisture_rnn import TimeWarpedFuelClassPredictors
import json
from fmda.fuel_moisture_rnn import RNNMoistureModel, read_states_from_netcdf
from geo.write_geogrid import read_geogrid
import pandas as pd
import joblib
import numpy as np
import warnings

test_path = "wksp/test_rnn_output"
os.makedirs(test_path, exist_ok=True)
rnn_cfg = Dict(json.load(open("etc/rnn_cycler.json")))

params = Dict(read_yml(osp.join(rnn_cfg.rnn_model_path, "params.yaml")))
median_seed = pd.read_csv(osp.join(rnn_cfg.rnn_model_path, "median_seed.csv")).seed[0]
scaler = joblib.load(osp.join(rnn_cfg.rnn_model_path, f"seed_{median_seed}", "scaler.joblib"))
fm10_weights_path = osp.join(rnn_cfg.rnn_model_path, f"seed_{median_seed}", "rnn.weights.h5")
# Transfer learning, get twarps
twarp_summary = pd.read_csv(osp.join(rnn_cfg.transfer_path, "all_seeds_summary.csv"))
row = twarp_summary[twarp_summary.seed == median_seed].iloc[0]

warps = {
    "fm1": (row["bi_1"], row["bf_1"]),
    "fm100": (row["bi_100"], row["bf_100"]),
    "fm1000": (row["bi_1000"], row["bf_1000"]),
}

# Building a model will raise a warning if parameters related to training are provided,
# since the model class does not support training. The `timesteps` parameter raises a 
# warning, since it fixes sequence length during training.
# We turn this off and raise an informative warning.
# In practice, we copy the params dictionary from the training m
# model to create the prediction model, so the warning gets raised. This is intended, so in this test we capture the warning
with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="timesteps should be None for flexible prediction",
        category=UserWarning,
    )
    rnn0 = TimeWarpedFuelClassPredictors(
        params=params,
        weights_path=fm10_weights_path,
        warps=warps
    )


with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="timesteps should be None for flexible prediction",
        category=UserWarning,
    )
    rnn = RNNMoistureModel(
        params=params,
        weights_path=fm10_weights_path,
        warps=warps
    ) 

# Predict, then get states, should still be None since not set with vanilla predict
ny = nx = 10
index = {
        "projection": "lambert",
        "dx" : 3000.0,
        "dy" : -3000.0,
        "truelat1" : 38.5,
        "truelat2" : 38.5,
        "stdlon" : 262.5,
        "radius" : 6370000.0
    }
lats = np.linspace(38.5, 38.5 - 0.27 * (ny - 1), ny)[:, None] * np.ones((1, nx))
lons = np.ones((ny, 1)) * np.linspace(262.5, 262.5 + 0.35 * (nx - 1), nx)[None, :]

X = np.zeros((ny*nx, 24, len(params.features_list)))
def test_initial_predictions_match():

    p0 = rnn0.predict(X)
    p1 = rnn.predict(X)

    assert np.allclose(p0, p1, equal_nan=False)

def test_predict_cycle_matches():
    p0 = rnn0.predict_cycle(X, reset_states=True)
    p1 = rnn.predict_cycle(X, reset_states=True)
    assert np.allclose(p0, p1, equal_nan=False)

def test_predict_cycle_states_match():
    rnn0.predict_cycle(X, reset_states=True)
    rnn.predict_cycle(X, reset_states=True)

    states0 = rnn0.get_states()
    states1 = rnn.get_states()

    assert states0.keys() == states1.keys()

    for fuel_class in states0:
        assert len(states0[fuel_class]) == len(states1[fuel_class])
        for s0, s1 in zip(states0[fuel_class], states1[fuel_class]):
            assert np.allclose(s0, s1, equal_nan=False)

def test_netcdf_roundtrip(tmp_path):
    preds = rnn.predict_cycle(X, reset_states=True)

    preds_grid = preds.reshape(ny, nx, preds.shape[1], preds.shape[2])
    pred_slice = preds_grid[:, :, 0, :].squeeze()

    path = tmp_path / "test_rnn_output.nc"
    rnn.to_netcdf(
        path=str(path),
        preds=pred_slice,
        grid_shape=(ny, nx),
        data_vars=None,
    )

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="timesteps should be None for flexible prediction",
            category=UserWarning,
        )
        rnn2 = RNNMoistureModel.from_netcdf(str(path), fm10_weights_path)

    states1 = rnn.get_states()
    states2 = rnn2.get_states()

    assert states1.keys() == states2.keys()
    for fuel_class in rnn.FUEL_CLASSES:
        for state1, state2 in zip(states1[fuel_class], states2[fuel_class]):
            np.testing.assert_array_equal(state1, state2)

    p1 = rnn.predict_cycle(X, reset_states=False)
    p2 = rnn2.predict_cycle(X, reset_states=False)

    np.testing.assert_allclose(p1, p2, rtol=1e-6, atol=1e-6)

def test_netcdf_state_restoration(tmp_path):
    p1 = rnn.predict_cycle(X, reset_states=True)

    preds_grid = p1.reshape(ny, nx, p1.shape[1], p1.shape[2])
    pred_slice = preds_grid[:, :, 0, :].squeeze()

    path = tmp_path / "test_rnn_output.nc"
    rnn.to_netcdf(
        path=str(path),
        preds=pred_slice,
        grid_shape=(ny, nx),
        data_vars=None,
    )

    X2 = X + 0.5
    p2 = rnn.predict_cycle(X2, reset_states=False)

    rnn_states = read_states_from_netcdf(path=str(path))
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="timesteps should be None for flexible prediction",
            category=UserWarning,
        )
        rnn2 = RNNMoistureModel(
            params=params,
            weights_path=fm10_weights_path,
            warps=warps,
        )

    rnn2.set_states(rnn_states)

    p3 = rnn2.predict_cycle(X2, reset_states=False)

    np.testing.assert_allclose(p2, p3, rtol=1e-6, atol=1e-6)


def test_to_geogrid_roundtrip(tmp_path):
    preds = rnn.predict_cycle(X, reset_states=True)

    preds_grid = preds.reshape(ny, nx, preds.shape[1], preds.shape[2])
    pred_slice = preds_grid[:, :, 0, :]

    path = tmp_path / "test_rnn_output.geo"
    rnn.to_geogrid(
        preds=pred_slice,
        path=str(path),
        index=index,
        lats=lats,
        lons=lons,
    )

    geo_fmc = read_geogrid(osp.join(str(path), "FMC_GC"))
    geo_rnn = read_geogrid(osp.join(str(path), "RNN_STATES"))

    np.testing.assert_allclose(geo_fmc[0][..., :4], pred_slice)

    rnn_states = rnn.states_to_grid((ny, nx))
    _, _, k, n_rnn_vars, nunits = rnn_states.shape
    rnn_states = rnn_states.reshape(ny, nx, k * n_rnn_vars * nunits)

    np.testing.assert_allclose(
        geo_rnn[0],
        rnn_states,
        rtol=1e-6,
        atol=geo_rnn[1]["scale_factor"] / 2,
    )


# Test WPS format

def test_to_wps_format(tmp_path):
    preds = rnn.predict_cycle(X, reset_states=True)

    preds_grid = preds.reshape(ny, nx, preds.shape[1], preds.shape[2])
    pred_slice = preds_grid[:, :, 0, :]

    rnn.to_wps_format(
        preds=pred_slice,
        path=str(tmp_path),
        index=index,
        lats=lats,
        lons=lons,
        time_tag="00:00",
    )

    # TODO: Verify the generated WPS output by reading it back and checking
    # variable names, dimensions, and values against the expected predictions
    # and recurrent states.

