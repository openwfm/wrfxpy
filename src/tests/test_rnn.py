import sys
sys.path.insert(0, "src")
import os.path as osp
from ml_fmda.utils import read_yml, Dict
from ml_fmda.moisture_rnn import TimeWarpedFuelClassPredictors
import json
from fmda.fuel_moisture_rnn import RNNMoistureModel
import pandas as pd
import joblib
import numpy as np

test_path = "wksp/test_rnn_output"

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
rnn0 = TimeWarpedFuelClassPredictors(
    params=params,
    weights_path=fm10_weights_path,
    warps=warps
)


rnn = RNNMoistureModel(
    params=params,
    weights_path=fm10_weights_path,
    warps=warps
) 

# Get states, should be None (TODO: switch to zeros?)
rnn0.get_states()
rnn.get_states()

# Predict, then get states, should still be None since not set with vanilla predict
ny = nx = 10
X = np.zeros((ny*nx, 24, len(params.features_list)))
p0 = rnn0.predict(X)
p1 = rnn.predict(X)
assert np.allclose(p0, p1)



# Predict cycle, then get states, should be valid list of states
p0 = rnn0.predict_cycle(X)
p1 = rnn.predict_cycle(X)
assert np.allclose(p0, p1)

states0 = rnn0.get_states()
states1 = rnn.get_states()

preds_grid = p1.reshape(ny, nx, p1.shape[1], p1.shape[2])
rnn_states_grid = rnn.states_to_grid((ny, nx))

# Test writing out with one time slice, then reading
pred_slice = preds_grid[:,:,0,:].squeeze()
rnn.to_netcdf(
    path=f"{test_path}.nc",
    preds=pred_slice,
    grid_shape=(ny, nx),
    data_vars=None
)

rnn2 = RNNMoistureModel.from_netcdf(f"{test_path}.nc", fm10_weights_path) 

states1 = rnn.get_states()
states2 = rnn2.get_states()

for fuel_class in rnn.FUEL_CLASSES:
    for state1, state2 in zip(states1[fuel_class], states2[fuel_class]):
        np.testing.assert_array_equal(state1, state2)

assert states1.keys() == states2.keys()

p1 = rnn.predict_cycle(X, reset_states=False)
p2 = rnn2.predict_cycle(X, reset_states=False)

np.testing.assert_allclose(p1, p2, rtol=1e-6, atol=1e-6)

breakpoint()
