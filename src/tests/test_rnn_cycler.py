"""Exercise cycler control flow with synthetic weather and no external services."""

import copy
import importlib
import json
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from test_rnn import params, rnn, rnn_cfg, scaler
from test_rnn_integration import rectangular_case
from utils import Dict


@pytest.fixture
def cycler(tmp_path, monkeypatch):
    """Import the actual cycler using temporary configurations without credentials."""
    # Absolute import paths still work after moving into the temporary case.
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]))
    monkeypatch.chdir(tmp_path)
    etc = tmp_path / "etc"
    etc.mkdir()
    configs = {
        "conf.json": {},
        "fmda_cycler.json": {"workspace_path": "wksp", "regions": {}},
        "rnn_cycler.json": dict(rnn_cfg),
        "tokens.json": {"mesowest": ""},
    }
    for name, value in configs.items():
        (etc / name).write_text(json.dumps(value))
    module = importlib.import_module("hrrr_cycler")
    model_config = Dict(dict(rnn_cfg))
    model_config.params = params
    model_config.scaler = {"mean": scaler.mean_, "scale": scaler.scale_}
    monkeypatch.setattr(module, "rnn_cfg", model_config)
    monkeypatch.setattr(module, "rnn", rnn, raising=False)
    monkeypatch.setattr(module, "sys_cfg", Dict({}))
    return module


def test_feature_order_and_scaling(cycler):
    """Check the operational model's fourteen input columns and moisture units."""
    source_names = [
        "Ed", "Ew", "rh", "t2", "swdown", "ws", "hgt", "lons", "lats",
        "lograin", "hod_sin", "hod_cos", "doy_sin", "doy_cos",
    ]
    assert list(params.features_list) == [
        "Ed", "Ew", "rh", "temp", "solar", "wind", "elev", "lon", "lat",
        "lograin", "hod_sin", "hod_cos", "doy_sin", "doy_cos",
    ]
    data = {}
    expected = np.empty((15, 1, 14))
    for column, name in enumerate(source_names):
        values = np.arange(15).reshape(3, 5) / 100 + column + 1
        data[name] = values.copy()
        # The model was trained with equilibrium moisture in percent.
        factor = 100 if name in ("Ed", "Ew") else 1
        expected[:, 0, column] = (
            values.ravel() * factor - scaler.mean_[column]
        ) / scaler.scale_[column]
    original = copy.deepcopy(data)
    np.testing.assert_allclose(cycler.format_ml_data(data, cycler.rnn_cfg), expected)
    for name in data:
        np.testing.assert_array_equal(data[name], original[name])


def test_synthetic_analysis_forecast_restart(cycler, tmp_path, monkeypatch):
    """Advance two analysis hours and repeat a forecast from saved states."""
    _, lats, lons = rectangular_case()
    data = {"lats": lats, "lons": lons}
    for name, value in {
        "t2": 290, "rh": 30, "hgt": 1500, "swdown": 400, "ws": 3,
        "rain": 0, "snowh": 0, "u10": 2, "v10": 1,
    }.items():
        data[name] = np.ma.array(np.full((3, 5), value, dtype=float), mask=False)
    observations_requested = []

    def no_observations(*args, **kwargs):
        """Record analysis requests and return an explicitly empty observation set."""
        observations_requested.append(args[1])
        return {}

    monkeypatch.setattr(cycler, "load_hrrr_data", lambda *args: copy.deepcopy(data))
    monkeypatch.setattr(cycler, "retrieve_mesowest_observations", no_observations)
    # These tests isolate evolution and file writing; rendering is a separate check.
    monkeypatch.setattr(cycler, "postprocess_cycle", lambda *args: None)
    grib = tmp_path / "synthetic-input-placeholder"
    grib.touch()
    region = Dict({"code": "TEST", "region_id": "Synthetic test", "bbox": [38, -98, 39, -97]})
    workspace = str(tmp_path / "wksp")
    start = datetime(2026, 10, 19)
    steps = [(start, 0), (start + timedelta(hours=1), 0), (start + timedelta(hours=1), 1)]
    for cycle, hour in steps:
        cycler.fmda_advance_region(cycle, region, [str(grib)], workspace, 0, hour, "")
        with netCDF4.Dataset(cycler.compute_model_path(cycle, "TEST", workspace, hour)) as ds:
            timelag = ds["FMC_GC"][:]
        with netCDF4.Dataset(cycler.compute_rnn_path(cycle, "TEST", workspace, hour)) as ds:
            moisture = ds["FMC_GC"][:]
            states = ds["RNN_STATES"][:]
        with netCDF4.Dataset(cycler.compute_blend_path(cycle, "TEST", workspace, hour)) as ds:
            blend = ds["FMC_GC"][:]
        assert np.isfinite(moisture).all()
        np.testing.assert_allclose(blend[..., :4], (timelag[..., :4] + moisture) / 2, rtol=1e-6)
        np.testing.assert_array_equal(blend[..., 4:], timelag[..., 4:])
    assert len(observations_requested) == 2

    # A repeated forecast must reload the analysis states, even after other predictions.
    inputs, _, _ = rectangular_case()
    rnn.predict_cycle(inputs, reset_states=True, verbose=0)
    cycle, hour = steps[-1]
    cycler.fmda_advance_region(cycle, region, [str(grib)], workspace, 0, hour, "")
    with netCDF4.Dataset(cycler.compute_rnn_path(cycle, "TEST", workspace, hour)) as ds:
        np.testing.assert_allclose(ds["FMC_GC"][:], moisture, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(ds["RNN_STATES"][:], states, rtol=1e-6, atol=1e-6)
    assert len(observations_requested) == 2


def test_incomplete_cycle_is_not_complete(cycler, tmp_path):
    """Require missing RNN and blend products to prevent a completed-cycle result."""
    region = Dict({"code": "TEST"})
    cycle = datetime(2026, 10, 19)
    path = Path(cycler.compute_model_path(cycle, "TEST", str(tmp_path)))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    assert not cycler.is_cycle_computed(cycle, region, str(tmp_path))
