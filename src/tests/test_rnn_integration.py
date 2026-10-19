"""Check RNN output formats and continuation using the configured model assets."""

import sys

sys.path.insert(0, "src")

import netCDF4
import numpy as np
import pytest

from fmda.fuel_moisture_model import BlendMoistureModel
from geo.write_geogrid import read_geogrid
from ml_fmda.moisture_rnn import predict_auto_batch
from wrf.wps_format import WPSFormat

# Reuse the supplied test module's model setup and weight selection.
from test_rnn import index, params, rnn


def rectangular_case():
    """Return varying inputs and coordinates for a 3 by 5 grid."""
    ny, nx = 3, 5
    rows, cols = np.indices((ny, nx))
    lats = 38.5 - rows * 0.03
    lons = 262.5 + cols * 0.04
    inputs = np.random.default_rng(20261019).normal(
        scale=0.2, size=(ny * nx, 6, len(params.features_list))
    )
    return inputs, lats, lons


def assert_wps_fields(path, expected):
    """Check field names, grid dimensions, and values after WPS readback."""
    records = WPSFormat.from_file(str(path))
    assert len(records) == len(expected)
    assert [record["field"].decode() for record in records] == list(expected)
    for record, values in zip(records, expected.values()):
        assert (record["ny"], record["nx"]) == values.shape
        # WPS slabs store x before y; the model arrays store y before x.
        np.testing.assert_allclose(record["slab"].T, values, rtol=1e-6, atol=1e-6)


@pytest.mark.parametrize("output_format", ["netcdf", "geogrid", "wps"])
def test_rectangular_rnn_outputs(tmp_path, output_format):
    """Preserve cell, fuel, hidden-state, and cell-state ordering in exports."""
    inputs, lats, lons = rectangular_case()
    predictions = rnn.predict_cycle(inputs, reset_states=True, verbose=0)
    # Export the last prediction, matching the saved end-of-sequence states.
    moisture = predictions[:, -1, :].reshape(3, 5, 4) / 100
    states = rnn.states_to_grid((3, 5))
    assert np.isfinite(moisture).all()
    assert np.isfinite(states).all()

    # Run formats independently so a failure in one does not hide another.
    if output_format == "netcdf":
        nc_path = tmp_path / "rnn.nc"
        rnn.to_netcdf(str(nc_path), moisture, (3, 5), data_vars=None)
        with netCDF4.Dataset(nc_path) as dataset:
            np.testing.assert_allclose(dataset["FMC_GC"][:], moisture, rtol=1e-6)
            np.testing.assert_array_equal(dataset["RNN_STATES"][:], states)
        return

    if output_format == "geogrid":
        geo_path = tmp_path / "rnn.geo"
        rnn.to_geogrid(moisture, str(geo_path), index.copy(), lats, lons)
        fields, metadata = read_geogrid(str(geo_path / "FMC_GC"))
        expected = np.concatenate((moisture, lons[..., None], lats[..., None]), axis=-1)
        np.testing.assert_allclose(
            fields, expected, rtol=1e-6, atol=metadata["scale_factor"] / 2
        )
        fields, metadata = read_geogrid(str(geo_path / "RNN_STATES"))
        np.testing.assert_allclose(
            fields, states.reshape(3, 5, -1), rtol=1e-6,
            atol=metadata["scale_factor"] / 2,
        )
        return

    time_tag = "2026-10-19_00"
    rnn.to_wps_format(moisture, str(tmp_path), index.copy(), lats, lons, time_tag)
    expected = {}
    for fuel_index, hours in enumerate((1, 10, 100, 1000)):
        expected[f"{hours}H"] = moisture[..., fuel_index]
    expected.update(FMXLON=lons, FMXLAT=lats)
    for fuel_index, hours in enumerate((1, 10, 100, 1000)):
        for state_index, state_name in enumerate(("H", "C")):
            for unit in range(states.shape[-1]):
                expected[f"FM{hours}_{state_name}{unit + 1}"] = states[
                    :, :, fuel_index, state_index, unit
                ]
    assert_wps_fields(tmp_path / f"RNN:{time_tag}", expected)


@pytest.mark.parametrize("batch_size", [4, 64])
def test_split_sequence_and_batching(batch_size):
    """Match a full sequence with continuation across two prediction calls."""
    inputs, _, _ = rectangular_case()
    full = rnn.predict_cycle(inputs, reset_states=True, batch_size=64, verbose=0)
    expected_states = rnn.states_to_grid((3, 5)).copy()
    # Fifteen cells and batches of four exercise an incomplete final batch.
    first = predict_auto_batch(
        rnn, inputs[:, :3], batch_sizes=(batch_size,), reset_states=True, verbose=0
    )
    second = predict_auto_batch(
        rnn, inputs[:, 3:], batch_sizes=(batch_size,), reset_states=False, verbose=0
    )
    np.testing.assert_allclose(
        np.concatenate((first, second), axis=1), full, rtol=1e-6, atol=1e-6
    )
    np.testing.assert_allclose(
        rnn.states_to_grid((3, 5)), expected_states, rtol=1e-6, atol=1e-6
    )


def test_blend_outputs(tmp_path):
    """Preserve the arithmetic blend and correction fields in each format."""
    _, lats, lons = rectangular_case()
    # Distinct cell and field values expose transposition and field swaps.
    timelag = np.arange(3 * 5 * 6, dtype=float).reshape(3, 5, 6) / 1000
    moisture = np.arange(3 * 5 * 4, dtype=float).reshape(3, 5, 4) / 200
    original = timelag.copy()
    blend = BlendMoistureModel(timelag, moisture)
    expected = timelag.copy()
    expected[..., :4] = (timelag[..., :4] + moisture) / 2
    np.testing.assert_array_equal(timelag, original)
    np.testing.assert_array_equal(blend.m_ext, expected)

    nc_path = tmp_path / "blend.nc"
    blend.to_netcdf(str(nc_path))
    with netCDF4.Dataset(nc_path) as dataset:
        np.testing.assert_allclose(dataset["FMC_GC"][:], expected, rtol=1e-6)

    geo_path = tmp_path / "blend.geo"
    blend.to_geogrid(str(geo_path), index.copy(), lats, lons)
    fields, metadata = read_geogrid(str(geo_path / "FMC_GC"))
    geo_expected = np.concatenate(
        (expected[..., :4], lons[..., None], lats[..., None]), axis=-1
    )
    np.testing.assert_allclose(
        fields, geo_expected, rtol=1e-6, atol=metadata["scale_factor"] / 2
    )
    corrections, metadata = read_geogrid(str(geo_path / "FMEP"))
    np.testing.assert_allclose(
        corrections, expected[..., 4:], rtol=1e-6,
        atol=metadata["scale_factor"] / 2,
    )

    time_tag = "2026-10-19_00"
    blend.to_wps_format(str(tmp_path), index.copy(), lats, lons, time_tag)
    wps_expected = {}
    for fuel_index, hours in enumerate((1, 10, 100, 1000)):
        wps_expected[f"{hours}H"] = expected[..., fuel_index]
    wps_expected.update(
        FMXLON=lons, FMXLAT=lats, FMEP0=expected[..., 4], FMEP1=expected[..., 5]
    )
    assert_wps_fields(tmp_path / f"BLEND:{time_tag}", wps_expected)
