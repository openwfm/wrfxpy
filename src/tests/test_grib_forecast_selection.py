from datetime import datetime, timezone
import pickle
import re
import sys
from types import SimpleNamespace

import pytest

# utils.py imports dill, which is absent from the portable test environment.
# Pickle provides the dump/load calls used by utils.py if another test needs them.
sys.modules.setdefault('dill', pickle)

import ingest.grib_forecast as grib_forecast_module
from ingest.HRRR import HRRR
from ingest.grib_source import GribError
import utils


def hrrr_without_local_files(tmp_path, monkeypatch, status_for):
    hrrr = HRRR({
        'ingest_path': str(tmp_path / 'ingest'),
        'cache_path': str(tmp_path / 'cache'),
        'sys_install_path': str(tmp_path),
    })
    hrrr.remote_url = ['s3://archive.example/']
    hrrr.browse_aws = 'https://archive.example/index.html'
    monkeypatch.setattr(hrrr, 'colmet_missing', lambda prefix, files: files)
    monkeypatch.setattr(hrrr, 'grib_available_locally', lambda path: False)
    monkeypatch.setattr(
        grib_forecast_module,
        'readhead',
        lambda url, msg_level=0: SimpleNamespace(status_code=status_for(url)),
    )
    return hrrr


def test_hrrr_selects_one_extended_cycle_after_short_cycles(tmp_path, monkeypatch):
    checked = []

    def status_for(url):
        cycle_hour, forecast_hour = map(
            int, re.search(r't(\d\d)z\.wrfprsf(\d\d)', url).groups())
        checked.append((cycle_hour, forecast_hour))
        return 200 if cycle_hour == 6 or forecast_hour <= 18 else 404

    hrrr = hrrr_without_local_files(tmp_path, monkeypatch, status_for)
    downloaded = []
    monkeypatch.setattr(
        hrrr, 'download_grib_many',
        lambda url, paths, workers: downloaded.extend(paths))

    manifest = hrrr.retrieve_gribs(
        datetime(2026, 6, 27, 8, tzinfo=timezone.utc),
        datetime(2026, 6, 28, 10, tzinfo=timezone.utc),
        ref_utc=datetime(2026, 7, 30, tzinfo=timezone.utc),
    )

    assert checked[0] == (6, 28)
    assert all(cycle_hour == 6 for cycle_hour, _ in checked)
    assert downloaded[0].endswith('hrrr.t06z.wrfprsf02.grib2')
    assert downloaded[-1].endswith('hrrr.t06z.wrfprsf28.grib2')
    assert all('hrrr.t06z.' in path for path in manifest.grib_files)


def test_hrrr_rejects_explicit_f00_cycle(tmp_path, monkeypatch):
    hrrr = hrrr_without_local_files(tmp_path, monkeypatch, lambda url: 200)
    with pytest.raises(GribError, match='starts before f01'):
        hrrr.retrieve_gribs(
            datetime(2026, 6, 27, 8, tzinfo=timezone.utc),
            datetime(2026, 6, 27, 12, tzinfo=timezone.utc),
            cycle_start=datetime(2026, 6, 27, 8, tzinfo=timezone.utc),
        )


def test_access_error_does_not_select_an_older_cycle(tmp_path, monkeypatch):
    checked = []

    def inaccessible(url):
        checked.append(url)
        return -1

    hrrr = hrrr_without_local_files(tmp_path, monkeypatch, inaccessible)
    with pytest.raises(GribError, match='availability check returned -1'):
        hrrr.retrieve_gribs(
            datetime(2026, 6, 27, 8, tzinfo=timezone.utc),
            datetime(2026, 6, 27, 12, tzinfo=timezone.utc),
        )
    assert len(checked) == 1


def test_normal_cycle_search_stops_after_three_attempts(tmp_path, monkeypatch):
    checked = []

    def missing(url):
        checked.append(int(re.search(r't(\d\d)z', url).group(1)))
        return 404

    hrrr = hrrr_without_local_files(tmp_path, monkeypatch, missing)
    hrrr.cycle_search_attempts = 3
    with pytest.raises(GribError, match='no complete HRRR cycle'):
        hrrr.retrieve_gribs(
            datetime(2026, 6, 27, 8, tzinfo=timezone.utc),
            datetime(2026, 6, 27, 12, tzinfo=timezone.utc),
        )
    assert checked == [7, 6, 5]


def test_readhead_returns_404_without_retrying(monkeypatch):
    calls = []

    def not_found(url):
        calls.append(url)
        return SimpleNamespace(status_code=404)

    monkeypatch.setattr(utils.requests, 'head', not_found)
    assert utils.readhead('https://archive.example/missing').status_code == 404
    assert len(calls) == 1
