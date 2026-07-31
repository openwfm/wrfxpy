from pathlib import Path

from ingest.grib_source import GribSource
from ingest import downloader


def test_zero_length_grib_is_not_available(tmp_path):
    grib_path = tmp_path / 'forecast.grib2'
    grib_path.touch()
    Path(str(grib_path) + '.size').write_text('0')

    source = object.__new__(GribSource)
    assert not source.grib_available_locally(str(grib_path))


def test_failed_availability_check_raises(tmp_path, monkeypatch):
    def unavailable(*args, **kwargs):
        raise OSError('unavailable')

    monkeypatch.setattr(downloader, 'request_url', unavailable)
    monkeypatch.setattr(downloader.time, 'sleep', lambda seconds: None)
    monkeypatch.setattr(downloader.random, 'random', lambda: 0)

    try:
        downloader.download_url('https://example.invalid/file',
                                str(tmp_path / 'forecast.grib2'),
                                max_retries=0)
    except downloader.DownloadError:
        pass
    else:
        raise AssertionError('missing remote file did not raise DownloadError')
