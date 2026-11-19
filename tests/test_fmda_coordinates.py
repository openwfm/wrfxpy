"""Run from wrfxpy: python -B -m unittest discover -s tests -p test_fmda_coordinates.py."""
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
# Existing ingest imports read site configuration; these tests need none of it.
with patch('utils.load_sys_cfg', return_value={}):
    from fmda.fuel_moisture_da import check_fmda_coordinates


class CoordinateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def make_input(self, domain=2, offset=0.0, classes=6, missing=False):
        path = self.root / ('wrfinput_d%02d' % domain)
        with netCDF4.Dataset(path, 'w') as ds:
            for name, size in [('Time', 1), ('fuel_moisture_classes_stag', classes),
                               ('south_north', 2), ('west_east', 3)]:
                ds.createDimension(name, size)
            dims = ('Time', 'south_north', 'west_east')
            ds.createVariable('XLONG', 'f8', dims)[:] = -110.0
            ds.createVariable('XLAT', 'f8', dims)[:] = 40.0
            if not missing:
                fmc = ds.createVariable('FMC_GC', 'f8', (dims[0], 'fuel_moisture_classes_stag') + dims[1:])
                fmc[:] = 0.1
                fmc[0, -2] = -110.0
                fmc[0, -1] = 40.0 + offset
        return path

    def test_all_fire_domains_and_read_only(self):
        paths = [self.make_input(2), self.make_input(3, classes=5)]
        before = [hashlib.sha256(p.read_bytes()).digest() for p in paths]
        with self.assertLogs(level='INFO') as logs:
            check_fmda_coordinates(self.root, [0, 1, 1])
        self.assertEqual(len(logs.output), 2)
        self.assertIn('RMS=0.000000 m max=0.000000 m', logs.output[0])
        self.assertEqual(before, [hashlib.sha256(p.read_bytes()).digest() for p in paths])

    def test_submeter_error_and_scalar_ifire(self):
        self.make_input(1, offset=0.5 / 6371000 * 180 / np.pi)
        with self.assertLogs(level='INFO') as logs:
            check_fmda_coordinates(self.root, 1)
        self.assertIn('RMS=0.500000 m max=0.500000 m', logs.output[0])

    def test_displaced_coordinates_fail(self):
        self.make_input(offset=0.002)
        with self.assertRaisesRegex(ValueError, 'require max < 1 m'):
            check_fmda_coordinates(self.root, [0, 1])

    def test_maximum_not_just_rms_is_checked(self):
        path = self.make_input()
        with netCDF4.Dataset(path, 'r+') as ds:
            ds['FMC_GC'][0, -1, 0, 0] += 1.5 / 6371000 * 180 / np.pi
        with self.assertRaisesRegex(ValueError, 'max=1.500000 m'):
            check_fmda_coordinates(self.root, [0, 1])

    def test_equivalent_longitude_wrap(self):
        path = self.make_input()
        with netCDF4.Dataset(path, 'r+') as ds:
            ds['FMC_GC'][0, -2] = 250.0
        check_fmda_coordinates(self.root, [0, 1])

    def test_missing_or_invalid_diagnostics_fail(self):
        for value in [np.nan, np.inf, netCDF4.default_fillvals['f8'], 91.0]:
            with self.subTest(value=value):
                path = self.make_input()
                with netCDF4.Dataset(path, 'r+') as ds:
                    ds['FMC_GC'][0, -1, 0, 0] = value
                with self.assertRaises(ValueError):
                    check_fmda_coordinates(self.root, [0, 1])
        self.make_input(missing=True)
        with self.assertRaisesRegex(ValueError, 'missing or incorrectly dimensioned FMC_GC'):
            check_fmda_coordinates(self.root, [0, 1])
        self.make_input(classes=2)
        with self.assertRaisesRegex(ValueError, 'moisture plus lon/lat'):
            check_fmda_coordinates(self.root, [0, 1])

    def test_no_fire_or_missing_file_fails(self):
        with self.assertRaisesRegex(ValueError, 'no fire domains'):
            check_fmda_coordinates(self.root, [0, 0])
        with self.assertRaises(OSError):
            check_fmda_coordinates(self.root, [0, 1])


if __name__ == '__main__':
    unittest.main()
