#!/usr/bin/env python
"""Run a fire's WRF-driven and GRIB-driven ForeFire ensembles to a common end time.

Both ensembles have to stop at the same wall-clock instant or their members cannot
be compared -- areas are read at whatever time each one happens to reach, and a
couple of hours of extra burn swamps the difference being measured.  The common end
here is the earlier of the two data ends, the last wrfout and the last GRIB valid
time, and *both* timing tables are trimmed to it.  Red Bank's first pass trimmed
only one side and ended two hours apart, which cost a full re-run.

The WRF side is re-run with overwrite=True rather than reused, because any run made
before the 2026-09-21 clock fix carries wrong `valid_at` stamps, and any run made
before the cron lock may have mostly-empty members.  Geometry in those runs is fine;
labels and completeness are not.

Both take the ForeFire cron lock, so a caller should hold it (a blocking flock, so
this queues rather than skips) and they are sequential by construction.

    PYTHONPATH=src:src/ingest python src/ingest/run_both_ensembles.py <wksp> <nc_dir>
"""

from __future__ import print_function

import argparse
import os
import os.path as osp

import pandas as pd

import forefire as ff


def run_one(tag, wksp, run_dir, tt, dest, ign_latlon, grid_code, cfg, suffix):
    print('\n===== %s: %d rows -> %s' % (tag, len(tt), dest))
    if not osp.isdir(run_dir):
        os.makedirs(run_dir)
    tt2 = ff.make_script_set(run_dir, tt, ign_latlon, grid_code, cfg=cfg, params={})
    tt2.to_csv(osp.join(run_dir, 'timing_table.csv'), index=False)
    ff.run_timing_table(tt2, overwrite=True, cfg=cfg, run_dir=run_dir)
    ff.cleanup_ff_run(run_dir, wksp, grid_code, cfg=cfg, dest_dir=dest, suffix=suffix)
    print('===== %s done' % tag)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('wksp', help='wrfxpy workspace with wrf/ and input.json')
    p.add_argument('nc_dir', help='directory of FF_*.nc built by forefire_grib.py')
    p.add_argument('--skip-wrf', action='store_true',
                   help='only run the GRIB-driven side')
    p.add_argument('--skip-hrrr', action='store_true')
    a = p.parse_args()

    cfg = ff.load_ff_cfg()
    start_utc, ign_utc, ign_latlon, grid_code = ff.read_input(a.wksp)
    print('ignition %s at %s' % (ign_utc, ign_latlon))

    wrf_tt = ff.make_timing_table(a.wksp, ign_utc)
    hrrr_tt = pd.read_csv(osp.join(a.nc_dir, 'timing_table_grib.csv'))
    #make_script_set only dereferences 'wrfout' when a step's netcdf is missing,
    #which it never is on the GRIB side, but the column has to exist
    hrrr_tt['wrfout'] = hrrr_tt['grib']

    common_end = min(wrf_tt['UTC_str'].max(), hrrr_tt['UTC_str'].max())
    print('WRF data ends  %s' % wrf_tt['UTC_str'].max())
    print('HRRR data ends %s' % hrrr_tt['UTC_str'].max())
    print('common end     %s' % common_end)
    wrf_tt = wrf_tt[wrf_tt['UTC_str'] <= common_end].reset_index(drop=True)
    hrrr_tt = hrrr_tt[hrrr_tt['UTC_str'] <= common_end].reset_index(drop=True)
    print('rows: WRF %d, HRRR %d' % (len(wrf_tt), len(hrrr_tt)))

    if not a.skip_wrf:
        run_one('WRF', a.wksp, cfg['run_dir'], wrf_tt, osp.join(a.wksp, 'forefire'),
                ign_latlon, grid_code, cfg, None)
    if not a.skip_hrrr:
        run_one('HRRR', a.wksp, a.nc_dir, hrrr_tt, osp.join(a.wksp, 'forefire_hrrr'),
                ign_latlon, grid_code, cfg, 'hrrr')
    print('\nBOTH ENSEMBLES COMPLETE')


if __name__ == '__main__':
    main()
