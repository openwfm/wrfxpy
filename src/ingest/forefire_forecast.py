#!/usr/bin/env python
"""Run a ForeFire forecast from a wrfxpy job file, with no WRF-SFIRE in the loop.

    PYTHONPATH=src:src/ingest python src/ingest/forefire_forecast.py <job.json> \\
        [--steps 8] [--issue-utc 2026-09-28T12:00:00Z] [--out-dir DIR] [--overwrite]

The whole point of the GRIB path is latency: a WRF-SFIRE run has to finish before
ForeFire can start, and that wait is the blocker.  Everything here comes from an
ignition point, an ignition time, and data that is already on disk.

Three sources feed a run, and all three are now HRRR:

  winds     HRRR **forecast** cycles (00/06/12/18, leads out to 33 h), downscaled by
            WindNinja onto the fire mesh.  `forefire_grib.find_gribs` picks one file
            per valid time -- smallest lead, so the freshest cycle -- and honours
            `issue_utc` so a hindcast cannot use a cycle published after the fire
            started.
  fuels     LANDFIRE, remapped to the Anderson 13 categories.
  moisture  FMDA driven by HRRR (`wksp_fmda/CONUS_HRRR`), read per step and applied
            **per fuel category**: each category mixes the 1-h/10-h/100-h classes with
            its own CAWFE weights, so Md changes at every hourly step rather than
            being one constant for the run.  Forecast hours come from the cycler's
            `-fNN` output; past hours prefer the assimilated analysis.

This is the orchestration that was missing.  `build_step_ncs` built the netcdfs and
`make_script_set`/`run_timing_table` ran ForeFire, but nothing joined them except a
hand-written scratch script per fire.

**Scope note.** WindNinja reads HRRR correctly and *silently garbles* NAM218 -- it
exits 0, prints a plausible valid time, and emits speeds up to 45,632 m/s.  A job whose
`grib_source` is NAM is rejected here rather than run, because the failure does not
announce itself.
"""

from __future__ import absolute_import
from __future__ import print_function

import argparse
import glob
import json
import logging
import os
import os.path as osp
import re
import sys

import pandas as pd

import forefire as ff
import forefire_grib as fg


#The forecast cache cache_grib_files.py fills, as opposed to forefire_grib.HRRR_CACHE,
#which is the *analysis* cache (every file lead 3).  Same product and layout, so only
#the directory differs; see the 2026-09-25 handoff section 8.
HRRR_FORECAST_CACHE = '/data/jhaley/wrfxpy/ingest/HRRR'

#WindNinja's measured ceiling lives in forefire_grib.MAX_CELLS.  A job's domain is
#sized for WRF-SFIRE and can exceed it, so the domain is trimmed rather than the fire
#mesh coarsened: an 8 h fire does not reach the edge of a 30 km box, whereas coarsening
#the solve mesh degrades the ridge-top winds that drive the fastest spread.
DEFAULT_DOMAIN_KM = 30.0


def parse_job_time(value, now=None):
    """A wrfxpy job time as a naive UTC Timestamp.

    Three spellings are in circulation and all three turn up in real files:
    `2026-09-24_21:19:00` in job files and workspace input.json, ISO
    `2026-09-24T21:19:00Z` from anything that has been through pandas, and the
    relative `T+30` / `T-60` (minutes from now) that the job *templates* under tests/
    use.  Accepting only the first is how a template silently becomes a 1970 ignition.
    """
    if isinstance(value, pd.Timestamp):
        return fg.as_naive_utc(value)
    text = str(value).strip()
    m = re.match(r'^T([+-]\d+)$', text)
    if m:
        base = pd.Timestamp.utcnow().tz_localize(None) if now is None else fg.as_naive_utc(now)
        return base + pd.Timedelta(minutes=int(m.group(1)))
    #the job spelling uses '_' between date and time; pandas wants a space or a 'T'
    return fg.as_naive_utc(pd.Timestamp(text.replace('_', 'T')))


def read_job(job_path):
    """(ign_latlon, ign_utc, grid_code, domain_km, fire_dx, grib_source) from a job file.

    Reads the same keys wrfxpy itself reads, so a job that would have driven a
    WRF-SFIRE run drives a ForeFire one unchanged.  The domain is derived from the
    job's own grid -- `cell_size / subgrid_ratio` is the fire mesh WRF-SFIRE would have
    used, and `cell_size * domain_size` its extent -- so the two runs are comparable
    rather than accidentally different.
    """
    with open(job_path) as f:
        job = json.load(f)

    ign = job.get('ignitions', {}).get('1')
    if not ign:
        raise ForeFireJobError('%s has no ignitions["1"]' % job_path)
    ign = ign[0]
    latlon = ign.get('latlon')
    if not latlon or len(latlon) != 2:
        raise ForeFireJobError('%s ignition has no usable latlon' % job_path)
    ign_utc = parse_job_time(ign['time_utc'])

    grid_code = job.get('grid_code') or osp.splitext(osp.basename(job_path))[0]

    dom = (job.get('domains') or {}).get('1', {})
    cell = (dom.get('cell_size') or [None])[0]
    size = (dom.get('domain_size') or [None])[0]
    sub = (dom.get('subgrid_ratio') or [None])[0]
    fire_dx = float(cell) / float(sub) if cell and sub else 30.0
    domain_km = float(cell) * float(size) / 1000.0 if cell and size else DEFAULT_DOMAIN_KM

    return (list(latlon), ign_utc, grid_code, domain_km, fire_dx,
            job.get('grib_source'))


class ForeFireJobError(Exception):
    pass


def fit_domain(domain_km, fire_dx):
    """Shrink domain_km until the fire mesh fits under WindNinja's cell ceiling."""
    n = int(round(domain_km * 1000.0 / fire_dx))
    if n * n <= fg.MAX_CELLS:
        return domain_km
    n_ok = int(fg.MAX_CELLS ** 0.5)
    trimmed = n_ok * fire_dx / 1000.0
    logging.warning('job domain %.0f km at %.0f m is %.2f M cells, past the WindNinja '
                    'ceiling of %.1f M; trimming the domain to %.0f km and keeping the '
                    'fire mesh', domain_km, fire_dx, n * n / 1e6,
                    fg.MAX_CELLS / 1e6, trimmed)
    return trimmed


def run_job(job_path, out_dir=None, steps=8, issue_utc=None, grib_dir=None,
            overwrite=False, cfg=None, domain_km=None, fire_dx=None, params=None):
    """Build the netcdfs, write the .ff scripts, run ForeFire, collect the outputs."""
    if cfg is None:
        cfg = ff.load_ff_cfg()
    latlon, ign_utc, grid_code, job_km, job_dx, grib_source = read_job(job_path)

    domain_km = fit_domain(domain_km or job_km, fire_dx or job_dx)
    fire_dx = fire_dx or job_dx
    grib_dir = grib_dir or cfg.get('hrrr_forecast_cache', HRRR_FORECAST_CACHE)

    #What *does* need checking is the directory, because WindNinja reads HRRR correctly
    #and silently garbles NAM218 -- exit 0, plausible valid time, speeds to 45,632 m/s.
    #find_gribs would simply match nothing in a NAM cache and report "no GRIB valid
    #times", which is safe but tells you the wrong thing.
    if not glob.glob(osp.join(grib_dir, 'hrrr.*')):
        raise ForeFireJobError(
            '%s holds no hrrr.* cycle directories. ForeFire winds must come from HRRR: '
            'WindNinja garbles NAM218 without failing.' % grib_dir)

    run_dir = out_dir or osp.join(cfg['run_dir'], 'job_%s' % grid_code)
    dest_dir = osp.join(run_dir, 'results')
    for d in (run_dir, dest_dir):
        if not osp.isdir(d):
            os.makedirs(d)

    logging.info('%s', grid_code)
    logging.info('  ignition %s at %.5f, %.5f', ign_utc, latlon[0], latlon[1])
    logging.info('  domain %.0f km at %.0f m fire mesh, %d steps', domain_km, fire_dx, steps)
    logging.info('  winds  %s%s', grib_dir,
                 '' if issue_utc is None else ' (as issued %s)' % issue_utc)
    logging.info('  fuel moisture %s', cfg.get('moisture', {}).get('fmda_dir', 'off'))

    #The job's grib_source describes the WRF-SFIRE run it was written for, not this
    #one: ForeFire reads grib_dir and nothing else, so a NAM-sourced job is fine here.
    #Saying so is worth a line, because "grib_source: NAM" in the file people are
    #looking at otherwise implies these winds came from NAM.
    if grib_source and str(grib_source).upper().startswith('NAM'):
        logging.info("  (job grib_source is %r -- that was its WRF-SFIRE run; these "
                     "winds come from the HRRR cache above)", grib_source)


    tt = fg.build_step_ncs(latlon[0], latlon[1], ign_utc, run_dir, grid_code,
                           steps=steps, domain_km=domain_km, fire_dx=fire_dx,
                           grib_dir=grib_dir, issue_utc=issue_utc,
                           overwrite=overwrite)
    if tt is None or not len(tt):
        raise ForeFireJobError('no GRIB valid times cover this ignition')

    #make_script_set dereferences 'wrfout' only when a step's netcdf is missing, which
    #it never is here -- but the column has to exist for the frame to be interchangeable
    #with the one make_timing_table builds.
    tt['wrfout'] = tt['grib']

    #make_script_set applies the per-step, per-fuel-category moisture itself, via
    #md_series -> write_fuel_table, so the fuel table changes at every step.
    tt = ff.make_script_set(run_dir, tt, latlon, grid_code, cfg=cfg,
                            params=dict(params or {}))
    tt.to_csv(osp.join(run_dir, 'timing_table.csv'), index=False)

    ff.run_timing_table(tt, overwrite=overwrite, cfg=cfg, run_dir=run_dir,
                        dest_dir=dest_dir)
    ff.cleanup_ff_run(run_dir, None, grid_code, cfg=cfg, dest_dir=dest_dir)
    logging.info('done -> %s', dest_dir)
    return dest_dir


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('job', help='a wrfxpy job json, or a workspace input.json')
    p.add_argument('--out-dir', default=None,
                   help='run directory; default is <run_dir>/job_<grid_code>')
    p.add_argument('--steps', type=int, default=8,
                   help='hourly steps to forecast (default 8, the low-latency target)')
    p.add_argument('--issue-utc', default=None,
                   help='use only HRRR cycles published by this time, so a hindcast '
                        'cannot see the future. Omit for the best files on disk.')
    p.add_argument('--grib-dir', default=None,
                   help='default: the HRRR forecast cache, %s' % HRRR_FORECAST_CACHE)
    p.add_argument('--domain-km', type=float, default=None,
                   help='override the domain size the job implies')
    p.add_argument('--fire-dx', type=float, default=None,
                   help='override the fire mesh the job implies')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--verbose', action='store_true')
    a = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format='%(levelname)s %(message)s')
    try:
        run_job(a.job, out_dir=a.out_dir, steps=a.steps, issue_utc=a.issue_utc,
                grib_dir=a.grib_dir, overwrite=a.overwrite,
                domain_km=a.domain_km, fire_dx=a.fire_dx)
    except (ForeFireJobError, fg.ForeFireGribError) as e:
        raise SystemExit('%s' % e)


if __name__ == '__main__':
    main()
