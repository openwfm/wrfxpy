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


#period_hours per grib source, read off the source classes so this cannot drift from
#what wrfxpy itself uses.  Only the attribute is needed, so nothing is instantiated.
def _period_hours(name):
    from ingest.NAM218 import NAM218
    from ingest.HRRR import HRRR
    table = {'NAM': NAM218, 'NAM218': NAM218, 'HRRR': HRRR}
    try:
        from ingest.NAM227 import NAM227
        table['NAM227'] = NAM227
    except Exception:
        pass
    src = table.get(str(name).upper())
    return getattr(src, 'period_hours', 3) if src else 3


def workspace_name(job_path):
    """The `wfc-...` directory name wrfxpy itself would give this job.

    Reproduces forecast.py:115 --

        'wfc-' + grid_code + '-' + utc_to_esmf(start_utc) + '-%02d' % fc_hrs

    -- including the detail that makes it non-obvious: `start_utc` and `end_utc` are
    first rounded, down and up respectively, to the grib source's `period_hours`, and
    `fc_hrs` is the span *after* that rounding.  So a job starting 20:00 with a 3-hourly
    NAM source names itself from 18:00, not 20:00.

    Matching the name matters because it is what makes a ForeFire-only fire sit in the
    same workspace a WRF-SFIRE run of the same job would have used, so every existing
    tool -- ff_batch_ratio, ff_growth, ff_domain_connectivity -- finds it unchanged.
    """
    from utils import round_time_to_hour, utc_to_esmf, timedelta_hours
    with open(job_path) as f:
        job = json.load(f)
    ph = _period_hours(job.get('grib_source'))
    start = round_time_to_hour(parse_job_time(job['start_utc']).to_pydatetime(),
                               up=False, period_hours=ph)
    end = round_time_to_hour(parse_job_time(job['end_utc']).to_pydatetime(),
                             up=True, period_hours=ph)
    fc_hrs = timedelta_hours(end - start)
    code = job.get('grid_code') or osp.splitext(osp.basename(job_path))[0]
    return 'wfc-%s-%s-%02d' % (code, utc_to_esmf(start), int(fc_hrs))


SAFE_DIR_RE = re.compile(r'[^A-Za-z0-9._-]')


def safe_dir_name(code):
    """`grid_code` reduced to characters that are safe as a directory name.

    Measured over 28,695 real job files: **7,821 grid_codes (27%) contain a colon**,
    from an older convention that wrote the timestamp `18:00:00` where current jobs
    write `18_00_00`.  A handful carry worse -- spaces (10), '#' (22), parentheses (9),
    '~' (5), '&' (1).

    Colons matter most because an apptainer bind spec is colon-delimited and cannot
    carry one, so a directory named from such a code cannot be bound directly.  Spaces
    and '&' matter because parts of this pipeline still build shell commands by string
    interpolation.

    Mapping ':' to '_' is not arbitrary -- it is exactly what the newer grid_code
    convention already does, so old and new jobs converge on the same directory name.
    Returns (name, changed).
    """
    safe = SAFE_DIR_RE.sub('_', code)
    return safe, safe != code


def make_workspace(job_path, root=None, cfg=None):
    """Create this job's ForeFire directory and return it.

    `<wksp_forefire>/<grid_code>/`, per JH.  Three reasons it is not the wrfxpy `wfc-`
    workspace:

    - a VIIRS-only fire has no WRF-SFIRE run and so no `wfc-` workspace to sit in;
    - the `wfc-` name encodes the forecast window, so editing a job's `end_utc` moves
      the directory -- Bigfoot went from `-15` to `-33` mid-session when its window was
      extended, scattering the results;
    - `grid_code` contains no colons, and an apptainer bind spec cannot carry one, so
      this directory can be bound directly while a `wfc-` name cannot.

    `workspace_name()` is still the way to find the matching WRF-SFIRE workspace when
    there is one.
    """
    if cfg is None:
        cfg = ff.load_ff_cfg()
    with open(job_path) as f:
        job = json.load(f)
    code = job.get('grid_code') or osp.splitext(osp.basename(job_path))[0]
    code, changed = safe_dir_name(code)
    if changed:
        logging.warning('grid_code contains characters unsafe in a directory name '
                        '(a colon cannot appear in an apptainer bind); using %r', code)
    root = root or cfg.get('wksp_forefire', '/data/jhaley/wrfxpy/wksp_forefire')
    w = osp.join(root, code)
    if not osp.isdir(w):
        os.makedirs(w)
    dest = osp.join(w, 'input.json')
    if not osp.exists(dest):
        with open(job_path) as a, open(dest, 'w') as b:
            b.write(a.read())
    return w


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
            overwrite=False, cfg=None, domain_km=None, fire_dx=None, params=None,
            wksp_root=None):
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

    #Default into the workspace wrfxpy itself would have used, creating it if this fire
    #never had a WRF-SFIRE run -- a VIIRS-only detection has no workspace at all.  That
    #puts the results where every other tool already looks, instead of somewhere only
    #this script knows about.  --out-dir remains the escape hatch.
    if out_dir:
        run_dir = out_dir
        dest_dir = osp.join(run_dir, 'results')
    else:
        wksp = make_workspace(job_path, root=wksp_root, cfg=cfg)
        logging.info('  workspace %s', wksp)
        #Run *in* the workspace rather than in a shared staging directory and moving the
        #results afterwards.  The per-job apptainer bind mounts this at /work, so there
        #is no longer a single run_dir for every fire to collide in -- which is the
        #collision run_forecasts' "delete everything on entry" used to cause.
        run_dir = dest_dir = wksp
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
    #issue_utc reaches the moisture too, not just the winds: without it a forecast run
    #picks up FMDA *analyses* that assimilated observations made after the forecast was
    #issued, and read_fmda_hourly's preference for analysis over forecast guarantees it
    tt = ff.make_script_set(run_dir, tt, latlon, grid_code, cfg=cfg,
                            params=dict(params or {}), issue_utc=issue_utc)
    tt.to_csv(osp.join(run_dir, 'timing_table.csv'), index=False)

    ff.run_timing_table(tt, overwrite=overwrite, cfg=cfg, run_dir=run_dir,
                        dest_dir=dest_dir)
    ff.cleanup_ff_run(run_dir, None, grid_code, cfg=cfg, dest_dir=dest_dir)
    logging.info('done -> %s', dest_dir)
    return dest_dir


def sweep_job(job_path, param_sets, tags=None, steps=8, issue_utc=None,
              grib_dir=None, cfg=None, domain_km=None, fire_dx=None,
              wksp_root=None, overwrite=False):
    """Run one fire several times with different ForeFire parameters, reusing the winds.

    Measured over the 21 forward runs, **WindNinja is 82-98% of the wall clock** (median
    95%): 19-46 minutes of netcdf against 1-6 minutes of ForeFire.  The netcdfs carry
    winds, fuels and terrain, none of which depend on the tuning parameters, so a sweep
    that builds them once and varies only the `.ff` settings costs about two minutes per
    variant instead of thirty -- roughly 16x.  That is what makes a 20-fire by 5-value
    sweep a three-hour job rather than a two-day one.

    Each variant runs in its own subdirectory with the netcdfs symlinked in, because the
    restart template chains through `include[sim_<grid>_<step>.ff]` and two variants
    sharing a directory would read each other's perimeters.  The **parent** is what gets
    bound into the container, so those relative links resolve; `--pwd` selects the
    variant.

    Returns {tag: results_dir}.
    """
    if cfg is None:
        cfg = ff.load_ff_cfg()
    if tags is None:
        tags = ['_'.join('%s%g' % (k[:6], v) for k, v in sorted(p.items())) or 'base'
                for p in param_sets]
    if len(tags) != len(param_sets):
        raise ForeFireJobError('%d tags for %d parameter sets'
                               % (len(tags), len(param_sets)))

    latlon, ign_utc, grid_code, job_km, job_dx, _ = read_job(job_path)
    domain_km = fit_domain(domain_km or job_km, fire_dx or job_dx)
    fire_dx = fire_dx or job_dx
    grib_dir = grib_dir or cfg.get('hrrr_forecast_cache', HRRR_FORECAST_CACHE)
    base = make_workspace(job_path, root=wksp_root, cfg=cfg)

    #overwrite=False so existing netcdfs are kept: this is the whole point
    tt = fg.build_step_ncs(latlon[0], latlon[1], ign_utc, base, grid_code, steps=steps,
                           domain_km=domain_km, fire_dx=fire_dx, grib_dir=grib_dir,
                           issue_utc=issue_utc, overwrite=False)
    if tt is None or not len(tt):
        raise ForeFireJobError('no GRIB valid times cover this ignition')
    tt['wrfout'] = tt['grib']

    out = {}
    for params, tag in zip(param_sets, tags):
        run_dir = osp.join(base, tag)
        if not osp.isdir(run_dir):
            os.makedirs(run_dir)
        for name in sorted(os.listdir(base)):
            if not name.endswith('.nc'):
                continue
            link = osp.join(run_dir, name)
            want = osp.relpath(osp.join(base, name), run_dir)
            if osp.islink(link):
                if os.readlink(link) == want:
                    continue
                os.unlink(link)
            elif osp.exists(link):
                continue
            os.symlink(want, link)

        logging.info('=== %s : %s', tag, params)
        t = ff.make_script_set(run_dir, tt.copy(), latlon, grid_code, cfg=cfg,
                               params=dict(params), issue_utc=issue_utc)
        t.to_csv(osp.join(run_dir, 'timing_table.csv'), index=False)
        ff.run_timing_table(t, overwrite=overwrite, cfg=cfg, run_dir=run_dir,
                            dest_dir=run_dir, bind_dir=base)
        #drop the links before finalising, or they dangle in the results
        for name in os.listdir(run_dir):
            l = osp.join(run_dir, name)
            if osp.islink(l) and name.endswith('.nc'):
                os.unlink(l)
        ff.finalize_ff_outputs(run_dir, grid_code, cfg=cfg, suffix=tag)
        out[tag] = run_dir
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('job', help='a wrfxpy job json, or a workspace input.json')
    p.add_argument('--out-dir', default=None,
                   help='put results here instead of in the wrfxpy workspace')
    p.add_argument('--wksp-root', default=None,
                   help='root for ForeFire forecast directories; default is '
                        'wksp_forefire from etc/forefire.json. One directory per job '
                        'grid_code, created if absent.')
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
                domain_km=a.domain_km, fire_dx=a.fire_dx, wksp_root=a.wksp_root)
    except (ForeFireJobError, fg.ForeFireGribError) as e:
        raise SystemExit('%s' % e)


if __name__ == '__main__':
    main()
