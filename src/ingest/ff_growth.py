#!/usr/bin/env python
"""Fire area against time, for one or more ForeFire runs side by side.

An ensemble comparison reports one number per run at one valid time, which hides
*when* two runs diverged. This reads the chained step track -- the non-final
geojsons, which are the sequence of perimeters through the forecast -- and reports
area and growth rate hour by hour.

It exists for two questions that a single endpoint cannot answer:

  - **Does growth stall?** WRF-SFIRE stalls as fuel moisture climbs overnight.
    ForeFire could not reproduce that until Md was allowed to vary per step, because
    a constant Md of 0.1 sits below the 0.12 moisture of extinction and the fire can
    always spread. A stall shows up here as a run of near-zero growth rates.
  - **Where did two runs part company?** Two runs can reach similar totals by
    different paths, and the hour they separate usually says why.

Areas are computed in the local UTM zone, so they are real areas rather than
latitude-scaled approximations, and every run is projected through the same zone so
the columns are comparable.

    PYTHONPATH=src:src/ingest python src/ingest/ff_growth.py \\
        label=<dir> [label=<dir> ...] [--pattern '*.geojson']
"""

from __future__ import print_function

import argparse
import glob
import json
import os.path as osp
from datetime import datetime, timedelta
import re

import numpy as np
import netCDF4 as nc
from shapely.geometry import shape
from shapely.ops import unary_union

from ff_score_perim import _to_utm, _utm_proj, _valid_at


def track(model_dir, pattern='*.geojson'):
    """{valid_at: area_ha} for a run's chained step perimeters, in a UTM zone
    chosen from the first perimeter seen."""
    rows = []
    proj = None
    for f in sorted(glob.glob(osp.join(model_dir, pattern))):
        base = osp.basename(f)
        if 'final' in base:
            continue
        fc = json.load(open(f))
        geoms = [shape(x['geometry']) for x in fc.get('features', []) if x.get('geometry')]
        geoms = [g for g in geoms if not g.is_empty]
        if not geoms:
            continue
        g = unary_union(geoms)
        if not g.is_valid:
            g = g.buffer(0)
        t = _valid_at(fc)
        if t is None:
            continue
        #ForeFire stamps a perimeter at whatever second the front step landed on, so
        #two runs of the same fire differ by a few seconds at the same nominal time.
        #Snapping to the minute puts them on one row instead of two half-empty ones.
        t = (t + timedelta(seconds=30)).replace(second=0, microsecond=0)
        if proj is None:
            c = g.centroid
            proj, _ = _utm_proj(c.x, c.y)
        gu = _to_utm(g, proj)
        if not gu.is_valid:
            gu = gu.buffer(0)
        rows.append((t, gu.area / 1e4))
    rows.sort()
    return rows, proj


def wrf_track(wksp_dir, step_minutes=60, cell_m=None):
    """WRF-SFIRE's own burned area through time, from TIGN_G.

    `TIGN_G` in the **retained final wrfout** is each fire-grid cell's ignition time
    in seconds since the simulation start, with a large sentinel for cells that never
    burned.  Area at time t is therefore just the count of cells whose TIGN_G is at
    or below t, which recovers the whole progression from **one file**.

    That matters because cleaned workspaces lose `FIRE_AREA` -- their wrfouts are
    replaced by saveouts, which do not carry it -- so summing FIRE_AREA per output
    time only works on fires that have not been cleaned yet.  TIGN_G works on both,
    and evaluates at any instant rather than only at output times.

    Validated on Hot Spring 226, where both survive: 442.8 ha against FIRE_AREA's
    440.9 at the final time, tracking at every hour.  It runs **0.5-1% high** because
    TIGN_G marks a cell burned at its ignition instant while FIRE_AREA ramps
    fractionally -- a consistent bias, not noise.

    Returns [(datetime, area_ha)].
    """
    outs = sorted(glob.glob(osp.join(wksp_dir, 'wrf', 'wrfout_d01_*')))
    if not outs:
        return []
    d = nc.Dataset(outs[-1])
    if 'TIGN_G' not in d.variables:
        return []
    #drop the edge strip make_FF_nc drops, so areas match the ForeFire domain
    sr = int(d.dimensions['west_east_subgrid'].size /
             (d.dimensions['west_east'].size + 1))
    t = np.array(d['TIGN_G'][0, :-sr, :-sr])
    if cell_m is None:
        cell_m = float(d.DX) / sr
    cell_ha = cell_m * cell_m / 1e4
    burned = t < t.max() * 0.999          # the rest carry the never-burned sentinel
    if not burned.any():
        return []

    m = re.search(r'-(\d{4}-\d{2}-\d{2})_(\d{2}):(\d{2}):(\d{2})-\d+$',
                  osp.basename(wksp_dir.rstrip('/')))
    if m:
        start = datetime.strptime('%s %s:%s:%s' % m.groups(), '%Y-%m-%d %H:%M:%S')
    else:
        start = datetime.strptime(osp.basename(outs[0])[11:], '%Y-%m-%d_%H:%M:%S')

    end = float(t[burned].max())
    rows = []
    step = step_minutes * 60
    k = step
    while k <= end + step:
        rows.append((start + timedelta(seconds=k),
                     float(((t <= k) & burned).sum()) * cell_ha))
        k += step
    return rows


def compare(runs, pattern='*.geojson', stall_ha_per_h=1.0, wrf=None):
    """runs is [(label, dir), ...].  Prints area and growth rate per hour.

    wrf is [(label, workspace), ...] read from TIGN_G via wrf_track, so WRF-SFIRE's
    own fire appears as a column beside the ForeFire forecasts.
    """
    tracks = {}
    proj = None
    for label, d in runs:
        r, p = track(d, pattern)
        if proj is None:
            proj = p
        tracks[label] = dict(r)
        print('%-14s %3d perimeters  %s -> %s' % (
            label, len(r), r[0][0] if r else '-', r[-1][0] if r else '-'))
    for label, w in (wrf or []):
        r = wrf_track(w)
        tracks[label] = dict(r)
        print('%-14s %3d steps from TIGN_G  %s -> %s' % (
            label, len(r), r[0][0] if r else '-', r[-1][0] if r else '-'))
        runs = list(runs) + [(label, w)]
    times = sorted({t for v in tracks.values() for t in v})
    if not times:
        raise SystemExit('no step perimeters found; try --pattern')

    labels = [l for l, _ in runs]
    print('\n%-20s' % 'valid', end='')
    for l in labels:
        print('%12s' % l, end='')
    print('   |  growth ha/h')
    prev = {}
    for t in times:
        print('%-20s' % t.strftime('%Y-%m-%d %H:%M'), end='')
        rates = []
        for l in labels:
            a = tracks[l].get(t)
            print('%12s' % ('%.1f' % a if a is not None else '-'), end='')
            if a is not None and l in prev:
                dt = (t - prev[l][0]).total_seconds() / 3600.0
                rates.append((a - prev[l][1]) / dt if dt > 0 else float('nan'))
            else:
                rates.append(float('nan'))
            if a is not None:
                prev[l] = (t, a)
        print('   | ' + '  '.join(
            ('%7.1f' % r if np.isfinite(r) else '      -') for r in rates), end='')
        #flag hours where a run has effectively stopped growing
        stalled = [labels[i] for i, r in enumerate(rates)
                   if np.isfinite(r) and r < stall_ha_per_h]
        print('   <= stalled: ' + ','.join(stalled) if stalled else '')
    print()
    for l in labels:
        v = tracks[l]
        if len(v) < 2:
            continue
        ts = sorted(v)
        hrs = (ts[-1] - ts[0]).total_seconds() / 3600.0
        rate = [(v[b] - v[a]) / ((b - a).total_seconds() / 3600.0)
                for a, b in zip(ts, ts[1:]) if (b - a).total_seconds() > 0]
        n_stall = sum(1 for r in rate if r < stall_ha_per_h)
        print('%-14s final %8.1f ha over %4.1f h   mean %6.1f ha/h   peak %6.1f   '
              'stalled %d of %d hours'
              % (l, v[ts[-1]], hrs, (v[ts[-1]] - v[ts[0]]) / hrs if hrs else 0,
                 max(rate) if rate else 0, n_stall, len(rate)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs', nargs='+', help='label=/path/to/forefire/output')
    p.add_argument('--pattern', default='*.geojson')
    p.add_argument('--stall', type=float, default=1.0,
                   help='growth below this many ha/h counts as stalled (default 1)')
    a = p.parse_args()
    runs = []
    wrf = []
    for spec in a.runs:
        if '=' not in spec:
            raise SystemExit('expected label=dir, got %r' % spec)
        label, d = spec.split('=', 1)
        #a 'wrf:' prefix reads WRF-SFIRE's own fire from TIGN_G instead of ForeFire
        #geojsons, so the reference curve sits in the same table as the forecasts
        if label.startswith('wrf:'):
            wrf.append((label[4:] or 'WRF-SFIRE', d))
        else:
            runs.append((label, d))
    compare(runs, a.pattern, a.stall, wrf=wrf)


if __name__ == '__main__':
    main()
