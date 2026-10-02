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


def _perimeter(path, out_tol=0.01):
    """One step's burned area as a geometry, with ForeFire's rings read correctly.

    `shape()` on the raw geometry is WRONG and silently so.  ForeFire's
    dumpMode=geojson packs every separate front as a RING of one polygon, and the
    GeoJSON convention makes ring 0 the shell and the rest holes -- so a second front
    is SUBTRACTED from the first.  The step where one appears reads as a sudden shrink
    and every step after it is wrong too: the cron run of SAN_FRANCISCO_2 read 64.2 ha,
    then 36.3, and ended at 132.4 against a true 1621.8.

    `forefire.read_ff_geojson` classifies rings by nesting depth under the even-odd
    rule -- even is burning, odd is an unburned island.  That is right for the KML and
    not enough for an area, because **not every nested ring is a hole**.  It decides
    containment from a representative point, so a ring that merely overlaps the
    perimeter is taken for a hole inside it.  Two different artefacts land in that trap:

      - a front dumped as a thin band, where the inner edge reads as an island nearly
        as large as the burn (1621.8 against 1607.4 ha on SAN_FRANCISCO_2, a mean width
        of 7.25 m on a 25 m mesh);
      - a **second front ignited partly within the first**, whose boundary reads as a
        hole.  Hatch_Grande step 27 grew its outer ring from 1493.3 to 1601.1 ha and
        gained a 147.7 ha "island", so the measured area FELL from 1483.8 to 1457.8.
        Most of that 147.7 ha was ground the perimeter already enclosed one step
        earlier, and burned ground does not become unburned.

    The discriminator is geometric: **a hole lies entirely inside the ring it
    punctures.**  Measured on both artefacts and on the genuine islands beside them,
    real holes sit at 0.000% of their area outside the burn while both artefacts poke
    out by over 8%.  So a nested ring is differenced out only when it is enclosed, and
    unioned in otherwise.  `out_tol` is slack for reprojection slivers, far below the
    8% the artefacts show.

    This replaces an area-ratio threshold used briefly on 2026-10-01, which caught the
    thin-band case and not the second-front case.  Containment is an invariant rather
    than a tuned constant and covers both.
    """
    import forefire as _ff
    gdf = _ff.read_ff_geojson(path)
    if len(gdf) == 0:
        return None
    def _clean(g):
        return g if g.is_valid else g.buffer(0)
    burn = [_clean(g) for g in gdf[gdf['burning']].geometry
            if g is not None and not g.is_empty]
    nest = [_clean(g) for g in gdf[~gdf['burning']].geometry
            if g is not None and not g.is_empty]
    if not burn and not nest:
        return None

    g = _clean(unary_union(burn)) if burn else None
    #a ring that is not enclosed by the burn is not a hole in it -- it is another front
    holes, fronts = [], []
    for x in nest:
        inside = x.intersection(g).area if g is not None else 0.0
        (holes if x.area > 0 and 1.0 - inside / x.area <= out_tol else fronts).append(x)
    if fronts:
        g = _clean(unary_union(([g] if g is not None else []) + fronts))
    if g is None:
        return None
    if holes:
        g = _clean(g.difference(_clean(unary_union(holes))))
    return None if g.is_empty else g


def shrink_events(rows, tol=0.01):
    """Steps where the burned area fell, as [(step, when, before, after, frac)].

    **A fire cannot shrink.**  That makes a falling track the one check on a ForeFire
    area that needs no reference data, and it is the check that was missing: the ring
    subtraction fixed in `_perimeter` ran through a 65-fire pSAF calibration, three
    handoffs and two scoring tools before anyone noticed, because an under-measured
    fire looks exactly like a fire that did not spread.  SAN_FRANCISCO_2 lost 44% of
    its area in one step and was scored as a 0.19 ratio without comment.

    `tol` is a fraction of the preceding area.  Real runs do lose a little at a chained
    restart, where ForeFire re-meshes the perimeter it was handed -- measured over 57
    runs the only survivors of the ring fix were Hatch_Grande at 1.75% and
    Lauderdale at 1.41%, both healthy.  So small losses are reported and large ones
    disqualify the run; see SHRINK_SKIP in the scorers.
    """
    out = []
    for i in range(1, len(rows)):
        prev, cur = rows[i - 1][1], rows[i][1]
        if prev > 0 and cur < prev * (1.0 - tol):
            out.append((i + 1, rows[i][0], prev, cur, 1.0 - cur / prev))
    return out


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
        g = _perimeter(f)
        if g is None:
            continue
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
