#!/usr/bin/env python
"""Score a perturbed ensemble against an observed perimeter.

    PYTHONPATH=src:src/ingest python src/ingest/ff_ensemble_score.py <grid_code> <irwin_id>

A deterministic forecast is scored by a ratio.  An ensemble is scored by **where the
observation falls in the distribution**, which is a different question and the one an
ensemble exists to answer: a spread that brackets the observation is useful even when
its median is wrong, and a spread that misses it entirely says the error is not in the
perturbed parameters.

Members live in <wksp_forefire>/<grid_code>/<tag>/, one directory per tag, written by
`sweep_job`.  The tag encodes the perturbation as pNNN_mSNNN -- pSAF in thousandths and
the moisture offset in thousandths, signed.
"""
from __future__ import print_function
import argparse, glob, json, os, os.path as osp, re, sys, warnings
from datetime import timedelta
warnings.filterwarnings('ignore')
sys.path[:0] = ['src', 'src/ingest']
import numpy as np
from shapely.geometry import shape
from shapely.ops import unary_union
from ff_growth import track, wrf_track, _utm_proj, _to_utm

FWD = '/data/jhaley/wrfxpy/wksp_forefire'
PERIMS = '/data/jhaley/wrfxpy/ngfs/perims'
TAG = re.compile(r'^p(\d{3})_m([-+]\d{3})$')


def hourly(rows):
    d = {}
    for t, a in rows:
        h = (t + timedelta(minutes=30)).replace(minute=0, second=0, microsecond=0)
        if h not in d or abs((t - h).total_seconds()) < abs((d[h][1] - h).total_seconds()):
            d[h] = (a, t)
    return dict((h, v[0]) for h, v in d.items())


def observed(irwin):
    """(area_ha, timestamp, vertex count, median vertex gap m) for the IR perimeter."""
    f = glob.glob(osp.join(PERIMS, '*%s*.geojson' % irwin))
    if not f:
        raise SystemExit('no perimeter for %s' % irwin)
    fc = json.load(open(f[0]))
    pr = fc['features'][0].get('properties', {})
    g = [shape(x['geometry']) for x in fc.get('features', []) if x.get('geometry')]
    u = unary_union([x for x in g if not x.is_empty])
    u = u if u.is_valid else u.buffer(0)
    c = u.centroid
    proj, _ = _utm_proj(c.x, c.y)
    uu = _to_utm(u, proj).buffer(0)
    polys = list(uu.geoms) if uu.geom_type == 'MultiPolygon' else [uu]
    gaps = []
    nv = 0
    for q in polys:
        xy = np.array(q.exterior.coords)
        nv += len(xy)
        d = np.sqrt(((xy[1:] - xy[:-1]) ** 2).sum(axis=1))
        gaps.extend(d[d > 0])
    import datetime as _dt
    ts = _dt.datetime.strptime(pr['poly_CreateDate'], '%Y/%m/%d %H:%M:%S')
    return uu.area / 1e4, ts, nv, float(np.median(gaps))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('name', help='fire name or grid_code substring')
    p.add_argument('irwin', help='IrwinID of the observed perimeter')
    a = p.parse_args()

    gcs = [x for x in os.listdir(FWD) if a.name in x and osp.isdir(osp.join(FWD, x))]
    if not gcs:
        raise SystemExit('no forward workspace matching %r' % a.name)
    gc = gcs[0]
    ir_ha, ir_ts, nv, gap = observed(a.irwin)
    th = (ir_ts + timedelta(minutes=30)).replace(minute=0, second=0, microsecond=0)
    print('%s' % gc.split('_20')[0])
    print('  observed %.1f ha at %s  (%d vertices, %.1f m median spacing -- %s)'
          % (ir_ha, ir_ts.strftime('%m-%d %H:%M'), nv, gap,
             'IR density' if gap < 40 else 'sparse, likely hand-drawn'))

    rows = []
    for d in sorted(glob.glob(osp.join(FWD, gc, 'p*_m*'))):
        m = TAG.match(osp.basename(d))
        if not m:
            continue
        try:
            t = hourly(track(d)[0])
        except Exception:
            continue
        if th in t:
            rows.append((int(m.group(1)) / 1000.0, int(m.group(2)) / 1000.0, t[th]))
    if not rows:
        raise SystemExit('no member reached %s' % th)

    v = np.array([r[2] for r in rows])
    q = np.percentile(v, [0, 10, 25, 50, 75, 90, 100])
    print('  %d members at %s' % (len(rows), th.strftime('%m-%d %H:%M')))
    print('     min %8.1f   p10 %8.1f   p25 %8.1f   median %8.1f' % (q[0], q[1], q[2], q[3]))
    print('     p75 %8.1f   p90 %8.1f   max %8.1f   spread %.1fx'
          % (q[4], q[5], q[6], q[6] / q[0] if q[0] > 0 else float('nan')))
    below = (v < ir_ha).sum()
    pct = 100.0 * below / len(v)
    print('  observation sits at the %.0fth percentile of the ensemble' % pct)
    if ir_ha < q[0]:
        print('  *** the observation is BELOW every member -- the error is not in these '
              'parameters ***')
    elif ir_ha > q[6]:
        print('  *** the observation is ABOVE every member ***')
    else:
        print('  the ensemble brackets the observation')

    # which axis moves it more
    ps = sorted(set(r[0] for r in rows))
    ms = sorted(set(r[1] for r in rows))
    byp = [np.median([r[2] for r in rows if r[0] == x]) for x in ps]
    bym = [np.median([r[2] for r in rows if r[1] == x]) for x in ms]
    print('  pSAF     %s' % '  '.join('%.3f:%.0f' % (x, y) for x, y in zip(ps, byp)))
    print('  moisture %s' % '  '.join('%+.3f:%.0f' % (x, y) for x, y in zip(ms, bym)))
    print('  range from pSAF %.1fx, from moisture %.1fx'
          % (max(byp) / min(byp) if min(byp) > 0 else float('nan'),
             max(bym) / min(bym) if min(bym) > 0 else float('nan')))


if __name__ == '__main__':
    main()
