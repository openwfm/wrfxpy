#!/usr/bin/env python
"""Score ForeFire perimeters against an observed IR perimeter.

Replaces `dry_ff_compare.py`, which the 2026-09-18 session wrote in a scratchpad and
lost; its open item 6 named it as worth rewriting into the repo.

Metrics, following that session's choices:

    area        hectares, computed in the local UTM zone rather than on lon/lat, so
                it is a real area and not a latitude-scaled approximation
    obs x       modelled area / observed area
    IoU         intersection over union against the observation
    elongation  long side / short side of the minimum-area rotated rectangle

**Elongation, not reach.** Per JH the estimated ignition point in an observation
file is always suspect, and reach measured from it inherits that error. Elongation
is a property of the polygon alone, so it is the shape metric that survives a bad
ignition. The 09-18 session established this after finding `reach ~ pSAF^1` only
measurable where the origin is trustworthy.

The rotated-rectangle definition is stated because 09-18's own implementation is
gone. It does, however, **reproduce that session's numbers**: on the Dry River IR
perimeter this tool gives area 2078.5 ha and elongation 2.90 against the handoff's
recorded 2078.5 ha and 2.90, and the area is 100.0% of the stated acreage. So
elongation here is comparable to the 09-18 and 09-11 figures, not merely
self-consistent.

**Suppression.** An estimated 98% of fires are fought, and neither model represents
that, so `obs_x > 1` on a suppressed fire is expected rather than a model error. Where
`attr_ContainmentDateTime` falls inside the run, pass it to `--at`: scoring at the
perimeter's own timestamp otherwise measures however many hours the model kept
growing after the real fire was held. Even the containment-time score stays biased,
because suppression acts throughout rather than only at the end.

Observed perimeters are NIFC/IRWIN geojson; the observation time is
`poly_PolygonDateTime`. ForeFire geojsons carry `valid_at`. Modelled perimeters are
grouped by `valid_at` and the group nearest the observation time is scored, each
member separately so an ensemble reports a distribution rather than one number.

    PYTHONPATH=src:src/ingest python src/ingest/ff_score_perim.py \\
        <observed.geojson> <dir of ForeFire geojsons> [--pattern '*final*.geojson']
"""

from __future__ import print_function

import argparse
import glob
import json
import os.path as osp
from datetime import datetime

import numpy as np
from pyproj import Proj
from shapely.geometry import shape
from shapely.ops import unary_union


def _load(path):
    """Geojson -> a single shapely geometry in lon/lat, plus its properties."""
    fc = json.load(open(path))
    geoms = [shape(f['geometry']) for f in fc.get('features', [])
             if f.get('geometry')]
    geoms = [g for g in geoms if not g.is_empty]
    if not geoms:
        return None, fc
    g = unary_union(geoms)
    if not g.is_valid:
        g = g.buffer(0)
    return g, fc


def _utm_proj(lon, lat):
    zone = int((lon + 180.0) / 6.0) + 1
    epsg = (32600 if lat >= 0 else 32700) + zone
    return Proj(init='epsg:{}'.format(epsg)), epsg


def _to_utm(geom, proj):
    """Reproject a (multi)polygon's coordinates with a pyproj 1.x Proj."""
    from shapely.geometry import Polygon, MultiPolygon

    def ring(coords):
        xs, ys = proj([c[0] for c in coords], [c[1] for c in coords])
        return list(zip(xs, ys))

    def poly(p):
        return Polygon(ring(p.exterior.coords),
                       [ring(i.coords) for i in p.interiors])

    if geom.geom_type == 'Polygon':
        return poly(geom)
    return MultiPolygon([poly(p) for p in geom.geoms])


def elongation(geom):
    """Long side / short side of the minimum-area rotated rectangle."""
    mrr = geom.minimum_rotated_rectangle
    if mrr.geom_type != 'Polygon':
        return float('nan')
    c = np.array(mrr.exterior.coords)[:4]
    sides = [float(np.hypot(*(c[(i + 1) % 4] - c[i]))) for i in range(4)]
    a, b = sides[0], sides[1]
    if min(a, b) <= 0:
        return float('nan')
    return max(a, b) / min(a, b)


def _obs_time(props):
    for k in ('poly_PolygonDateTime', 'poly_DateCurrent'):
        v = props.get(k)
        if v:
            for fmt in ('%Y/%m/%d %H:%M:%S', '%Y-%m-%dT%H:%M:%SZ'):
                try:
                    return datetime.strptime(v, fmt)
                except ValueError:
                    pass
    return None


def _valid_at(fc):
    v = fc.get('valid_at')
    if not v:
        return None
    for fmt in ('%Y-%m-%dT%H:%M:%SZ', '%Y-%m-%d %H:%M:%S'):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            pass
    return None


def score(obs_path, model_dir, pattern='*final*.geojson', at=None):
    obs, obs_fc = _load(obs_path)
    if obs is None:
        raise SystemExit('observed perimeter has no geometry: %s' % obs_path)
    props = obs_fc['features'][0].get('properties', {})
    t_obs = _obs_time(props)
    #Scoring target, which is not always the observation time.  On a suppressed fire
    #the model keeps growing after the real fire was held, so scoring at the
    #perimeter's own timestamp measures that overrun rather than the forecast.  Pass
    #--at with attr_ContainmentDateTime to score where the comparison is fairest.
    t_target = at or t_obs
    c = obs.centroid
    proj, epsg = _utm_proj(c.x, c.y)
    obs_u = _to_utm(obs, proj)
    obs_ha = obs_u.area / 1e4

    print('observed  %s' % osp.basename(obs_path))
    print('  time %s   area %.1f ha   elongation %.2f   (EPSG:%d)'
          % (t_obs, obs_ha, elongation(obs_u), epsg))
    stated = props.get('poly_GISAcres')
    if stated:
        print('  stated %.1f acres = %.1f ha  -> computed is %.1f%% of stated'
              % (stated, stated * 0.404686, 100.0 * obs_ha / (stated * 0.404686)))

    #group candidates by valid_at, keep the group nearest the observation
    cand = {}
    for f in sorted(glob.glob(osp.join(model_dir, pattern))):
        g, fc = _load(f)
        if g is None:
            continue
        t = _valid_at(fc)
        cand.setdefault(t, []).append((osp.basename(f), g))
    if not cand:
        raise SystemExit('no non-empty perimeters matching %s in %s'
                         % (pattern, model_dir))
    times = [t for t in cand if t is not None]
    if times and t_target is not None:
        t_pick = min(times, key=lambda t: abs((t - t_target).total_seconds()))
        dt_min = (t_pick - t_target).total_seconds() / 60.0
    else:
        t_pick = list(cand)[0]
        dt_min = float('nan')
    members = cand[t_pick]
    print('\nmodelled  %d perimeters at %s  (%+.0f min from target %s)'
          % (len(members), t_pick, dt_min, t_target))

    rows = []
    print('  %-40s %9s %7s %7s %7s' % ('member', 'area_ha', 'obs_x', 'IoU', 'elong'))
    for name, g in members:
        gu = _to_utm(g, proj)
        if not gu.is_valid:
            gu = gu.buffer(0)
        a = gu.area / 1e4
        inter = gu.intersection(obs_u).area
        union = gu.union(obs_u).area
        iou = inter / union if union > 0 else float('nan')
        e = elongation(gu)
        rows.append((a, a / obs_ha, iou, e))
        print('  %-40s %9.1f %7.3f %7.3f %7.2f' % (name[:40], a, a / obs_ha, iou, e))

    if len(rows) > 1:
        arr = np.array(rows)
        print('\n  ensemble  area %.1f-%.1f ha (mean %.1f)   obs_x %.3f-%.3f (mean %.3f)'
              % (arr[:, 0].min(), arr[:, 0].max(), arr[:, 0].mean(),
                 arr[:, 1].min(), arr[:, 1].max(), arr[:, 1].mean()))
        print('            IoU %.3f-%.3f (mean %.3f)   elongation %.2f-%.2f (mean %.2f)'
              % (arr[:, 2].min(), arr[:, 2].max(), arr[:, 2].mean(),
                 arr[:, 3].min(), arr[:, 3].max(), arr[:, 3].mean()))
        print('  observed elongation %.2f' % elongation(obs_u))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('observed', help='NIFC/IRWIN perimeter geojson')
    p.add_argument('model_dir', help='directory of ForeFire geojson output')
    p.add_argument('--pattern', default='*final*.geojson',
                   help="glob for modelled perimeters (default '*final*.geojson'; "
                        "use '*.geojson' to include the chained step track)")
    p.add_argument('--at', default=None,
                   help='score the modelled perimeters nearest this UTC time instead '
                        'of the observation time, e.g. the fire\'s '
                        'attr_ContainmentDateTime. Format YYYY-MM-DDTHH:MM:SS')
    a = p.parse_args()
    at = datetime.strptime(a.at, '%Y-%m-%dT%H:%M:%S') if a.at else None
    score(a.observed, a.model_dir, a.pattern, at=at)


if __name__ == '__main__':
    main()
