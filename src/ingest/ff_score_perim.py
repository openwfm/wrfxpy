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
import re
import os.path as osp
from datetime import datetime
from email.utils import parsedate_to_datetime

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
    """Observation time, from either date format NIFC ships.

    Per-fire files use '2026/09/15 16:09:00'; the year-to-date collections use
    RFC-822, '"Wed, 16 Sep 2026 13:31:00 GMT"'.  Both appear in ngfs/perims.
    """
    for k in ('poly_PolygonDateTime', 'poly_DateCurrent'):
        v = props.get(k)
        if not v:
            continue
        for fmt in ('%Y/%m/%d %H:%M:%S', '%Y-%m-%dT%H:%M:%SZ'):
            try:
                return datetime.strptime(v, fmt)
            except ValueError:
                pass
        try:
            return parsedate_to_datetime(v).replace(tzinfo=None)
        except Exception:
            pass
    return None


def extract_from_ytd(path, irwin):
    """Pull one fire's feature out of a year-to-date perimeter collection.

    `ngfs/perims/perims_ytd_<date>.geojson` holds every fire's latest perimeter in a
    single ~170 MB line, so these are read as text and the one feature is cut out
    rather than parsing the whole collection.

    Two traps.  `poly_IRWINID` is `"{AACFD673-...}"` — **literal braces inside a
    string** — so naive brace matching starts from the wrong place and never closes;
    the feature start is found with rfind on '{"type":"Feature"' instead, and the
    forward scan tracks whether it is inside a string.  And the dates here are
    RFC-822, not the per-fire files' '%Y/%m/%d' (handled in _obs_time).
    """
    with open(path) as fh:
        txt = fh.read()
    i = txt.find(irwin)
    if i < 0:
        return None
    j = txt.rfind('{"type":"Feature"', 0, i)
    if j < 0:
        return None
    depth = 0
    k = j
    instr = False
    esc = False
    while k < len(txt):
        c = txt[k]
        if esc:
            esc = False
        elif c == '\\':
            esc = True
        elif c == '"':
            instr = not instr
        elif not instr:
            if c == '{':
                depth += 1
            elif c == '}':
                depth -= 1
                if depth == 0:
                    break
        k += 1
    try:
        return json.loads(txt[j:k + 1])
    except ValueError:
        return None


def ytd_series(ytd_dir, irwin, since=None, until=None):
    """Every distinct perimeter for one fire across the ytd collections.

    Returns [(time, feature, source_file)] sorted by time, de-duplicated on
    poly_PolygonDateTime, since consecutive daily files usually repeat the same
    perimeter until a new one is flown.
    """
    out = {}
    for f in sorted(glob.glob(osp.join(ytd_dir, 'perims_ytd_*.geojson'))):
        base = osp.basename(f)
        if since or until:
            #cheap date filter on the filename before reading 170 MB
            m = re.search(r'perims_ytd_(\d{4}-\d{2}-\d{2})', base)
            if m:
                d = m.group(1)
                if since and d < since:
                    continue
                if until and d > until:
                    continue
        feat = extract_from_ytd(f, irwin)
        if not feat:
            continue
        t = _obs_time(feat.get('properties', {}))
        if t is None or t in out:
            continue
        out[t] = (feat, base)
    return [(t, out[t][0], out[t][1]) for t in sorted(out)]


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
    p.add_argument('observed', nargs='?', default=None,
                   help='NIFC/IRWIN perimeter geojson; omit when using --irwin')
    p.add_argument('--irwin', default=None,
                   help='IRWIN id (with or without braces) to pull from the '
                        'year-to-date collections instead of a single file')
    p.add_argument('--ytd-dir', default='/data/jhaley/wrfxpy/ngfs/perims',
                   help='directory holding perims_ytd_*.geojson')
    p.add_argument('--since', default=None, help='earliest ytd file date, YYYY-MM-DD')
    p.add_argument('--until', default=None, help='latest ytd file date, YYYY-MM-DD')
    p.add_argument('--list', action='store_true',
                   help='with --irwin, list the perimeters found and stop')
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

    obs_path = a.observed
    if a.irwin:
        irwin = a.irwin.strip('{}')
        series = ytd_series(a.ytd_dir, irwin, a.since, a.until)
        if not series:
            raise SystemExit('no perimeters for %s in %s' % (irwin, a.ytd_dir))
        if a.list:
            for t, feat, src in series:
                pr = feat.get('properties', {})
                print('%s  %10s acres  %s' % (t, pr.get('poly_GISAcres'), src))
            return
        #pick the perimeter nearest the scoring target, so an in-window observation
        #is used rather than whichever file happened to be read last
        target = at or max(t for t, _, _ in series)
        t, feat, src = min(series, key=lambda r: abs((r[0] - target).total_seconds()))
        print('using the %s perimeter from %s (%d of %d available)\n'
              % (t, src, [r[0] for r in series].index(t) + 1, len(series)))
        import tempfile
        fh = tempfile.NamedTemporaryFile('w', suffix='.geojson', delete=False)
        json.dump({'type': 'FeatureCollection', 'features': [feat]}, fh)
        fh.close()
        obs_path = fh.name
    if not obs_path:
        raise SystemExit('give an observed geojson or --irwin')
    score(obs_path, a.model_dir, a.pattern, at=at)


if __name__ == '__main__':
    main()
