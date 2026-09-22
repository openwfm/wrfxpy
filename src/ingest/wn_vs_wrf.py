#!/usr/bin/env python
"""Compare a WindNinja-derived ForeFire wind field against WRF-SFIRE's UF/VF.

Both fields carry the per-fuel `windrf` reduction -- the wrfout's by construction,
the ForeFire netcdf's because `forefire_grib.apply_windrf` put it there -- so they
are directly comparable without dividing anything out.

The two live on different grids: WRF's fire grid is Lambert, the netcdf's is UTM.
WRF's fire-grid cell centres (FXLONG/FXLAT) are projected into the netcdf's UTM and
the netcdf is sampled bilinearly there, so every statistic is computed on WRF's grid.

Reports, per hour, on cells that are burnable, above 0.5 m/s in both fields and
inside the netcdf domain:

    R         circular resultant length of the WRF field.  Below ~0.3 a direction
              comparison cannot distinguish anything and the direction columns
              should be ignored -- the trap recorded in the 2026-09-18 handoff.
    B/A       mean speed ratio, WindNinja over WRF
    bias/RMS  direction difference, degrees
    >90       fraction of cells whose direction differs by more than 90 degrees
    anomR     correlation of the speed anomalies, i.e. whether the spatial
              structure agrees once the means are removed
    peakB/A   ratio of peak speeds, which is where terrain acceleration shows up
"""

from __future__ import print_function

import argparse
import glob
import os
import os.path as osp
import re

import numpy as np
import netCDF4 as nc
from pyproj import Proj

import forefire_grib as fg

_UTC_RE = re.compile(r'_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)\.nc$')


def _bilinear(A, fx, fy):
    ny, nx = A.shape
    x0 = np.floor(fx).astype(int)
    y0 = np.floor(fy).astype(int)
    ok = (x0 >= 0) & (y0 >= 0) & (x0 < nx - 1) & (y0 < ny - 1)
    xs = np.clip(x0, 0, nx - 2)
    ys = np.clip(y0, 0, ny - 2)
    wx = fx - xs
    wy = fy - ys
    v = (A[ys, xs] * (1 - wx) * (1 - wy) + A[ys, xs + 1] * wx * (1 - wy)
         + A[ys + 1, xs] * (1 - wx) * wy + A[ys + 1, xs + 1] * wx * wy)
    return v, ok


def compare(wksp, nc_dir, lat, lon, domain_km=30.0, fire_dx=25.0, min_speed=0.5):
    d0 = nc.Dataset(sorted(glob.glob(osp.join(wksp, 'wrf', 'wrfout_d01_*')))[0])
    sr = int(d0.dimensions['west_east_subgrid'].size / (d0.dimensions['west_east'].size + 1))
    lonf = np.array(d0['FXLONG'][0, :-sr, :-sr])
    latf = np.array(d0['FXLAT'][0, :-sr, :-sr])
    fuel = np.array(nc.Dataset(osp.join(wksp, 'wrf', 'wrfinput_d01'))['NFUEL_CAT'][0, :-sr, :-sr])
    burn = (fuel >= 1) & (fuel <= 13)

    dom = fg.build_domain(lat, lon, domain_km=domain_km, fire_dx=fire_dx)
    proj = Proj(init='epsg:{}'.format(dom['epsg']))
    X, Y = proj(lonf, latf)
    fx = (X - dom['x0']) / dom['dx'] - 0.5
    fy = (Y - dom['y0']) / dom['dx'] - 0.5

    rows = []
    print('%-18s %6s %6s %7s %7s %7s %7s %8s %9s'
          % ('valid', 'R', 'B/A', 'bias', 'dirRMS', '>90deg', 'anomR', 'peakB/A', 'cells'))
    for f in sorted(glob.glob(osp.join(nc_dir, 'FF_*.nc'))):
        m_utc = _UTC_RE.search(f)
        if not m_utc:
            continue
        utc = m_utc.group(1)
        wf = osp.join(wksp, 'wrf', 'wrfout_d01_' + utc[:-1].replace('T', '_'))
        if not osp.exists(wf):
            continue
        d = nc.Dataset(wf)
        if 'UF' not in d.variables:
            continue
        UF = np.array(d['UF'][0, :-sr, :-sr])
        VF = np.array(d['VF'][0, :-sr, :-sr])
        g = nc.Dataset(f)
        wu, ok1 = _bilinear(np.array(g['windU'][0, 0]), fx, fy)
        wv, ok2 = _bilinear(np.array(g['windV'][0, 0]), fx, fy)

        sA = np.hypot(UF, VF)
        sB = np.hypot(wu, wv)
        m = ok1 & ok2 & burn & (sA > min_speed) & (sB > min_speed) & np.isfinite(sB)
        if m.sum() < 1000:
            continue
        th = np.arctan2(UF[m], VF[m])
        R = float(np.hypot(np.cos(th).mean(), np.sin(th).mean()))
        dd = (np.degrees(np.arctan2(wu, wv)) - np.degrees(np.arctan2(UF, VF)) + 180) % 360 - 180
        a = sA[m] - sA[m].mean()
        b = sB[m] - sB[m].mean()
        anom = float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))
        row = dict(utc=utc[:16], R=R, ba=float(sB[m].mean() / sA[m].mean()),
                   bias=float(dd[m].mean()), rms=float(np.sqrt((dd[m] ** 2).mean())),
                   gt90=float((np.abs(dd[m]) > 90).mean()), anom=anom,
                   peak=float(sB[m].max() / sA[m].max()), n=int(m.sum()))
        rows.append(row)
        print('%-18s %6.3f %6.3f %+7.1f %7.1f %6.1f%% %7.3f %8.2f %9d'
              % (row['utc'], R, row['ba'], row['bias'], row['rms'],
                 100 * row['gt90'], anom, row['peak'], row['n']))

    if rows:
        good = [r for r in rows if r['R'] >= 0.3]
        print('\n%d hours, %d with R >= 0.3 (direction meaningful)' % (len(rows), len(good)))
        for key, label, fmt in (('ba', 'B/A', '%.3f'), ('rms', 'dir RMS', '%.1f deg'),
                                ('anom', 'anomaly corr', '%.3f'), ('peak', 'peak B/A', '%.2f')):
            src = good if key in ('rms',) else rows
            v = np.array([r[key] for r in src])
            print('  %-14s mean ' % label + fmt % v.mean()
                  + '   range ' + fmt % v.min() + ' to ' + fmt % v.max())
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('wksp', help='wrfxpy workspace holding wrf/wrfout_d01_*')
    p.add_argument('nc_dir', help='directory of FF_*.nc built by forefire_grib.py')
    p.add_argument('lat', type=float)
    p.add_argument('lon', type=float)
    p.add_argument('--domain-km', type=float, default=30.0)
    p.add_argument('--fire-dx', type=float, default=25.0)
    p.add_argument('--min-speed', type=float, default=0.5)
    a = p.parse_args()
    compare(a.wksp, a.nc_dir, a.lat, a.lon, a.domain_km, a.fire_dx, a.min_speed)


if __name__ == '__main__':
    main()
