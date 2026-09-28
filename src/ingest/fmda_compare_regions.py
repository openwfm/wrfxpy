"""RTMA-driven vs HRRR-driven FMDA, as box means over fire-sized domains.

    PYTHONPATH=src:src/ingest python src/ingest/fmda_compare_regions.py

Written 2026-09-28 to answer whether moving ForeFire's moisture from wksp_fmda/CONUS
(RTMA) to wksp_fmda/CONUS_HRRR invalidates the pSAF calibration, which was fitted on
RTMA-driven output.  Result over two weeks and two domains: the mean difference in the
quantity that drives spread is 0.2-0.4 percentage points, so the calibration travels --
but see the threshold count, which is the part that does not.

The two products are on different grids -- RTMA 1258x2145, HRRR 1038x1799 -- so there
is no cell-by-cell comparison to make.  A box mean over a fire-sized domain is the
honest comparison and is also the quantity ForeFire actually consumes: md_series takes
a radius_km mean around the ignition, not a point value.

Reports each dead class and, more usefully, the **fuel 2 Md** that the CAWFE weights
produce, because that is what gets written into the fuel table and compared against the
moisture of extinction.  A small mean difference still matters if it lands on the wrong
side of me=0.15 for one product and not the other, so the disagreement rate across that
threshold is reported separately.
"""
from __future__ import print_function
import glob, os.path as osp, re, sys, warnings
import numpy as np, netCDF4 as nc, pandas as pd
warnings.filterwarnings('ignore')
sys.path[:0] = ['src', 'src/ingest']

ROOT = '/data/jhaley/clean_wrfxpy/wrfxpy/wksp_fmda'
DOMAINS = {'Sierra / DOME':        (37.568215, -119.615585),
           'Oklahoma panhandle':   (36.75,     -101.60)}
HALF_KM = 15.0                      # 30 x 30 km box
FUEL2_W = (0.5714, 0.2857, 0.1429)  # CAWFE 1-h/10-h/100-h weights for Anderson fuel 2
ME_FUEL2 = 0.15

RE = re.compile(r'^fmda-(?P<code>.+)-(?P<date>\d{8})-(?P<hour>\d{2})\.nc$')


def box_mask(geo_file, lat, lon):
    g = nc.Dataset(geo_file)
    la = np.array(g['XLAT'][:]).squeeze(); lo = np.array(g['XLONG'][:]).squeeze()
    dy = (la - lat) * 110.54
    dx = (lo - lon) * 111.32 * np.cos(np.deg2rad(lat))
    return (np.abs(dx) <= HALF_KM) & (np.abs(dy) <= HALF_KM)


def series(region, geo, lat, lon):
    m = box_mask(geo, lat, lon)
    out = {}
    for p in sorted(glob.glob('%s/%s/*/fmda-*.nc' % (ROOT, region))):
        mt = RE.match(osp.basename(p))
        if not mt:
            continue                     # skips the -fNN forecast files
        t = pd.Timestamp(mt.group('date') + 'T' + mt.group('hour') + ':00:00')
        try:
            d = nc.Dataset(p)
            a = np.ma.filled(np.ma.masked_invalid(d['FMC_GC'][:]), np.nan)
            v = [np.nanmean(a[:, :, c][m]) for c in range(3)]
            d.close()
        except Exception:
            continue
        if all(np.isfinite(x) for x in v):
            out[t] = v
    return out

print('box %.0f x %.0f km\n' % (2 * HALF_KM, 2 * HALF_KM))
for name, (lat, lon) in DOMAINS.items():
    A = series('CONUS',      ROOT + '/CONUS-geo.nc',      lat, lon)   # RTMA
    B = series('CONUS_HRRR', ROOT + '/CONUS_HRRR-geo.nc', lat, lon)   # HRRR
    common = sorted(set(A) & set(B))
    print('=== %s  (%.3f, %.3f)   %d common hours, %s .. %s'
          % (name, lat, lon, len(common), common[0], common[-1]))
    if not common:
        continue
    a = np.array([A[t] for t in common]); b = np.array([B[t] for t in common])
    print('  %-8s %9s %9s %9s %9s %7s' % ('class', 'RTMA', 'HRRR', 'HRRR-RTMA', 'RMSdiff', 'corr'))
    for i, cl in enumerate(('1-h', '10-h', '100-h')):
        r = np.corrcoef(a[:, i], b[:, i])[0, 1]
        print('  %-8s %9.4f %9.4f %+9.4f %9.4f %7.3f'
              % (cl, a[:, i].mean(), b[:, i].mean(), (b[:, i] - a[:, i]).mean(),
                 np.sqrt(((b[:, i] - a[:, i]) ** 2).mean()), r))
    ma = a.dot(FUEL2_W); mb = b.dot(FUEL2_W)
    print('  %-8s %9.4f %9.4f %+9.4f %9.4f %7.3f'
          % ('fuel2 Md', ma.mean(), mb.mean(), (mb - ma).mean(),
             np.sqrt(((mb - ma) ** 2).mean()), np.corrcoef(ma, mb)[0, 1]))
    da = ma < ME_FUEL2; db = mb < ME_FUEL2
    print('  hours able to burn (Md < %.2f):  RTMA %d/%d,  HRRR %d/%d,  disagree %d (%.0f%%)'
          % (ME_FUEL2, da.sum(), len(da), db.sum(), len(db),
             (da != db).sum(), 100.0 * (da != db).sum() / len(da)))
    print()
