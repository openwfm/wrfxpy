"""Which fuel-moisture product reproduces the stall structure WRF-SFIRE shows?

    PYTHONPATH=src:src/ingest python src/ingest/ff_moisture_stall_test.py

Generalises the BOON case (handoff 2026-09-28 section 7) to every fire whose area moved
most when ForeFire's moisture switched from RTMA- to HRRR-driven FMDA.  Run over those
seven it found the answer is **neither**: RTMA matches on two, HRRR on one, neither on a
fourth, and each blocks spuriously once where WRF-SFIRE shows no stall at all.  BOON is
the extreme case, not a representative one.

For each fire, reconstruct the RTMA and HRRR fuel-moisture series over its own timing
table, compare each fuel's Md against that fuel's moisture of extinction, and ask which
product's blocked-step pattern matches the stall WRF-SFIRE independently shows.

md_series is deterministic given the timing table and cfg, so this works even though the
RTMA run's perimeters were overwritten.
"""
from __future__ import print_function
import sys, copy, glob, os.path as osp, warnings
warnings.filterwarnings('ignore')
sys.path[:0] = ['src', 'src/ingest']
import numpy as np, pandas as pd, netCDF4 as nc
import forefire as ff
from ff_growth import wrf_track

ME = {1:.12, 2:.15, 3:.25, 4:.2, 5:.2, 6:.25, 7:.4, 8:.3, 9:.25, 10:.25, 11:.15, 12:.2, 13:.25}
FIRES = ['Merit_Creek', 'BOON', 'Menogyn_Piles_RX', 'HOP-PATTERSON_PILE_RX',
         'D6_RX_Treasured_Landscapes', 'Little_Blue_Stem_RX', 'Silver']
DELTA = {'Merit_Creek': +93, 'BOON': +61, 'Menogyn_Piles_RX': +45,
         'HOP-PATTERSON_PILE_RX': -64, 'D6_RX_Treasured_Landscapes': -51,
         'Little_Blue_Stem_RX': -49, 'Silver': -37}

base = ff.load_ff_cfg()
CFG = {}
for nm, reg in (('RTMA', 'CONUS'), ('HRRR', 'CONUS_HRRR')):
    c = copy.deepcopy(base)
    c['moisture']['fmda_dir'] = '/data/jhaley/clean_wrfxpy/wrfxpy/wksp_fmda/' + reg
    c['moisture']['fmda_geo_file'] = '/data/jhaley/clean_wrfxpy/wrfxpy/wksp_fmda/%s-geo.nc' % reg
    CFG[nm] = c
NML = base['moisture']['class_weight_namelist']
AF = base['moisture']['class_weight_anderson_fuels']


def dominant_fuels(w, n=3):
    f = sorted(glob.glob(osp.join(w, 'forefire', 'FF_*.nc')))
    if not f:
        return []
    d = nc.Dataset(f[0]); a = np.array(d['fuel'][:]).squeeze(); d.close()
    ny, nx = a.shape
    c = a[ny//2-25:ny//2+25, nx//2-25:nx//2+25].astype(int)
    u, k = np.unique(c, return_counts=True)
    order = [(int(x), int(y)) for x, y in sorted(zip(u, k), key=lambda p: -p[1])
             if 1 <= x <= 13]
    return order[:n]


def stall_window(w):
    """(fraction of hours below a quarter of peak growth, the hours themselves)."""
    t = wrf_track(w)
    if len(t) < 4:
        return None
    a = np.array([x[1] for x in t]); hrs = [x[0] for x in t]
    g = np.diff(a)
    if g.max() <= 0:
        return None
    low = g < 0.25 * g.max()
    return dict(frac=low.mean(), n=len(g), final=a[-1],
                lo=[hrs[i+1] for i in range(len(g)) if low[i]], gmax=g.max())


print('%-27s %5s %7s %7s | %-22s | %s'
      % ('fire', 'dArea', 'WRFha', 'stall%', 'blocked steps (fuel: n)', 'verdict'))
for name in FIRES:
    w = sorted(glob.glob('wksp/wfc-%s_*' % name))
    if not w:
        print('%-27s NOT FOUND' % name); continue
    w = w[0]
    try:
        _, ign, latlon, _ = ff.read_input(w)
        tt = ff.make_timing_table(w, ign)
        if not len(tt):
            print('%-27s no timing table' % name); continue
        fuels = dominant_fuels(w)
        sw = stall_window(w)
        ser = {nm: ff.md_series(tt, latlon, c) for nm, c in CFG.items()}
    except Exception as e:
        print('%-27s %s: %s' % (name, type(e).__name__, e)); continue
    blocked = {}
    for nm in ('RTMA', 'HRRR'):
        s = ser[nm]
        parts = []
        for fu, _cnt in fuels:
            wgt = ff.dead_class_weights(fu, mode='cawfe', namelist_path=NML, anderson_fuels=AF)
            if not wgt: continue
            n = sum(1 for i, (_, r) in enumerate(tt.iterrows())
                    if r['ign_seconds'] >= 0 and s[i] is not None
                    and sum(x*y for x, y in zip(wgt, s[i])) >= ME[fu])
            parts.append('%d:%d' % (fu, n))
        blocked[nm] = ','.join(parts)
    stallpct = '%.0f%%' % (100*sw['frac']) if sw else '  -'
    wrfha = '%7.0f' % sw['final'] if sw else '      -'
    verdict = ''
    if sw:
        rb = sum(int(x.split(':')[1]) for x in blocked['RTMA'].split(',') if ':' in x)
        hb = sum(int(x.split(':')[1]) for x in blocked['HRRR'].split(',') if ':' in x)
        if sw['frac'] > 0.3 and rb > hb: verdict = 'WRF stalls; RTMA blocks, HRRR less'
        elif sw['frac'] > 0.3 and hb > rb: verdict = 'WRF stalls; HRRR blocks, RTMA less'
        elif sw['frac'] > 0.3: verdict = 'WRF stalls; neither blocks'
        else: verdict = 'no WRF stall'
    print('%-27s %+4d%% %s %7s | R %-9s H %-9s | %s'
          % (name, DELTA[name], wrfha, stallpct, blocked['RTMA'], blocked['HRRR'], verdict))
