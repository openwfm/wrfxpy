"""Is a fire's domain connectivity-gated -- can a small parameter change bottle it up?

Per JH on HOP-PATTERSON: desert, much of the domain carries no fuel, and the three
forecasts share a common western perimeter.  The moisture difference was enough to stop
the HRRR fire reaching ground where it could run uphill, while the RTMA fire got through.
That is the Mud_Bayou failure again -- WRF-SFIRE crossed a bridge and ForeFire did not --
a **connectivity** difference, not a rate difference, and it corrupts a spread
calibration because the area change is binary rather than proportional.

Three measures off the ForeFire fuel field, which is what the model actually spreads on:

  burnable%   fraction of the domain with a burnable Anderson category
  reach%      of all burnable cells, the fraction connected to the ignition -- a fire
              confined to an island cannot reach the rest however fast it spreads
  fragility   how much that reachable area shrinks when the burnable mask is eroded by
              two cells.  A blob barely notices; a domain joined by narrow necks loses a
              large share, and it is exactly those necks a small moisture change closes.
"""
from __future__ import print_function
import sys, glob, os.path as osp, warnings
warnings.filterwarnings('ignore')
import numpy as np, netCDF4 as nc
from scipy import ndimage

import argparse


def measure(ncfile):
    """(burnable_frac, reach_frac, fragility, relief) from a ForeFire step netcdf."""
    d = nc.Dataset(ncfile)
    fuel = np.array(d['fuel'][:]).squeeze()
    alt = np.array(d['altitude'][:]).squeeze() if 'altitude' in d.variables else None
    d.close()
    burn = (fuel >= 1) & (fuel <= 13)
    if not burn.any():
        return 0.0, 0.0, 1.0, None
    ny, nx = burn.shape

    def comp(mask):
        lab, _ = ndimage.label(mask)
        c = lab[ny//2, nx//2]
        if c == 0:
            if not mask.any():
                return 0
            yy, xx = np.where(mask)
            k = np.argmin((yy-ny//2)**2 + (xx-nx//2)**2)
            c = lab[yy[k], xx[k]]
        return int((lab == c).sum())

    lab, _ = ndimage.label(burn)
    c = lab[ny//2, nx//2]
    if c == 0:
        yy, xx = np.where(burn)
        k = np.argmin((yy-ny//2)**2 + (xx-nx//2)**2)
        c = lab[yy[k], xx[k]]
    reach_mask = (lab == c)
    reach = int(reach_mask.sum())
    #Erode by two cells: a blob barely notices, a domain joined by narrow necks loses a
    #large share, and those necks are what a small moisture change closes.
    #
    #Intersect with the ORIGINAL reachable component first.  Without that, a domain whose
    #ignition cell falls out of every component after erosion picks the nearest surviving
    #one instead -- which can be a different and larger component, giving reach2 > reach
    #and a negative fragility.  That produced a minimum of -0.583 over the batch.
    er = ndimage.binary_erosion(burn, np.ones((5, 5))) & reach_mask
    reach2 = comp(er) if er.any() else 0
    relief = float(alt.max()-alt.min()) if alt is not None else None
    frag = 1.0 - reach2/max(reach, 1)
    return (float(burn.mean()), reach/max(burn.sum(), 1),
            min(max(frag, 0.0), 1.0), relief)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--wksp', default='wksp', help='workspace root')
    p.add_argument('--burnable', type=float, default=0.80,
                   help='flag domains below this burnable fraction (default 0.80)')
    p.add_argument('--reach', type=float, default=0.90,
                   help='flag when less than this fraction of burnable is reachable')
    p.add_argument('--fragility', type=float, default=0.50,
                   help='flag above this fragility (default 0.50)')
    a = p.parse_args()

    rows = []
    for w in sorted(glob.glob(osp.join(a.wksp, 'wfc-*'))):
        f = sorted(glob.glob(osp.join(w, 'forefire', 'FF_*.nc')))
        if not f:
            continue
        try:
            b, r, fr, rel = measure(f[0])
        except Exception:
            continue
        rows.append((osp.basename(w)[4:].split('_20')[0], b, r, fr, rel))
    if not rows:
        raise SystemExit('no workspaces with retained FF netcdfs under %s' % a.wksp)

    B = np.array([x[1] for x in rows]); R = np.array([x[2] for x in rows])
    F = np.array([x[3] for x in rows])
    print('%d fire domains\n' % len(rows))
    print('%-12s %8s %8s %8s %8s %8s' % ('', 'min', '10th', 'median', '90th', 'max'))
    for nm, v in (('burnable', B), ('reach', R), ('fragility', F)):
        q = np.percentile(v, [10, 50, 90])
        print('%-12s %8.3f %8.3f %8.3f %8.3f %8.3f'
              % (nm, v.min(), q[0], q[1], q[2], v.max()))

    flag = [x for x in rows if x[1] < a.burnable or x[2] < a.reach or x[3] > a.fragility]
    print('\n%d of %d flagged (burnable<%.2f or reach<%.2f or fragility>%.2f)'
          % (len(flag), len(rows), a.burnable, a.reach, a.fragility))
    print('%-30s %9s %7s %10s %s' % ('fire', 'burnable', 'reach', 'fragility', 'why'))
    for n, b, r, fr, rel in sorted(flag, key=lambda x: -x[3]):
        why = ','.join(w for w, c in (('unburnable', b < a.burnable),
                                      ('unreachable', r < a.reach),
                                      ('fragile', fr > a.fragility)) if c)
        print('%-30s %8.1f%% %6.1f%% %10.2f %s' % (n[:30], 100*b, 100*r, fr, why))


if __name__ == '__main__':
    main()
