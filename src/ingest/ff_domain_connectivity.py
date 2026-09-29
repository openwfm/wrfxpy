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

FIRES = ['Merit_Creek', 'BOON', 'Menogyn_Piles_RX', 'HOP-PATTERSON_PILE_RX',
         'D6_RX_Treasured_Landscapes', 'Little_Blue_Stem_RX', 'Silver', 'Mud_Bayou']
DELTA = {'Merit_Creek': +93, 'BOON': +61, 'Menogyn_Piles_RX': +45,
         'HOP-PATTERSON_PILE_RX': -64, 'D6_RX_Treasured_Landscapes': -51,
         'Little_Blue_Stem_RX': -49, 'Silver': -37, 'Mud_Bayou': None}

print('%-27s %6s %9s %8s %10s %9s' %
      ('fire', 'dArea', 'burnable%', 'reach%', 'fragility', 'relief m'))
for name in FIRES:
    w = sorted(glob.glob('wksp/wfc-%s_*' % name))
    if not w:
        print('%-27s not found' % name); continue
    f = sorted(glob.glob(osp.join(w[0], 'forefire', 'FF_*.nc')))
    if not f:
        print('%-27s no FF netcdf retained' % name); continue
    d = nc.Dataset(f[0])
    fuel = np.array(d['fuel'][:]).squeeze()
    alt = np.array(d['altitude'][:]).squeeze() if 'altitude' in d.variables else None
    d.close()
    burn = (fuel >= 1) & (fuel <= 13)
    ny, nx = burn.shape
    lab, _ = ndimage.label(burn)
    ig = lab[ny//2, nx//2]
    if ig == 0:                      # ignition on non-fuel: take the nearest burnable
        yy, xx = np.where(burn)
        k = np.argmin((yy-ny//2)**2 + (xx-nx//2)**2)
        ig = lab[yy[k], xx[k]]
    reach = (lab == ig)
    er = ndimage.binary_erosion(burn, np.ones((5, 5)))
    lab2, _ = ndimage.label(er)
    ig2 = lab2[ny//2, nx//2]
    if ig2 == 0 and er.any():
        yy, xx = np.where(er)
        k = np.argmin((yy-ny//2)**2 + (xx-nx//2)**2)
        ig2 = lab2[yy[k], xx[k]]
    reach2 = (lab2 == ig2).sum() if ig2 else 0
    frag = 1.0 - reach2 / max(reach.sum(), 1)
    dz = '%9.0f' % (alt.max()-alt.min()) if alt is not None else '        -'
    dl = '%+5d%%' % DELTA[name] if DELTA[name] is not None else '    -'
    print('%-27s %6s %8.1f%% %7.1f%% %9.2f %s'
          % (name, dl, 100*burn.mean(), 100*reach.sum()/max(burn.sum(), 1), frag, dz))
