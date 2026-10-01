#!/usr/bin/env python
"""Score forward ForeFire forecasts against the WRF-SFIRE run of the same fire.

    PYTHONPATH=src:src/ingest python src/ingest/ff_fwd_score.py [--hours 8] [NAME ...]

`ff_batch_ratio.py` scores the *cron* runs, which live in `<wksp>/forefire` and are
chained to WRF-SFIRE's own winds and output times.  A forward run from
`forefire_forecast.py` lives somewhere else entirely -- `wksp_forefire/<grid_code>`,
created from the job file, with no WRF-SFIRE in the loop -- so nothing in the repo
could score one.  This pairs the two by grid_code and reports both horizons:

  full   the latest valid time BOTH models reach, as ff_batch_ratio does
  8 h    the low-latency product's horizon, which is the number that actually
         matters for an 8-hour forecast and is not the full-window number: skill
         degrades with lead, so a batch calibrated at 27 h is mis-calibrated at 8

The two horizons are scored from one run, not two.  A separate 8 h run would use the
same winds and the same moisture for its first 8 steps and could only differ by noise.
"""
from __future__ import print_function
import argparse, glob, math, os, os.path as osp, sys
from datetime import timedelta
sys.path[:0] = ['src', 'src/ingest']
from ff_growth import track, wrf_track

WKSP = '/data/jhaley/wrfxpy/wksp'
FWD = '/data/jhaley/wrfxpy/wksp_forefire'


def _geom_mean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float('nan')


def _median(xs):
    s = sorted(xs)
    if not s:
        return float('nan')
    return s[len(s) // 2] if len(s) % 2 else 0.5 * (s[len(s) // 2 - 1] + s[len(s) // 2])


def _at(ffd, wrd, t):
    return ffd[t], wrd[t], (ffd[t] / wrd[t] if wrd[t] > 0 else float('nan'))


def score(grid_code, hours=8.0, tol_min=30):
    """(full_row, window_row, note) for one forward run, either row None if unscorable."""
    ffdir = osp.join(FWD, grid_code)
    #run_job writes the perimeters into the workspace root; --out-dir puts them in
    #results/.  Take whichever holds geojsons so both layouts score.
    if not glob.glob(osp.join(ffdir, '*.geojson')) and osp.isdir(osp.join(ffdir, 'results')):
        ffdir = osp.join(ffdir, 'results')
    ffr, _ = track(ffdir)
    if not ffr:
        return None, None, 'no ForeFire perimeters'
    #the WRF-SFIRE workspace carries the grid_code plus its own start time and length
    cand = sorted(glob.glob(osp.join(WKSP, 'wfc-%s-*' % grid_code)))
    if not cand:
        return None, None, 'no WRF-SFIRE workspace'
    wr = wrf_track(cand[-1])
    if not wr:
        return None, None, 'no TIGN_G'
    #Snap ForeFire to the hour before intersecting.  track() rounds to the minute,
    #which is right for a cron run -- those are driven to WRF-SFIRE's own output times
    #and land on them exactly.  A forward run steps hourly from the issue hour, but
    #ForeFire stamps each perimeter at the second the front step actually landed on, so
    #the times drift to :01 or :02 and an exact intersection throws those hours away.
    #It is not a small loss: SOUTH_PARK_RX matched 10 hours of 28 and had no 8 h point
    #at all, LAC-352011 17 and STARR_33 19, while fires whose steps happened to land on
    #:00 matched 26.  Which hours survive is then an artefact of rounding, not of the
    #forecast.
    ffd = {}
    for t, a in ffr:
        h = (t + timedelta(minutes=30)).replace(minute=0, second=0, microsecond=0)
        #keep the perimeter nearest the hour if two round onto it
        if h not in ffd or abs((t - h).total_seconds()) < abs((ffd[h][1] - h).total_seconds()):
            ffd[h] = (a, t)
    ffd = dict((h, v[0]) for h, v in ffd.items())
    wrd = dict(wr)
    common = sorted(set(ffd) & set(wrd))
    if not common:
        return None, None, 'no common valid time'

    #same truncation guard as ff_batch_ratio: a run that died two hours in still yields
    #a ratio, and it looks entirely normal because both sides are read at the same
    #early hour.  Count only what track() reads -- it skips the 'final' geojsons.
    n_steps = len([f for f in glob.glob(osp.join(ffdir, '*.geojson'))
                   if 'final' not in osp.basename(f)])
    if n_steps and len(ffr) < 0.5 * n_steps:
        return None, None, ('truncated: %d perimeters from %d step files'
                            % (len(ffr), n_steps))

    full = _at(ffd, wrd, common[-1]) + (common[-1],)
    #the window is measured from ForeFire's own first perimeter, which is the ignition
    #it was given, not from WRF-SFIRE's start_utc -- those differ by about an hour and
    #an hour is a sixth of the horizon being reported
    want = ((ffr[0][0] + timedelta(minutes=30)).replace(minute=0, second=0, microsecond=0)
            + timedelta(hours=hours))
    near = [t for t in common if abs((t - want).total_seconds()) <= tol_min * 60]
    win = None
    if near:
        t = min(near, key=lambda x: abs((x - want).total_seconds()))
        win = _at(ffd, wrd, t) + (t,)
    return full, win, None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('names', nargs='*', help='grid codes, or a substring of one; '
                   'default every run in wksp_forefire')
    p.add_argument('--hours', type=float, default=8.0, help='short horizon (default 8)')
    a = p.parse_args()

    codes = sorted(os.listdir(FWD))
    if a.names:
        codes = [c for c in codes if any(n in c for n in a.names)]

    rows, skipped = [], []
    for c in codes:
        if not osp.isdir(osp.join(FWD, c)):
            continue
        try:
            full, win, why = score(c, hours=a.hours)
        except Exception as e:
            skipped.append((c, '%s: %s' % (type(e).__name__, e)))
            continue
        if why:
            skipped.append((c, why))
            continue
        rows.append((c.split('_20')[0], full, win))

    rows.sort(key=lambda r: r[1][2])
    print('%-30s %-14s %9s %9s %7s   %-14s %9s %9s %7s'
          % ('fire', 'full valid', 'FF ha', 'WRF ha', 'ratio',
             '%gh valid' % a.hours, 'FF ha', 'WRF ha', 'ratio'))
    for name, full, win in rows:
        line = ('%-30s %-14s %9.1f %9.1f %7.2f'
                % (name[:30], full[3].strftime('%m-%d %H:%M'), full[0], full[1], full[2]))
        if win:
            line += ('   %-14s %9.1f %9.1f %7.2f'
                     % (win[3].strftime('%m-%d %H:%M'), win[0], win[1], win[2]))
        else:
            line += '   %-14s %9s %9s %7s' % ('(no %gh)' % a.hours, '-', '-', '-')
        print(line)

    for label, idx in (('full window', 1), ('%g h' % a.hours, 2)):
        good = [r[idx][2] for r in rows if r[idx] and r[idx][2] == r[idx][2] and r[idx][2] > 0]
        if not good:
            continue
        print('\n%-11s n=%d  geometric mean %.2f  median %.2f  min %.2f  max %.2f'
              % (label, len(good), _geom_mean(good), _median(good), min(good), max(good)))
        print('%-11s below 1.0x %d of %d (%.0f%%)   within 2x %d of %d (%.0f%%)'
              % ('', sum(1 for x in good if x < 1.0), len(good),
                 100.0 * sum(1 for x in good if x < 1.0) / len(good),
                 sum(1 for x in good if 0.5 <= x <= 2.0), len(good),
                 100.0 * sum(1 for x in good if 0.5 <= x <= 2.0) / len(good)))
    if skipped:
        print('\nskipped %d:' % len(skipped))
        for n, why in skipped:
            print('  %-50s %s' % (n[:50], why))


if __name__ == '__main__':
    main()
