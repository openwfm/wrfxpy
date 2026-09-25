#!/usr/bin/env python
"""ForeFire/WRF-SFIRE area ratio across a batch, at a common valid time.

    PYTHONPATH=src:src/ingest python src/ingest/ff_batch_ratio.py [batch_script_or_hours]

This is what every pSAF calibration number rests on, and it exists because the obvious
way to do the comparison is wrong in two ways that both inflate ForeFire.

Comparing final areas alone is wrong when the two runs end at different hours, which
they routinely do -- ForeFire chains to the last GRIB, WRF-SFIRE stops at its own
last output.  So each fire is scored at the **latest time both models reach**, found
by intersecting the two tracks, and fires where that intersection is empty are
reported rather than silently dropped.

WRF-SFIRE's curve comes from TIGN_G, which survives workspace cleaning; ForeFire's
from the chained step geojsons.
"""
from __future__ import print_function
import glob, os.path as osp, sys, math, time
sys.path[:0] = ['src', 'src/ingest']
from ff_growth import track, wrf_track

#Only score workspaces re-run in THIS batch.  Every wksp keeps its last ForeFire
#output forever, so scanning them all silently mixes the new runs with older output at
#a different pSAF and the mean means nothing.
#
#A 24 h window is NOT tight enough and produced two false outliers: Union_400257 and
#Ouachita_259 were skipped by the wrfout completeness guard, so their forefire dirs
#still held pSAF 0.6 output from the previous night -- 1438 and 1432 minutes old, just
#inside 24 h -- and they scored 4.73x and 7.77x against a set centred near 0.58x.
#Keying off the batch script's own mtime instead ties the window to the run being
#scored rather than to a guessed number of hours.
#argv[1] is the batch script whose mtime bounds the window, or a plain number of hours.
#Try the number FIRST.  Checking osp.exists first looks harmless and is not: a bare '1'
#matched a stray script named `1` sitting in the repo root, so `ff_batch_ratio.py 1`
#silently took that file's mtime -- 2025-12-04 -- and scored the entire wksp archive,
#125 fires at three different pSAF values, reported as if it were the last hour.
ARG = sys.argv[1] if len(sys.argv) > 1 else '6'
try:
    SINCE, BATCH = time.time() - float(ARG) * 3600, 'last %s h' % ARG
except ValueError:
    if not osp.exists(ARG):
        raise SystemExit('%r is neither a number of hours nor an existing path' % ARG)
    SINCE, BATCH = osp.getmtime(ARG), osp.basename(ARG)
print('scoring ForeFire output written after %s (%s)\n'
      % (time.strftime('%Y-%m-%d %H:%M', time.localtime(SINCE)), osp.basename(BATCH)))

root = '/data/jhaley/wrfxpy/wksp'
rows, skipped = [], []
for w in sorted(glob.glob(osp.join(root, 'wfc-*'))):
    #Cotton 2 and Hot Spring 226 are being rewritten with per-fuel moisture right now;
    #their pSAF-0.175 scalar-Md output was copied aside first, so read that instead of
    #a directory mid-write
    ffdir = osp.join(w, 'forefire_psaf175_scalarmd')
    if not osp.isdir(ffdir):
        ffdir = osp.join(w, 'forefire')
    if not osp.isdir(ffdir):
        continue
    gj = glob.glob(osp.join(ffdir, '*.geojson'))
    if not gj or max(osp.getmtime(f) for f in gj) < SINCE:
        continue
    try:
        ffr, _ = track(ffdir)
    except Exception as e:
        skipped.append((osp.basename(w), 'ff: %s' % e)); continue
    if not ffr:
        skipped.append((osp.basename(w), 'no ForeFire perimeters')); continue
    try:
        wr = wrf_track(w)
    except Exception as e:
        skipped.append((osp.basename(w), 'wrf: %s' % e)); continue
    if not wr:
        skipped.append((osp.basename(w), 'no TIGN_G')); continue
    ffd, wrd = dict(ffr), dict(wr)
    #WRF's TIGN_G grid is hourly; ForeFire steps land on the same minute marks, so an
    #exact intersection is the honest common time rather than nearest-neighbour fudge
    common = sorted(set(ffd) & set(wrd))
    if not common:
        skipped.append((osp.basename(w), 'no common valid time')); continue
    t = common[-1]
    name = osp.basename(w)[4:].split('_20')[0]
    rows.append((name, t, ffd[t], wrd[t], ffd[t] / wrd[t] if wrd[t] > 0 else float('nan')))

rows.sort(key=lambda r: r[4])
print('%-26s %-16s %10s %10s %8s' % ('fire', 'common valid', 'FF ha', 'WRF ha', 'ratio'))
for n, t, a, b, r in rows:
    print('%-26s %-16s %10.1f %10.1f %8.2f' % (n, t.strftime('%m-%d %H:%M'), a, b, r))

good = [r[4] for r in rows if r[4] == r[4] and r[4] > 0]
if good:
    gm = math.exp(sum(math.log(x) for x in good) / len(good))
    s = sorted(good)
    med = s[len(s)//2] if len(s) % 2 else 0.5*(s[len(s)//2-1]+s[len(s)//2])
    under = sum(1 for x in good if x < 1.0)
    within2 = sum(1 for x in good if 0.5 <= x <= 2.0)
    print('\nn=%d   geometric mean %.2f   median %.2f   min %.2f   max %.2f'
          % (len(good), gm, med, s[0], s[-1]))
    print('below 1.0x (underpredict): %d of %d (%.0f%%)' % (under, len(good), 100.0*under/len(good)))
    print('within a factor of 2:      %d of %d (%.0f%%)' % (within2, len(good), 100.0*within2/len(good)))
if skipped:
    print('\nskipped %d:' % len(skipped))
    for n, why in skipped:
        print('  %-52s %s' % (n[:52], why))
