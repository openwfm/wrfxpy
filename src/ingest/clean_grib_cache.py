#!/usr/bin/env python
"""Prune cached forecast GRIBs by age.

    python src/ingest/clean_grib_cache.py [SOURCE ...] [--days 30] [--apply]

Dry run unless --apply is given, matching src/cleanup_fmda.py.

`cache_grib_files.py` downloads the 00/06/12/18 cycles out to +33 h so a forecast can
be built without waiting on a download. Nothing ever removed them again: `ingest/` has
no retention of any kind, and the HRRRA cache alone had reached 108 GB before this was
written. HRRR turns that from untidy into urgent -- 34 hourly leads at ~420 MB is
~14 GB per cycle and ~56 GB a day at four cycles, where NAM218's 3-hourly leads were a
fraction of that.

**Age alone is the rule, and that is a deliberate choice.** Per JH these files are
archived on AWS and can be re-downloaded if a forecast ever has to be rerun, so losing
a month-old GRIB costs download time and nothing else.

The previous version instead refused to delete any file still symlinked from a
workspace's `wps/<SOURCE>/GRIBFILE.???`. That sounds safer and is a trap at this scale:
workspaces are never cleaned up either, so every old run pins its GRIBs forever and the
sweep reclaims steadily less until it reclaims nothing. The check survives as
`--keep-linked` for the case where an offline rebuild has to be guaranteed, but it is
not the default.

Also replaces that version's mechanics, which shelled out to `ls -lah | grep grib2`,
wrote a `grib_list.txt` into the working directory, and called `rm` immediately with no
dry run.
"""

from __future__ import absolute_import
from __future__ import print_function

import argparse
import glob
import os
import os.path as osp
import sys
import time


#Where each source's files live under the ingest directory.  HRRR nests a grid
#directory ('conus') that the NAM sources do not, hence the extra level.  The key is
#also the name of the per-source staging directory under <wksp>/*/wps/.
SOURCES = {
    'HRRR':    'HRRR/hrrr.*/*/*.grib2',
    'HRRRA':   'HRRRA/hrrr.*/*/*.grib2',
    'HRRR_AK': 'HRRR_AK/hrrr.*/*/*.grib2',
    'NAM218':  'NAM218/nam.*/*.grib2',
    'NAM198':  'NAM198/nam.*/*.grib2',
    'NAM196':  'NAM196/nam.*/*.grib2',
    'NAM227':  'NAM227/nam.*/*.grib2',
}


def human(n):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024.0:
            return '%.1f %s' % (n, unit)
        n /= 1024.0
    return '%.1f PB' % n


def linked_files(wksp_root, source):
    """Absolute paths of every cache file a workspace still symlinks to.

    One glob over `<wksp>/*/wps/<SOURCE>/*` rather than a walk of every workspace:
    there are thousands of them and a naive scan takes minutes.  Broken links resolve
    to a path simply absent from the cache, which is harmless here.
    """
    out = set()
    for link in glob.glob(osp.join(wksp_root, '*', 'wps', source, '*')):
        try:
            out.add(osp.realpath(link))
        except OSError:
            continue
    return out


def sweep(ingest_root, wksp_root, source, max_age_s, now, keep_linked=False):
    """(doomed, kept_recent, kept_linked) for one source."""
    cache = glob.glob(osp.join(ingest_root, SOURCES[source]))
    if not cache:
        return [], 0, 0
    in_use = linked_files(wksp_root, source) if keep_linked else set()
    doomed, recent, linked = [], 0, 0
    for path in cache:
        try:
            age = now - osp.getmtime(path)
        except OSError:
            continue
        if age <= max_age_s:
            recent += 1
        elif keep_linked and osp.realpath(path) in in_use:
            linked += 1
        else:
            doomed.append(path)
    return doomed, recent, linked


def main():
    p = argparse.ArgumentParser(
        description='Prune cached forecast GRIBs older than --days. '
                    'Dry run unless --apply is given.')
    p.add_argument('sources', nargs='*', default=None,
                   help='sources to sweep (default: every one present in the cache). '
                        'Known: ' + ', '.join(sorted(SOURCES)))
    p.add_argument('--days', type=float, default=30.0,
                   help='delete files older than this many days (default: 30). These '
                        'are re-downloadable from the AWS archive.')
    p.add_argument('--ingest', default='ingest', help='ingest directory (default: ingest)')
    p.add_argument('--wksp', default='wksp', help='workspace root (default: wksp)')
    p.add_argument('--keep-linked', action='store_true',
                   help='also spare any file still symlinked from a workspace wps '
                        'directory. Off by default: workspaces are never pruned either, '
                        'so this eventually pins the whole cache.')
    p.add_argument('--apply', action='store_true',
                   help='actually delete. Without it, nothing is removed.')
    args = p.parse_args()

    sources = args.sources or [s for s in sorted(SOURCES)
                               if osp.isdir(osp.join(args.ingest, s))]
    unknown = [s for s in sources if s not in SOURCES]
    if unknown:
        raise SystemExit('unknown source(s): %s\nknown: %s'
                         % (', '.join(unknown), ', '.join(sorted(SOURCES))))
    if not sources:
        print('no known GRIB sources found under %s' % args.ingest)
        return

    now = time.time()
    max_age_s = args.days * 24 * 3600
    grand_doomed, grand_bytes = [], 0
    print('older than %g days, under %s\n' % (args.days, osp.abspath(args.ingest)))
    print('%-9s %9s %9s %9s %12s' % ('source', 'delete', 'recent', 'linked', 'reclaim'))
    for s in sources:
        doomed, recent, linked = sweep(args.ingest, args.wksp, s, max_age_s, now,
                                       keep_linked=args.keep_linked)
        size = 0
        for f in doomed:
            try:
                size += osp.getsize(f)
            except OSError:
                pass
        print('%-9s %9d %9d %9d %12s' % (s, len(doomed), recent, linked, human(size)))
        grand_doomed.extend(doomed)
        grand_bytes += size

    if not grand_doomed:
        print('\nnothing to remove')
        return
    if not args.apply:
        print('\nDRY RUN: would remove %d files, reclaiming %s'
              % (len(grand_doomed), human(grand_bytes)))
        print('re-run with --apply to delete')
        return

    removed = 0
    for f in grand_doomed:
        try:
            os.remove(f)
            removed += 1
        except OSError as e:
            print('  could not remove %s: %s' % (f, e), file=sys.stderr)
    print('\nremoved %d of %d files, reclaimed %s'
          % (removed, len(grand_doomed), human(grand_bytes)))


if __name__ == '__main__':
    main()
