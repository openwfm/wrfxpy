#!/usr/bin/env bash
source /home/jhaley/.bashrc
conda activate wrf_test
export PYTHONPATH=src

# One ForeFire run at a time.  run_forecasts() deletes every file in the shared
# run_dir on entry, and that directory is one fixed path for every fire, so a
# second invocation landing on top of a running one wipes its netcdfs mid-step.
# ForeFire then fails with "wrong input file, check your settings...", writes a
# restart holding FireDomain t=0 with no nodes, and every later step in the chain
# inherits no fire.  Seen on Union_400257: cron fired at 15:47, 16:07 and 16:27 on
# 2026-09-20, the chain broke at the second, and 26 of 47 ensemble members came out
# empty.  A 47-step fire takes longer than the 20 minute cron interval, so the
# overlap was not a rare race.
#
# Runs once a day at 00:10 rather than every 20 minutes.  The every-20-minutes
# schedule existed before anyone needed it that often, and it caused real damage: it
# wiped an in-flight run's netcdfs on Union (2026-09-20), and on 2026-09-23 it
# regenerated a workspace's .ff scripts in the middle of a controlled experiment.
# Per JH the actual requirement is only to have the previous day's forecasts ready
# to look at on arrival, which one nightly run satisfies.
#
# The flock stays: a nightly run can still overlap a manual one, and that is exactly
# the collision that caused the damage above.  -n means a tick that finds a run still
# going simply skips, matching the flock the FMDA cyclers use.  Two details that
# matter:
#   - the python is no longer backgrounded.  With '&' the shell exited immediately
#     and flock released the lock straight away, so the lock would do nothing.
#   - forefire.log is truncated inside the lock, so a skipped tick cannot destroy
#     the log of the run that is still writing it.
/usr/bin/flock -n /tmp/cron_forefire.lock \
    bash -c 'exec python src/ingest/forefire.py &> forefire.log' \
  || echo "$(date -u +%FT%TZ) skipped, previous ForeFire run still in progress" >> forefire_skipped.log
