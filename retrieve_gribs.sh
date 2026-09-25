#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
# "$@" rather than $1 $2 $3 $4: the optional 5th argument (cycle_start) would
# otherwise be dropped silently, and the caller would get an unpinned cycle.
python src/ingest/retrieve_gribs.py "$@"

