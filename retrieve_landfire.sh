#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python -u src/ingest/retrieve_landfire.py $*
