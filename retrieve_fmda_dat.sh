#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python src/ingest/retrieve_fmda_dat.py $1 $2 $3 $4