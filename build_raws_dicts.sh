#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python src/ingest/build_raws_dicts.py $1 $2 $3 $4 $5 

