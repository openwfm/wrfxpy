#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python src/ingest/build_fmda_dict.py $1 $2 $3 $4

