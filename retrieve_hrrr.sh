#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python src/ingest/retrieve_hrrr.py $1 $2 $3
