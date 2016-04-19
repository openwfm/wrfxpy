#!/usr/bin/env bash
export PYTHONPATH=src
python src/ingest/grib_file.py $1 $2 $3 $4

