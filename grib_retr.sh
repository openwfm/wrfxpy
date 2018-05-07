#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python src/ingest/retrieve_gribs.py $1 $2 $3 $4

