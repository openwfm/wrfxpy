#!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
python src/make_job_file.py "$1"
