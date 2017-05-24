#!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
python src/vis/csv2kml.py "$1" "$2"
