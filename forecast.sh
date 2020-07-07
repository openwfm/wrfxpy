#!/usr/bin/env bash
if [ $# -eq 0 ]
  then
     echo usage: ./forecast.sh input.json
     exit 1
fi
cd $(dirname "$0")
PYTHONPATH=src
python src/forecast.py $1
