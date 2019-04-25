#!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
if [ $# -eq 0 ]
  then
     echo usage: ./cycling.sh job.json
     exit 1
fi
python src/cycling.py $1
