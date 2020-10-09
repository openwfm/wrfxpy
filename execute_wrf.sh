#!/usr/bin/env bash
if [ $# -eq 0 ]
  then
     echo usage: ./execute_wrf.sh job_id
     exit 1
fi
cd $(dirname "$0")
export PYTHONPATH=src
python src/execute_wrf.py $1
