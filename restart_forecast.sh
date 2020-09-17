#!/usr/bin/env bash
if [ $# -eq 0 ]
  then
     echo usage: ./restart_forecast.sh job_id
     exit 1
fi
cd $(dirname "$0")
export PYTHONPATH=src
python src/restart_forecast.py $1
