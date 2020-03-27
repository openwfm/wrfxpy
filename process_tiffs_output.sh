#!/usr/bin/env bash
if [ $# -ne 1 ]
  then
     echo usage: ./process_output.sh job_id
     exit 1
fi
cd $(dirname "$0")
export PYTHONPATH=src
# if this fails, install conda and run:
#   conda create -n gdal python=3.4 gdal netcdf4 jpeg=8 pyproj matplotlib
#   conda activate gdal
python src/process_tiffs.py $1
