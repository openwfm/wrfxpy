#!/usr/bin/env bash
if [ $# -eq 0 ]
  then
     echo usage to search grid coordinates: ./clamp2mesh.sh wrffile x y
     echo usage to add subgrid coordinates: ./clamp2mesh.sh wrffile
     exit 1
fi
cd $(dirname "$0")
PYTHONPATH=src
python src/clamp2mesh.py $*
