#!/usr/bin/env bash
if [ $# -eq 0 ]
  then
     echo usage: ./clamp2mesh.sh wrffile x y:
     exit 1
fi
cd $(dirname "$0")
PYTHONPATH=src
python src/clamp2mesh.py $*
