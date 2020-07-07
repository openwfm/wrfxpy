#!/usr/bin/env bash
if [ $# -eq 0 ]
  then
     echo usage: ./cycling.sh input.json
     exit 1
fi
cd $(dirname "$0")
PYTHONPATH=src
python src/cycling.py $1
