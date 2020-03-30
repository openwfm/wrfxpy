#!/usr/bin/env bash
cd $(dirname "$0")
file=$1
dir=$2
if [ "${file:0:1}" != "/" ]; then
    file=$(pwd)/$1
fi
if [ "${dir:0:1}" != "/" ]; then
    dir=$(pwd)/$2
fi
export PYTHONPATH=src
python src/geo/convert_geotiff.py $file $dir $3
