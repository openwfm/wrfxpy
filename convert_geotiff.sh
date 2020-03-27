#!/usr/bin/env bash
if [ $# -ne 3 ]
  then
    echo usage: $0 geotiff_file geogrid_folder var_name
    echo example: $0 ./fuel.tif ./fuel NFUEL_CAT
    exit 1
fi
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
python src/tif/convert_geotiff.py $file $dir $3
