#!/usr/bin/env bash
cd $(dirname "$0")
path=$1
if [ "${path:0:1}" != "/" ]; then
    path=$(pwd)/$1
fi
export PYTHONPATH=src
python src/tests/test_domains.py $path
