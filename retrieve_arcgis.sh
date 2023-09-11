#!/usr/bin/env bash
source ~/.bashrc
conda activate arcgis
cd $(dirname "$0")
export PYTHONPATH=src
python src/fire_init/acq_arcgis.py $*
