#!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
python src/vis/csv2kml.py $1 doc.kml
/bin/rm -f $2
zip -9 $2 doc.kml
