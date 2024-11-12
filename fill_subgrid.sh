 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
if [ $# -eq 0 ]
  then
     echo usage: ./fill_subgrid.sh wrf_path
     exit 1
fi
python src/fill_subgrid.py $1
