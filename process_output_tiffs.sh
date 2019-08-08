 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
if [ $# -eq 0 ]
  then
     echo usage: ./process_output.sh job_id
     exit 1
fi
# if this fails, install conda and run:
#   conda create -n gdal python=3.4 gdal netcdf4 jpeg=8 pyproj matplotlib
conda activate gdal
python src/process_tiffs.py $1
