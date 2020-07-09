 #!/usr/bin/env bash
if [ $# -ne 1 ]
  then
     echo usage: ./process_sat_output.sh job_id
     exit 1
fi
cd $(dirname "$0")
PYTHONPATH=src
python src/process_sat_output.py $1
./make_kmz.sh $1
