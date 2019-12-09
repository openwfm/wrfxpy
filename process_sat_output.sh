 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
if [ $# -eq 0 ]
  then
     echo usage: ./process_sat_output.sh job_id
     exit 1
fi
python src/process_sat_output.py $1
./make_kmz.sh $1
