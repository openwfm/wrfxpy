 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
if [ $# -eq 0 ]
  then
     echo usage: ./send_to_server.sh job_id
     exit 1
fi
python src/send_to_server.py $1
