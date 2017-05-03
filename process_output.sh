 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
python src/process_output.py $1
./make_kmz.sh $1
