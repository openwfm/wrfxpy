 #!/usr/bin/env bash
PYTHONPATH=src
python src/forecast.py $1
./make_kmz.sh $1
