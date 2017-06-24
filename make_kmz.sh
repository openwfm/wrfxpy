 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
python src/make_kmz.py "$*"
