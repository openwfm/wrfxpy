 #!/usr/bin/env bash
cd $(dirname "$0")
PYTHONPATH=src
python src/ssh_command.py $*
