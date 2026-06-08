#!/usr/bin/env /bin/bash

#SBATCH --job-name=hrrri
#SBATCH --partition=math-alderaan
#SBATCH --output=logs/hinterval_%j.out
#SBATCH --ntasks=4
#SBATCH --mem=64G

if [ "$#" -ne 1 ]; then
    echo "Error: Expected exactly 1 arguments, but got $#."
    echo "Usage: $0 <config_path>"
    exit 1
fi

pwd
CONFIG_PATH="$1"
# Set up environment
eval "$(conda shell.bash hook)"
conda activate fmda_ml

export PYTHONPATH=src
python src/hrrr_interval.py "$CONFIG_PATH"

