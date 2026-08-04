#!/usr/bin/env /bin/bash
#SBATCH --job-name=cycler
#SBATCH --partition=math-alderaan
#SBATCH --output=logs/cycler_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# Set up environment
eval "$(conda shell.bash hook)"
conda activate fmda_ml

pwd
export PYTHONPATH=src
python src/hrrr_cycler.py $*

