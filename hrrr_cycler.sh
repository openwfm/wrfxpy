#!/usr/bin/env /bin/bash
#SBATCH --job-name=cycler
#SBATCH --partition=math-alderaan
#SBATCH --output=logs/cycler_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# Set up environment
eval "$(conda shell.bash hook)"
#conda activate wrfx_ml
echo "Activating conda env: wrfx_ml"
conda activate wrfx_ml

pwd
export PYTHONPATH=src
echo "Executing: python src/hrrr_cycler.py $*"
python src/hrrr_cycler.py $*

