#!/bin/bash
#SBATCH --job-name=single_run
#SBATCH --output=output.log
#SBATCH --error=err.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
# Define useful bash variables
source activate rmg_env
python data_processing.py