#!/bin/bash
#SBATCH --job-name=single_run
#SBATCH --output=slurm_output.log
#SBATCH --error=slurm_err.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=20
# Define useful bash variables
source activate rmg_env
PYTHONPATH=/home/xu.chao/cantera/build/python python simulation.py