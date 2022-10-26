#!/bin/bash
#SBATCH --job-name=rmg_runs_1-107.sh
#SBATCH --output=logs_arr/cpox_surfaces.%a.log
#SBATCH --error=logs_arr/cpox_surfaces.%a.slurm.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-107
source activate rmg_env
# python sens_trends.py
# python sens_trends_avg.py
python sens_trends_cal.py
