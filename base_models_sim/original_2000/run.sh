#!/bin/bash
#SBATCH --job-name=normal_sim
#SBATCH --output=logs/normal_sim.out
#SBATCH --error=logs/normal_sim.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
# Define useful bash variables
source activate rmg_env
python simulation.py
# python sim_base.py
# python sim_ref_data.py
# python sim_ref_avg.py