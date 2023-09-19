#!/bin/bash
#SBATCH --job-name=rmg_runs_0-107.sh
#SBATCH --output=logs/cpox_surfaces.%a.log
#SBATCH --error=logs/cpox_surfaces.%a.slurm.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-107
# Define useful bash variables
source activate rmg_env
RUN_i=$SLURM_ARRAY_TASK_ID
dirs=()

for value2 in $(seq 0 1 106)
do
    dirs+=("rxn${value2}")
done

# RUN Simulations
((index=$RUN_i-1))
d="${dirs[$index]}"
cd "${d}"
#PYTHONPATH=/home/xu.chao/cantera/build/python python simulation.py
python simulation.py