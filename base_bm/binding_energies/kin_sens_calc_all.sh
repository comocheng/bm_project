#!/bin/bash
#SBATCH --job-name=rmg_runs_0-81.sh
#SBATCH --output=logs/cpox_surfaces.%a.log
#SBATCH --error=logs/cpox_surfaces.%a.slurm.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-81
# Define useful bash variables
source activate rmg_env
RUN_i=$SLURM_ARRAY_TASK_ID
b_energies=()
for value1 in $(seq 5.5 0.25 7.5)
do
    for value2 in $(seq 3.25 0.25 5.25)
    do
        b_energies+=("c-${value1}o-${value2}")
    done
done
# RUN Simulations
((index=$RUN_i-1))
surf="${b_energies[$index]}"
cd "${RUN_i}.0_${surf}"
python kin_sens_calc.py
