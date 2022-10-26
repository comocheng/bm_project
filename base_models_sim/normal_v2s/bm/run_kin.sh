#!/bin/bash
#SBATCH --job-name=kin_sens
#SBATCH --output=logs_arr/kin_sens.out
#SBATCH --error=logs_arr/kin_sens.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
# Define useful bash variables
source activate rmg_env
# python kin_sens_gen.py
# python kin_data_processing.py
# python kin_avg.py
# python kin_sens_calc.py
python sim_base.py 
python sim_ref_data.py
python sim_ref_avg.py
