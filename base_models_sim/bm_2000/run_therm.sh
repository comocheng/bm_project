#!/bin/bash
#SBATCH --job-name=therm_sens
#SBATCH --output=logs/therm_sens.out
#SBATCH --error=logs/therm_sens.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c[5003,3040]
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
# Define useful bash variables
source activate rmg_env
# python therm_sens_gen.py
python data_processing.py
python thermo_avg.py
python thermo_sens_calc.py
