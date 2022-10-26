import pandas as pd
import numpy as np
import os
from pathlib import Path 

# rxn_id = int(os.getenv('SLURM_ARRAY_TASK_ID', default='0'))
rxn_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1

# total 107 reactions
for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
    sens_data = np.zeros((7001, 17))
    output_dir = Path(f"sens_trends/plot_data/{ratio}/kin_sens_trend_plot{rxn_id}.csv")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    benchmark_f = f'base_data/{ratio}_avg_ref.csv'
    sens_f = f'sens_trends/{ratio}/kin_sens_trend_avg_{rxn_id}.csv'
    if os.path.exists(benchmark_f):
        benchmark = pd.read_csv(benchmark_f)
    if os.path.exists(sens_f):
        sens = pd.read_csv(sens_f)
    if len(benchmark) >= 7000:
        # exclude Ratios and Reaction ids
        diff = sens.iloc[:,2:] - benchmark[0:len(sens)].iloc[:,1:]
        sens_data = np.divide(diff, benchmark[0:len(sens)].iloc[:,1:]) / 0.01
        ratio_col = sens[0:len(sens)].iloc[:,1]
        rxn_id_col = sens[0:len(sens)].iloc[:,0]
        sens_data.insert(0, ratio_col.name, ratio_col)
        sens_data.insert(0, rxn_id_col.name, rxn_id_col)
        sens_data[np.isinf(sens_data)] = 1e-15
        sens_data = np.nan_to_num(sens_data)

    table_names = ['Reactions', 'Ratio', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
                    'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
                    'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv']
    pd.DataFrame(sens_data, columns=table_names).to_csv(output_dir, index=False)
