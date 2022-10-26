"""
This file is for averaging the simulation characters in different tolerances
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path 

# find the paths of all the csv files
f_paths = []
# rtols = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
# atols = [1.0e-20, 1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12, 1.0e-10]
rtols = [1.0e-9, 1.0e-8, 1.0e-7]
atols = [1.0e-18, 1.0e-16, 1.0e-14]

for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
    sens_data = np.zeros((7001, 16))
    denominator = 0
    for i in range(len(rtols)):
        data_path = f'base_data/rtol_{rtols[i]}_atol_{atols[i]}/ref_{ratio}.csv'
        if os.path.exists(data_path):
            df = pd.read_csv(data_path)
            if len(df) >= 7001:
                sens_data += df
                denominator += 1
            elif len(df) >= 2001:
                sens_data[0:2001] += df[0:2001]
                denominator += 1
    if denominator != 0:
        sens_data /= denominator

    output_dir = Path(f"base_data/{ratio}_avg_ref.csv")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    table_names = ['Ratio', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
                   'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
                   'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv']
    pd.DataFrame(sens_data, columns=table_names).to_csv(output_dir, index=False)
