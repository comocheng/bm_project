"""
This file is for averaging the simulation characters in different tolerances
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path 

# find the paths of all the csv files
f_paths = []
rtols = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
atols = [1.0e-20, 1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12, 1.0e-10]
            
# total 107 reactions
for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
    sens_data = np.zeros((107, 16))
    sens_data[:,0] = range(107)
    denominator = np.zeros(107)
    for i in range(len(rtols)):
        data_path = f'kinetic_sens/rtol_{rtols[i]}_atol_{atols[i]}/{ratio}/kin_sens_properties.csv'
        if os.path.exists(data_path):
            df = pd.read_csv(data_path)
            for j in df.loc[:,'Reactions']:
                increment = df.loc[df['Reactions']==j].iloc[:,2:]
                is_all_zero = np.all((increment == 0))
                if not is_all_zero:
                    denominator[int(j)] += 1
                    sens_data[int(j)][1:] += increment
        
    #avoid divided by zero error
    for i in range(len(denominator)):
        if denominator[i] == 0:
            denominator[i] = 1
            
    sens_data[:, 1:] = sens_data[:, 1:] / denominator[:,None]
    output_dir = Path(f"kinetic_sensitivity_data/{ratio}_avg_kin.csv")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    table_names = ['Reactions', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
                       'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
                       'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv']
    pd.DataFrame(sens_data, columns=table_names).to_csv(output_dir)
