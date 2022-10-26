import pandas as pd
import os
import numpy as np

ref_p = 'sim_data/therm_sens_ref.csv'
df_ref = pd.read_csv(ref_p)
for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
    sens_data = []
    path = f'kinetic_sensitivity_data/{ratio}_avg_kin.csv'
    df = pd.read_csv(path)
    ref_data = df_ref.loc[df_ref['Ratio'].round(2)==ratio].iloc[:, 2:-1] # data sliced from syngas selectivity to dist to 50 ch4 conv
    sens_id = df_ref.loc[df_ref['Ratio'].round(2)==ratio].iloc[:, -1]
    for i in df.loc[:,'Reactions']:
        sp_id = df.loc[df['Reactions']==i].iloc[0, 1]
        perturbed_data = df.loc[df['Reactions']==i].iloc[:,2:]
        is_all_zero = np.all((perturbed_data == 0))
        if not is_all_zero:
            sens_value = (perturbed_data.to_numpy()[0] - ref_data.to_numpy()[0]) / ref_data.to_numpy()[0] / 0.01
        else:
            sens_value = perturbed_data.to_numpy()[0]  # leave is there if the whole line is filled with zero
        sens_value = np.insert(sens_value, 0, sp_id)
        sens_value = np.insert(sens_value, len(sens_value), sens_id)
        sens_data.append(sens_value)
    table_names = ['Reactions', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
                   'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
                   'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv', 'Half Temp Pos']
        
    out_dir = f'kinetic_sensitivity_data/{ratio}_sens_data.csv'
    pd.DataFrame(np.array(sens_data), columns=table_names).to_csv(out_dir)
