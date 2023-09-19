"""
This file is for averaging the simulation characters in different tolerances
The script order to generate the thermo sensitivity data is sim_sens_gen.py -> data_processing.py -> thermo_avg.py -> thermo_sens_calc.py
"""

# import pandas as pd
# import numpy as np
# import os
# from pathlib import Path 

# # find the paths of all the csv files
# f_paths = []
# rtols = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
# atols = [1.0e-20, 1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12, 1.0e-10]
            
# # total 20 species
# for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
#     sens_data = np.zeros((20,16))
#     sens_data[:,0] = range(20)
#     denominator = np.zeros(20)
#     for i in range(len(rtols)):
#         data_path = f'rtol_{rtols[i]}_atol_{atols[i]}/{ratio}/thermo_sens_properties.csv'
#         if os.path.exists(data_path):
#             df = pd.read_csv(data_path)
#             for j in df.loc[:,'Species']:
#                 increment = df.loc[df['Species']==j].iloc[:, 2:]
#                 is_all_zero = np.all((increment == 0))
#                 if not is_all_zero:
#                     denominator[int(j)] += 1
#                     sens_data[int(j)][1:] += increment
        
#     #avoid divided by zero error
#     for i in range(len(denominator)):
#         if denominator[i] == 0:
#             denominator[i] = 1
            
#     sens_data[:, 1:] = sens_data[:, 1:] / denominator[:,None]
#     output_dir = Path(f"thermo_sensitivity_data/{ratio}_avg_thermo.csv")
#     output_dir.parent.mkdir(parents=True, exist_ok=True)
#     table_names = ['Species', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
#                        'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
#                        'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv']
#     pd.DataFrame(sens_data, columns=table_names).to_csv(output_dir)

import pandas as pd
import numpy as np
import os
from pathlib import Path 
from scipy.stats import iqr

def average_without_outliers(arr):
    arr = np.array(arr)
    if not np.all(arr==0):
        arr = arr[arr != 0]
        q1 = np.percentile(arr, 25)
        q3 = np.percentile(arr, 75)
        iqr_value = q3 - q1
        lower_bound = q1 - 1.5 * iqr_value
        upper_bound = q3 + 1.5 * iqr_value
        arr_no_outliers = arr[(arr >= lower_bound) & (arr <= upper_bound)]
    else:
        arr_no_outliers = arr
    return np.mean(arr_no_outliers)

# find the paths of all the csv files
f_paths = []
rtols = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
atols = [1.0e-20, 1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12, 1.0e-10]
# total 20 species
for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
    df_tols = []
    for i in range(len(rtols)):
        data_path = f'rtol_{rtols[i]}_atol_{atols[i]}/{ratio}/thermo_sens_properties.csv'
        if os.path.exists(data_path):
            df = pd.read_csv(data_path)
            rxns_not_calculated = set(range(20)) - set(df.iloc[:,1])
            result_list = sorted(list(rxns_not_calculated))
            new_df = df.iloc[:, 1:]
            if len(result_list) != 0:
                data_dict = dict()
                new_df = df.iloc[:,1:].values
                for row in range(new_df.shape[0]):
                    data_dict[new_df[row][0]] = new_df[row]     
                for i in result_list:
                    insert_zeros = np.hstack((np.array([i]), np.zeros(15)))
                    data_dict[insert_zeros[0]] = insert_zeros
                sorted_dict_ascending = dict(sorted(data_dict.items(), key=lambda item: item[0]))
                new_df = np.array(list(sorted_dict_ascending.values()))
            else:
                new_df = df.iloc[:, 1:].values
            df_tols.append(new_df)
    combined_array = np.stack((df_tols[0], df_tols[1], df_tols[2], df_tols[3], df_tols[4], df_tols[5]), axis=-1)
    processed_array = np.empty_like(df_tols[0], dtype=float)

# #     # Loop through each cell and replace values with the calculated averages
    for i in range(combined_array.shape[0]):
        for j in range(combined_array.shape[1]):
            cell_data = combined_array[i, j]
            avg = average_without_outliers(cell_data)
            processed_array[i, j] = avg
    table_names = ['Species', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
                   'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
                   'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv']
    output_dir = Path(f"thermo_sensitivity_data/{ratio}_avg_thermo.csv")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(processed_array, columns=table_names).to_csv(output_dir)
