"""
This file is to calcualte the properties based on the species ratio profile data from Cantera simulation
The properties will be used in thermo sensitivity calculation.
Use it only after the species ratio profile data has been generated.
"""


import pandas as pd
import numpy as np
import os

def deriv(gas_out):
    deriv = []
    for x in range(len(gas_out) - 1):
        deriv.append((gas_out[x+1] - gas_out[x])/.01)
    deriv.append(0.)
    return deriv


def calculate(data, type='sens'):
    """
    Calculate properties of interest from the raw data
    :param data: the data
    :param type: 'sens' for sensitivity analyses
                 'output' for saving the output csv
                 'ratio' for plotting
    :return:
    """
    gas_out_data, gas_names_data, dist_array_data, T_array_data, n_surf_reactions = data

    reference = []
    for a in range(len(gas_names_data)):
        reference.append([gas_names_data[a], [gas_out_data[:, a]]])

    # This is the index to choose 
#     sens_T = (np.array(T_array_data).max() - T_array_data[1000]) / 2 + T_array_data[1000]
#     difference_array = np.absolute(T_array_data-sens_T)
#     sens_id = difference_array.argmin()
#     print(sens_id)
    
#     if sens_id <= 1000 or sens_id >= 2000:
#         sens_id = 1010
    sens_id = 1045

    for x in reference:
        if x[0] == 'CH4(2)':
            ch4_in = x[1][0][0]
            ch4_out = x[1][0][sens_id]
            if ch4_out < 0:
                ch4_out = 0.
            ch4_depletion = ch4_in - ch4_out
            reference_ch4_conv = ch4_depletion / ch4_in  # Sensitivity definition 7: CH4 conversion

            d_ch4 = deriv(x[1][0])
            reference_max_ch4_conv = min(d_ch4)  # Sensitivity definition 15: maximum rate of CH4 conversion

            conv50_pos = 0
            for y in range(len(x[1][0])):
                if (ch4_in - x[1][0][y]) / ch4_in >= 0.5:
                    if conv50_pos == 0:
                        conv50_pos = y
                        reference_dist_to_50_ch4_conv = dist_array_data[conv50_pos] # Sensitivity definition 14: distance to 50% CH4 conversion
                else:
                    # never reached 50% conversion
                    reference_dist_to_50_ch4_conv = 510.
        if x[0] == 'Ar':
            ar = x[1][0][sens_id]
        if x[0] == 'O2(3)':
            o2_in = x[1][0][0]
            o2_out = x[1][0][sens_id]
            if o2_out < 0:
                o2_out = 0.  # O2 can't be negative
            elif o2_out > o2_in:
                o2_out = o2_in  # O2 can't be created, to make it equal to O2 in
            o2_depletion = o2_in - o2_out
            reference_o2_conv = o2_depletion / o2_in  # Sensitivity definition 13: O2 conversion
        if x[0] == 'CO(7)':
            co_out = x[1][0][sens_id]
        if x[0] == 'H2(6)':
            h2_out = x[1][0][sens_id]
        if x[0] == 'H2O(5)':
            h2o_out = x[1][0][sens_id]
        if x[0] == 'CO2(4)':
            co2_out = x[1][0][sens_id]

    ratio = ch4_in / (2 * o2_in)

    # negative sensitivity is higher selectivity
    reference_h2_sel = h2_out / (ch4_depletion * 2)  # Sensitivity definition 5: H2 selectivity
    if reference_h2_sel <= 0:
        reference_h2_sel = 1.0e-15  # selectivity can't be 0

    reference_co_sel = co_out / ch4_depletion  # Sensitivity definition 3: CO selectivity
    if reference_co_sel <= 0:
        reference_co_sel = 1.0e-15  # selectivity can't be 0

    reference_syngas_selectivity = reference_co_sel + reference_h2_sel  # Sensitivity definition 1: SYNGAS selectivity

    reference_syngas_yield = reference_syngas_selectivity * reference_ch4_conv  # Sensitivity definition 2: SYNGAS yield
    if reference_syngas_yield <= 0:
        reference_syngas_yield = 1.0e-15  # yield can't be 0

    reference_co_yield = co_out / ch4_in  # Sensitivity definition 4: CO % yield
    # reference_co_yield = reference_co_sel * reference_ch4_conv

    reference_h2_yield = h2_out / (2 * ch4_in)  # Sensitivity definition 6: H2 % yield
    # reference_h2_yield = reference_h2_sel * reference_ch4_conv

    # Sensitivity definition 8: H2O + CO2 selectivity
    reference_h2o_sel = h2o_out / (ch4_depletion * 2)
    reference_co2_sel = co2_out / ch4_depletion
    if reference_h2o_sel <= 0:
        reference_h2o_sel = 1.0e-15  # H2O selectivity can't be 0
    if reference_co2_sel <= 0:
        reference_co2_sel = 1.0e-15  # CO2 selectivity can't be 0
    reference_full_oxidation_selectivity = reference_h2o_sel + reference_co2_sel

    # Sensitivity definition 9: H2O + CO2 yield
    reference_full_oxidation_yield = reference_full_oxidation_selectivity * reference_ch4_conv

    # Sensitivity definition 10: exit temperature
    reference_exit_temp = T_array_data[-1]

    # Sensitivity definition 11: peak temperature
    reference_peak_temp = max(T_array_data)

    # Sensitivity definition 12: distance to peak temperautre
    reference_peak_temp_dist = dist_array_data[T_array_data.index(max(T_array_data))]
    
    sens_property = [reference_syngas_selectivity, reference_syngas_yield, reference_co_sel, reference_co_yield, reference_h2_sel, 
                     reference_h2_yield, reference_ch4_conv, reference_full_oxidation_selectivity, reference_full_oxidation_yield, 
                     reference_exit_temp, reference_peak_temp, reference_peak_temp_dist, reference_o2_conv, reference_max_ch4_conv, 
                     reference_dist_to_50_ch4_conv]
    
    ref_data = [ratio, ch4_in, ch4_out, co_out, h2_out, h2o_out, co2_out, reference_exit_temp, reference_peak_temp, 
                reference_peak_temp_dist, reference_o2_conv, reference_max_ch4_conv, reference_dist_to_50_ch4_conv]
    if type is 'sens':
        return sens_property
    elif type is 'ratio':
        return reference_co_sel, reference_h2_sel, reference_ch4_conv, reference_exit_temp, reference_o2_conv, reference_co2_sel, reference_h2o_sel
    elif type is 'gas_data':
        return ratio, reference
    else:
        return ref_data

if __name__ == "__main__":        
    # find the paths of all the csv files
    f_paths = []
    rtols = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
    atols = [1.0e-20, 1.0e-18, 1.0e-16, 1.0e-14, 1.0e-12, 1.0e-10]
                
    for i in range(len(rtols)):
        for ratio in [.6, 1., 1.1, 1.2, 1.6, 2., 2.6]:
            f_paths.append(f'therm_sens/rtol_{rtols[i]}_atol_{atols[i]}/{ratio}')
            
    for f in f_paths:
        data_to_csv = []
        if os.path.exists(f):
            for sp_id in range(20):
                file_path = os.path.join(f, f'therm_sens_{sp_id}.csv')
                if os.path.exists(file_path):
                    output = pd.read_csv(file_path)
                    if len(output) >= 7000:
                        output_data = []
                        gas_data = []
                        gas_out = output.iloc[:,3:25]
                        for i in range(len(gas_out)):
                            gas_data.append(np.array(gas_out.iloc[i, :]))
                        gas_data = np.array(gas_data)
                        gas_names = list(output.columns[3:25])
                        dist_arr = list(output.iloc[:,1])
                        T_arr = list(output.iloc[:,2])
                        output_data.append(gas_data)
                        output_data.append(gas_names)
                        output_data.append(dist_arr)
                        output_data.append(T_arr)
                        output_data.append(107)
                        sens_data = calculate(output_data, type='sens')
                        sens_data.insert(0, sp_id)
                        data_to_csv.append(sens_data)
                    else:
                        sens_data = [0] * 15
                        sens_data.insert(0, sp_id)
                        data_to_csv.append(sens_data)
            table_names = ['Species', 'SYNGAS Selec', 'SYNGAS Yield', 'CO Selectivity', 'CO % Yield', 'H2 Selectivity', 
                           'H2 % Yield', 'CH4 Conversion', 'H2O+CO2 Selectivity', 'H2O+CO2 yield', 'Exit Temp', 'Peak Temp',
                           'Dist to peak temp', 'O2 Conversion', 'Max CH4 Conv', 'Dist to 50 CH4 Conv']
            df = pd.DataFrame(np.array(data_to_csv), columns=table_names)
            df.to_csv(os.path.join(f, 'thermo_sens_properties.csv'))
