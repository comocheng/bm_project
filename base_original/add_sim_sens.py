import os
import re
import pandas as pd
import numpy as np
import shutil
import subprocess
import multiprocessing

with open(os.path.join('/work/westgroup/chao/sketches/cpox_sim/bm_models_final', 'simulation_sens.py')) as infile:
    input_file = infile.read()

base_directory = 'binding_energies'
def directory(carbon,oxygen,index):
    return os.path.join(base_directory, "{}_c{:.2f}o{:.2f}".format(index,carbon,oxygen))

def make_input(binding_energies):
    """
    Make an input file for the given (carbon,oxygen) tuple (or iterable) of binding energies
    and return the name of the directory in which it is saved.
    """
    carbon, oxygen, index = binding_energies
    output = input_file
    out_dir = directory(carbon, oxygen, index)
    new_dir = f"out_root = '/work/westgroup/chao/sketches/cpox_sim/bm_models_final/base_original/{out_dir}'"
    output = re.sub("out_root = \'/home/xu.chao/sketches/cpox_sim/rmg_models/base_cathub\'", new_dir, output)
    os.path.exists(out_dir) or os.makedirs(out_dir)
    out_file = os.path.join(out_dir, 'simulation_sens.py')
    with open(out_file,'w') as outfile:
        outfile.write(output)    
    return out_dir

carbon_range = (-7.5, -5.5)
oxygen_range = (-5.25, -3.25)
grid_size = 9
mesh  = np.mgrid[carbon_range[0]:carbon_range[1]:grid_size*1j, oxygen_range[0]:oxygen_range[1]:grid_size*1j]

experiments = mesh.reshape((2,-1)).T
input = []
idx = 80
while idx>=0:
    new_arr = np.append(experiments[idx], abs(idx - 81))
    input.append(new_arr)
    idx -= 1
list(map(make_input, input))