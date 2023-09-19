import cantera as ct
import os
import shutil
import pandas as pd
import numpy as np

# copy the reactions from original model to BM model
def read_files(or_p, bm_p):
    gas_bm = ct.Solution(bm_p, 'gas')
    surf_bm = ct.Interface(bm_p, 'surface1', [gas_bm])
    rxn_num = surf_bm.n_reactions
    for i in range(rxn_num):
        gas_or = ct.Solution(or_p, 'gas')
        surf_or = ct.Interface(or_p, 'surface1', [gas_or])
        surf_or.reaction(i).rate = surf_bm.reaction(i).rate
        os.path.exists(f'bm2or') or os.mkdir('bm2or')
        # print(1)
        os.path.exists(f'bm2or/rxn{i}') or os.mkdir(f'bm2or/rxn{i}')
        surf_or.write_yaml(f'bm2or/rxn{i}/cantera.yaml')
        shutil.copy('simulation.py', f'bm2or/rxn{i}/simulation.py')

or_p = os.path.join(os.getcwd(), 'cantera.yaml')
bm_p = or_p.replace('base_original', 'base_bm')
read_files(or_p, bm_p)
