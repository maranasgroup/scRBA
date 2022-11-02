from vary_flux_options import *
gams_name = 'vary_flux'

import sys
sys.path.append(path_pycore)
import json
from simulate import get_GAMS_modelStat

def get_objval_from_txt(filepath='./objval.txt'):
    with open(filepath) as f:
        objval = f.read()
    objval = objval.replace('\n', '')
    objval = objval.replace(' ', '')
    objval = float(objval)
    return objval

import shutil
shutil.copy(path_gams + 'vary_flux/' + gams_name + '.gms', './' + gams_name + '.gms');
shutil.copy(path_gams + 'vary_flux/soplex.opt', './soplex.opt');

import os
import pandas as pd
import numpy as np

# Set glucose uptake rate here
vglc0 = 13.2
vglc = vglc0

# Load details on metabolic reactions' components in RBA model
df_rxninfo = pd.read_excel('../input/rxn_details.xlsx')
df_rxninfo.index = df_rxninfo.rxn.to_list()

# Initiate results storing dataframe
# df_rxninfo['lb'] = [np.nan] * df_rxninfo.shape[0]
# df_rxninfo['ub'] = [np.nan] * df_rxninfo.shape[0]

# Initiate run log
runlog_path = './runlog.txt'
with open(runlog_path, 'w') as f: f.write('')

# Initiate result storage txt files
maxflux_path = './maxflux.txt'
with open(maxflux_path, 'w') as f: f.write('')

# Set start index (of reaction list), will run to the end
istart = 0
iend = df_rxninfo.shape[0]
rxns_all = df_rxninfo.rxn.to_list()

for i in range(istart,iend):
    rxn = rxns_all[i]

    # Load forward and reverse components of a metabolic flux
    fwds = df_rxninfo.rxn_comp_fwd[rxn]
    if pd.isnull(fwds) == False:
        fwds = fwds.split(',')
    else:
        fwds = []

    revs = df_rxninfo.rxn_comp_rev[rxn]
    if pd.isnull(revs) == False:
        revs = revs.split(',')
    else:
        revs = []
    
    # Write objective mathematical expression in GAMS
    #### Maximizing
    sumparts = []
    for p in fwds:
        sumparts.append("(1) * " + "v('" + p + "')")
    for p in revs:
        sumparts.append("(-1) * " + "v('" + p + "')")

    obj = 'Obj.. z =e= ' + ' + '.join(sumparts) + ';'
    with open('./obj.txt', 'w') as f:
        f.write(obj)
    
    # Execute GAMS and get result
    os.system('module load gams\n' + 'gams ' + gams_name + '.gms' + \
              ' --vglc=' + str(vglc) + ' o=/dev/null')
              #' --vglc=' + str(vglc))
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
    if stat == 'infeasible':
        objval = 'N/A'
        text = 'RBA optimization is infeasible for ' + rxn + '\n'
        print(text)
        with open(runlog_path, 'a') as f:
            f.write(text)
            
    elif stat == 'optimal':
        objval = get_objval_from_txt()
    elif stat == 'need_rerun':
        itermax = 100; iternum = 0;
        while stat == 'need_rerun' and iternum < itermax:
            iternum += 1
            vglc += 1e-3
            os.system('module load gams\n' + 'gams ' + gams_name + '.gms' + \
                      ' --vglc=' + str(vglc) + ' o=/dev/null')
                      #' --vglc=' + str(vglc))
            stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
            if stat == 'infeasible':
                objval = 'N/A'
                text = 'RBA optimization is infeasible for ' + rxn + '\n'
                print(text)
                with open(runlog_path, 'a') as f:
                    f.write(text)
                
            elif stat == 'optimal':
                objval = get_objval_from_txt()
                
        if stat == 'need_rerun':
            objval = 'N/A'
            text = 'RBA optimization stills need rerun after ' + str(itermax) + ' iterations for ' + rxn + '\n'
            print(text)
            with open(runlog_path, 'a') as f:
                f.write(text)
                
        # Reset vglc to original value from the slightly perturbed value in "need_rerun" rerun
        vglc = vglc0
            
    # Write to result file
    with open(maxflux_path, 'a') as f:
        text = [str(i), rxn, stat, str(objval)]
        f.write('\t'.join(text) + '\n')
    ####
