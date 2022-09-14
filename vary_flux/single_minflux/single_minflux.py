from vary_flux_options import *
gams_name = 'vary_flux'

import sys
sys.path.append(path_pycore)
import json
from simulate import get_GAMS_modelStat, RBA_result

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

# Set reaction to be vary
rxn = 'DUTPDP_c'

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
minflux_path = './minflux.txt'
with open(minflux_path, 'w') as f: f.write('')

# Prepare optimization
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
#### Minimizing
sumparts = []
for p in fwds:
    sumparts.append("(-1) * " + "v('" + p + "')")
for p in revs:
    sumparts.append("(1) * " + "v('" + p + "')")

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
    if objval != 0:
        objval = -objval
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
            if objval != 0:
                objval = -objval
            
    if stat == 'need_rerun':
        objval = 'N/A'
        text = 'RBA optimization stills need rerun after ' + str(itermax) + ' iterations for ' + rxn + '\n'
        print(text)
        with open(runlog_path, 'a') as f:
            f.write(text)
            
    # Reset vglc to original value from the slightly perturbed value in "need_rerun" rerun
    vglc = vglc0
        
# Write to result file
with open(minflux_path, 'a') as f:
    text = [rxn, stat, str(objval)]
    f.write('\t'.join(text) + '\n')
####

# Write result
with open(path_enz_mw) as f:
    text = f.read().split('\n')
text = [i for i in text if i != '']
enz_mw = dict()
for i in text:
    k,v = i.split('\t')
    enz_mw[k] = float(v)
    
with open(path_pro_mw) as f:
    text = f.read().split('\n')
text = [i for i in text if i != '']
pro_mw = dict()
for i in text:
    k,v = i.split('\t')
    pro_mw[k] = float(v)
    
res = RBA_result(biom_id=biom_id)
res.load_raw_flux('./runRBA.flux.txt')
res.enzyme_mw = enz_mw
res.protein_mw = pro_mw

res.calculate_all()
res.enzyme_mw = ''
res.protein_mw = ''

res.save_to_json('./RBA_result.json')
res.make_escher_csv('./flux.escher.csv')
