from runRBA_options import *

import sys
sys.path.append(path_pycore)
import json
from simulate import get_GAMS_modelStat, RBA_result

import shutil
shutil.copy(path_gams + 'application/runRBA_max_prod.gms', './runRBA_max_prod.gms');
shutil.copy(path_gams + 'application/soplex.opt', './soplex.opt');

import os
import pandas as pd

# Set growth and glucose uptake rates
mu = 0.1
vglc0 = 13.2

# Initiate report
report = {k:None for k in ['stat', 'vglc', 'vprod', 'yield']}

# Execute GAMS
vglc = vglc0
os.system('module load gams\n' + 'gams runRBA_max_prod.gms' + \
          ' --mu=' + str(mu) + \
          ' --vglc=' + str(vglc) + \
          ' --vprod="' + vprod + '"' + ' o=/dev/null')
    #     ' --vprod="' + vprod + '"')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
if stat == 'infeasible':
    report['stat'] = stat
    report['vglc'] = 0
    report['vprod'] = 0
    report['yield'] = 0
        
elif stat == 'optimal':
    res = RBA_result(biom_id=biom_id)
    res.load_raw_flux(filepath='./runRBA.flux.txt')
    res.calculate_metabolic_flux()
    if vprod.split('-')[0] == 'RXNADD':
        res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod]
    pflux = res.metabolic_flux[vprod_coreid]
    vglc_sim = - res.metabolic_flux['EX_glc__D_e']
        
    report['stat'] = stat
    report['vglc'] = vglc_sim
    report['vprod'] = pflux
    report['yield'] = pflux * prod_mw / vglc_sim / 180.156
        
elif stat == 'need_rerun':
    itermax = 100; iternum = 0;
    while stat == 'need_rerun' and iternum < itermax:
        iternum += 1
        vglc += 1e-3
        os.system('module load gams\n' + 'gams runRBA_max_prod.gms' + \
                  ' --mu=' + str(mu) + \
                  ' --vglc=' + str(vglc) + \
                  ' --vprod="' + vprod + '"' + ' o=/dev/null')
        stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
            
        if stat == 'infeasible':
            report['stat'] = stat
            report['vglc'] = 0
            report['vprod'] = 0
            report['yield'] = 0
            
        elif stat == 'optimal':
            res = RBA_result(biom_id=biom_id)
            res.load_raw_flux(filepath='./runRBA.flux.txt')
            res.calculate_metabolic_flux()
            if vprod.split('-')[0] == 'RXNADD':
                res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod]
            pflux = res.metabolic_flux[vprod_coreid]
            vglc_sim = - res.metabolic_flux['EX_glc__D_e']
            
            report['stat'] = stat
            report['vglc'] = vglc_sim
            report['vprod'] = pflux
            report['yield'] = pflux * prod_mw / vglc_sim / 180.156
        
        elif stat == 'need_rerun':
            report['stat'] = stat
        
        else:
            print('wtf')
        
# Write report text file
text = []
for k in ['stat', 'vglc', 'vprod', 'yield']:
    text.append(k + '\t' + str(report[k]))
with open('report.txt', 'w') as f:
    f.write('\n'.join(text))

# Write JSON results
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
    
res.enzyme_mw = enz_mw
res.protein_mw = pro_mw

res.calculate_all()
res.enzyme_mw = ''
res.protein_mw = ''

res.save_to_json('./RBA_result.json')
res.make_escher_csv('./flux.escher.csv')
