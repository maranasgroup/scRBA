import time
time0 = time.time()

from binary_search_options import *

if mu_min0 < 0:
    print('Invalid negative value is set as the lower bound of mu. Mu is always positive. Revert to zero.')
    mu_min0 = 0

import sys
sys.path.append(path_pycore)
import json
from simulate import get_GAMS_modelStat, RBA_result

import shutil
shutil.copy(path_gams + 'runRBA.gms', './runRBA.gms');
shutil.copy(path_gams + 'soplex.opt', './soplex.opt');

import os
# Test evaluation at zero
os.system('module load gams\n' + 'gams runRBA.gms --mu=0 o=/dev/null')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
if stat == 'infeasible':
    text = 'Model is infeasible at mu = 0, check model connectivity and other bugs. Terminate program.'
    print(text)
    with open(report_file, 'w') as f:
        f.write(text + '\n')
    quit()
    
elif stat == 'optimal':
    with open(report_file, 'w') as f:
        f.write('Pass mu = 0 test, now proceed to binary search.\n')
    
mu_min = mu_min0; mu_max = mu_max0; itercount = 0;

# Start binary search
class bcolors:
    GREEN = '\033[92m' #GREEN
    RED = '\033[91m' #RED
    RESET = '\033[0m' #RESET COLOR
    
# Evaluate min feasibility
mu = mu_min;
os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')

if stat == 'need_rerun':
    while stat == 'need_rerun':
        print(f"mu = {mu:.7f}, status = {stat}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
        mu += mu_tol
        os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
        stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
        
    if stat == 'optimal':
        mu_min = mu
        print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
            
if stat == 'infeasible':
    mu_max = mu_min
    mu_min = 0
    print(f"{bcolors.RED}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
    with open(report_file, 'a') as f:
        f.write(f"mu = {mu:.7f}, status = {stat}\n")

# Evaluate max infeasibility
mu = mu_max;
os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
while stat == 'optimal':
    print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}\n")
    with open(report_file, 'a') as f:
        f.write(f"mu = {mu:.7f}, status = {stat}\n")
    mu_max = 1.5*mu_max
    mu = mu_max;
    os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
# Update min-max
mu = float(mu_min + mu_max) / 2; final_res = RBA_result(biom_id=biom_id);
while mu_max - mu_min > mu_tol and itercount < maxiter:
    itercount += 1
        
    os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
    if stat == 'need_rerun':
        while stat == 'need_rerun':
            print(f"mu = {mu:.7f}, status = {stat}")
            with open(report_file, 'a') as f:
                f.write(f"mu = {mu:.7f}, status = {stat}\n")
            mu += mu_tol
            os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
            stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
    if stat == 'optimal':
        mu_min = mu
        final_res.load_raw_flux(filepath='./runRBA.flux.txt')
        print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
        
    elif stat == 'infeasible':
        mu_max = mu
        print(f"{bcolors.RED}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
            
    else:
        print("Error unknown!")
        quit()
        
    mu = float(mu_min + mu_max) / 2
    
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
    
final_res.enzyme_mw = enz_mw
final_res.protein_mw = pro_mw

final_res.calculate_all()
final_res.enzyme_mw = ''
final_res.protein_mw = ''

final_res.save_to_json('./RBA_result.json')
final_res.make_escher_csv('./flux.escher.csv')

dt = float(time.time() - time0) / 60
with open(report_file, 'a') as f:
    f.write(f"Run time = {dt:.2f} mins")
