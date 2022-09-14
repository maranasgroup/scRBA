path_gams = '../../../../GAMS/parameterization/min_flux_sum/'
path_B2 = './min_flux_violation/'
path_out = './min_flux_sum/'

run_setting_file_from = './GAMS_setting_files/min_flux_sum_GAMS_settings.txt'
run_setting_file_to = './min_flux_sum/min_flux_sum_GAMS_settings.txt'

#### Create directory and copy run settings
import os,shutil
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
shutil.copy(run_setting_file_from, run_setting_file_to);

#### Load proteomics data and write protein translation fluxes
# Import from min_flux_violation outputs in GAMS

#### Append essential inactive reactions found in B2_min_flux_violation.py
fname = os.path.join(path_B2, 'rxns_inactive.txt')
with open(fname) as f:
    rxns_inactive_B2 = f.read().split('\n')[1:-1]
rxns_inactive_B2 = [i[1:-1] for i in rxns_inactive_B2]

fname = os.path.join(path_B2, 'min_flux_violation.flux_essential_inactive_rxns.txt')
with open(fname) as f:
    rxns_essential_inactive_B2 = f.read().split('\n')
rxns_essential_inactive_B2 = [i.split('\t')[0] for i in rxns_essential_inactive_B2 if i != '']

rxns_inactive = [i for i in rxns_inactive_B2 if i not in rxns_essential_inactive_B2]
rxns_inactive = ["'" + i + "'" for i in rxns_inactive]
rxns_inactive = ['/'] + rxns_inactive + ['/']
fname = os.path.join(path_out, 'rxns_inactive.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(rxns_inactive))

#### Simulation
import shutil
shutil.copy(os.path.join(path_gams, 'min_flux_sum.gms'),
            os.path.join(path_out, 'min_flux_sum.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));
            
cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_sum.gms  o=/dev/null']
        #'gams min_flux_sum.gms']
os.system('\n'.join(cmds))

#### Convert GAMS-scaled flux to actual flux
# All fluxes
fname = os.path.join(path_out, 'min_flux_sum.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
fname = os.path.join(path_out, 'min_flux_sum.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
