path_gams = '../../../../GAMS/parameterization/enz_from_proteome/'
path_out = './enz_from_proteome/'
path_data = './Elsemman2022_chemoGlc030.xlsx'

run_setting_file_from = './GAMS_setting_files/enz_from_proteome_GAMS_settings.txt'
run_setting_file_to = './enz_from_proteome/enz_from_proteome_GAMS_settings.txt'
pycore_path = '../../../../pycore/'

#### Create directory and copy run settings
import os,shutil
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
shutil.copy(run_setting_file_from, run_setting_file_to);

#### Load data
import pandas as pd
df_data = pd.read_excel(path_data)
df_data.index = df_data['id'].to_list()
df_data = df_data[df_data['conc (g/gDW)'] > 0]
# Excluding ribosome protein subunit (conflicting if fit to both enzymatic and ribosomal protein data)
df_data = df_data[(df_data.type == 'truedata_enz') | (df_data.type == 'gapfill_subunit')]

#### Process data
import sys
sys.path.append(pycore_path)
from utils import metabolites_dict_from_reaction_equation_RBA

with open(os.path.join(path_gams, 'pro_and_enz.txt')) as f:
    pro_list = f.read().split('\n')
pro_list = pro_list[1:-1]
pro_list = [i[1:-1] for i in pro_list]
pro_list = [i for i in pro_list if i.split('-')[0] == 'PRO']

data = []; pro_data = []; pro_nodata = []
for met in pro_list:
    _,sid = met.split('-', maxsplit=1)
    if sid in df_data.index:
        pro_data.append("'PROIN-" + sid + "'")
        data.append("'PROIN-" + sid + "' " + str(df_data.loc[sid, 'vtrans (mmol/gDW/h)']))
    else:
        pro_nodata.append("'PROIN-" + sid + "'")
        
data = ['/'] + data + ['/']
pro_data = ['/'] + pro_data + ['/']
pro_nodata = ['/'] + pro_nodata + ['/']

# Write out run files
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data))
with open(os.path.join(path_out, 'rxns_pro_data.txt'), 'w') as f:
    f.write('\n'.join(pro_data))
with open(os.path.join(path_out, 'rxns_pro_nodata.txt'), 'w') as f:
    f.write('\n'.join(pro_nodata))
    
#### Simulation
shutil.copy(os.path.join(path_gams, 'enz_from_proteome.gms'),
            os.path.join(path_out, 'enz_from_proteome.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));
            
from simulate import get_GAMS_modelStat
from collections import OrderedDict

with open(os.path.join(path_gams, 'rxns_enz.txt')) as f:
    enz_list = f.read().split('\n')
enz_list = enz_list[1:-1]
enz_list = [i[1:-1] for i in enz_list]

enzdict = OrderedDict()
for enz in enz_list:
    #enzid = enz.split('-', maxsplit=1)[1]
    cmds = ['cd ' + path_out,
            'module load gams',
            'gams enz_from_proteome.gms --enzobj=' + enz + '  o=/dev/null']
            #'gams enz_from_proteome.gms --enzobj=' + enz]
    
    os.system('\n'.join(cmds))
    fname = os.path.join(path_out, 'enz_from_proteome.modelStat.txt')
    stat = get_GAMS_modelStat(fname)
    
    if stat == 'optimal':
        fname = os.path.join(path_out, 'enz_from_proteome.objval.txt')
        with open(fname) as f:
            v = float(f.read().split('\n')[0]) / 1e6
        enzdict[enz] = v
    else:
        enzdict[enz] = 0
        

enztext = [k+'\t'+str(v) for k,v in enzdict.items()]
with open(os.path.join(path_out, 'enz_flux_calculation.txt'), 'w') as f:
    f.write('\n'.join(enztext))
