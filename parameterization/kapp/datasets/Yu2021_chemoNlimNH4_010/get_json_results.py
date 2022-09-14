path_pycore = '../../../../pycore/'
path_enz_mw = '../../../../input/enz_mw_g_per_mmol.txt'
path_pro_mw = '../../../../input/pro_mw_g_per_mmol.txt'
biom_id = 'BIOSYN-BIODILSTARVE' # Set the ID of the biomass reaction

import sys
sys.path.append(path_pycore)
import json
from simulate import get_GAMS_modelStat, RBA_result

final_res = RBA_result(biom_id=biom_id);
final_res.load_raw_flux(filepath='./min_flux_sum/min_flux_sum.flux.txt')

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
