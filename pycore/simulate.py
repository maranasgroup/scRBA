def get_GAMS_modelStat(filepath='./runRBA.modelStat.txt'):
    with open(filepath) as f:
        modelStat = f.read()
    modelStat = modelStat.replace('\n', '')
    modelStat = modelStat.replace(' ', '')
    modelStat = int(float(modelStat))
    if modelStat == 11:
        print('Licensing error')
    elif modelStat in [12,13,14]:
        stat = 'need_rerun'
    elif modelStat in [4,10,19]:
        stat = 'infeasible'
    elif modelStat == 1:
        stat = 'optimal'
    else:
        print('Feasible but not globally optimal')
        
    return stat

class RBA_result:
    def __init__(self, biom_id, growth_rate='', raw_flux='', metabolic_flux='',
                 ribo_capacity_usage=0, proteome_capacity_usage=0,
                 proteome_allocation='', proteome_allocation_by_rxns='',
                 protein_mw='', enzyme_mw='', twocol_format=False, warning=True):
        self.growth_rate = growth_rate
        self.raw_flux = raw_flux
        self.metabolic_flux = metabolic_flux
        self.ribo_capacity_usage = ribo_capacity_usage
        self.proteome_capacity_usage = proteome_capacity_usage
        self.proteome_allocation = proteome_allocation
        self.proteome_allocation_by_rxns = proteome_allocation_by_rxns
        self.protein_mw = protein_mw
        self.enzyme_mw = enzyme_mw
        self.twocol_format = twocol_format
        self.biom_id = biom_id
        self.warning = warning
        
    def load_raw_flux(self, filepath='./runRBA.flux.txt'):
        with open(filepath) as f:
            text = f.read().split('\n')
        text = [i for i in text if i != '']
        
        fluxdict = dict()
        for i in text:
            if self.twocol_format:
                r,v = i.split('\t')
            else:
                r,_,v = i.split('\t')
            fluxdict[r] = float(v)
        self.raw_flux = fluxdict
        try:
            self.growth_rate = fluxdict[self.biom_id]
        except:
            if self.warning:
                print(self.biom_id + ' is not found in raw flux. No growth rate assigned')
        
    def calculate_metabolic_flux(self):
        from utils import extract_details_from_rxnid
        metfluxdict = dict()
        for k,v in self.raw_flux.items():
            if k[:4] == 'RXN-':
                _,rxn,rxn_dir,_ = extract_details_from_rxnid(k)
                if rxn_dir == 'FWD':
                    rval = v
                elif rxn_dir == 'REV':
                    rval = -v
                
                if rxn not in metfluxdict.keys():
                    metfluxdict[rxn] = rval
                else:
                    metfluxdict[rxn] += rval
        self.metabolic_flux = metfluxdict
        
    def calculate_ribo_capacity_usage(self):
        try:
            rrna_cap = 0.8 * self.raw_flux['BIOSYN-RNA']
            try:
                rrna_unused = self.raw_flux['BIOSYN-RNA7']
            except:
                rrna_unused = 0
            self.ribo_capacity_usage = (rrna_cap - rrna_unused) / rrna_cap
        except:
            self.ribo_capacity_usage = 'N/A'
        
    def calculate_proteome_capacity_usage(self):
        try:
            dummytoprot = self.raw_flux['BIOSYN-PROTDUMMY']
            if 'BIOSYN-PROTDUMMY2' in self.raw_flux.keys():
                dummy_modeled_load = self.raw_flux['BIOSYN-PROTDUMMY2']
            else:
                dummy_modeled_load = 0
            pload = self.raw_flux['BIOSYN-PROTMODELED']
            procap = self.raw_flux['BIOSYN-PROTTOBIO']
            self.proteome_capacity_usage = (pload - dummy_modeled_load) / (procap - dummytoprot)
        except:
            self.proteome_capacity_usage = 'N/A'
        
    def calculate_proteome_allocation(self):
        if self.protein_mw == '':
            print('Cannot calculate. Need to load in protein molecular weight data ' + \
                  'as self.protein_mw = dictionary_of_protein_mw')
            return None
            
        protfluxdict = dict()
        for k,v in self.raw_flux.items():
            if k[:7] == 'PROSYN-':
                prot_id = k[7:]
                mw = self.protein_mw[prot_id]
                if '_' in prot_id:
                    prot_id = prot_id.split('_')[0]
                if prot_id not in protfluxdict.keys():
                    protfluxdict[prot_id] = v*mw
                else:
                    protfluxdict[prot_id] += v*mw
        ptot = float(self.raw_flux['BIOSYN-PROTTOBIO'])
        protfluxdict = {k:v/ptot for k,v in protfluxdict.items()}
        self.proteome_allocation = protfluxdict
        
    def calculate_proteome_allocation_by_rxns(self):
        from utils import extract_details_from_rxnid
        if self.enzyme_mw == '':
            print('Cannot calculate. Need to load in enzyme molecular weight data ' + \
                  'as self.enzyme_mw = dictionary_of_enzyme_mw')
            return None
        
        enzreq = dict()
        for k,v in self.raw_flux.items():
            if k[:8] == 'ENZLOAD-' and len(k.split('-')):
                mw = self.enzyme_mw[k]
                _,rxn,_,_ = extract_details_from_rxnid(k)
                if rxn not in enzreq.keys():
                    enzreq[rxn] = v*mw
                else:
                    enzreq[rxn] += v*mw
        ptot = float(self.raw_flux['BIOSYN-PROTTOBIO'])
        enzreq = {k:v/ptot for k,v in enzreq.items()}
        self.proteome_allocation_by_rxns = enzreq
        
    def calculate_all(self):
        self.calculate_metabolic_flux()
        self.calculate_ribo_capacity_usage()
        self.calculate_proteome_capacity_usage()
        self.calculate_proteome_allocation()
        self.calculate_proteome_allocation_by_rxns()
        
    def load_and_calculate(self, filepath='./runRBA.flux.txt'):
        self.load_raw_flux(filepath)
        self.calculate_all()
    
    def save_to_json(self, filepath):
        import json
        with open(filepath, 'w') as f:
            json.dump(self.__dict__, f, indent=4, separators=(",", ": "), sort_keys=True)
            
    def load_from_json(self, filepath):
        import json
        with open(filepath) as f:
            resdict = json.load(f)
        for k,v in resdict.items():
            self.__setattr__(k,v)
            
    def make_escher_csv(self, filepath):
        import csv
        with open(filepath, 'w') as f:
            fcsv = csv.writer(f, delimiter=',')
            fcsv.writerow(['Rxn', 'Flux'])
            for rxn,val in self.metabolic_flux.items():
                fcsv.writerow([rxn, val])
