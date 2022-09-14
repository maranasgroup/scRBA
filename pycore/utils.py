def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def metabolites_dict_from_reaction_equation_RBA(eqn, split=False):
    import re
    rStr, pStr = re.split('-->|->|<--|<-|<=>|<->', eqn)
    rStr = rStr.strip(' ')
    pStr = pStr.strip(' ')
    rs = re.split(' \+ | \+|\+ |\+', rStr)
    ps = re.split(' \+ | \+|\+ |\+', pStr)

    r_dict = dict()
    for r in rs:
        if ' ' in r:
            val, met = re.split('\s+', r)
            if is_number(val):
                r_dict[met] = -float(val)
            else:
                r_dict[met] = '-' + val
        else:
            r_dict[r] = -1.0

    p_dict = dict()
    for p in ps:
        if ' ' in p:
            val, met = re.split('\s+', p)
            if is_number(val):
                p_dict[met] = float(val)
            else:
                p_dict[met] = val
        else:
            p_dict[p] = 1.0
            
    if split:
        return r_dict, p_dict
    else:
        return merge_two_dicts(r_dict, p_dict)

def build_reaction_equation_from_metabolites_dict_RBA(met_dict, arrow='<=>', floatdecimal=6):
    lhs = []; rhs = [];
    for k,v in met_dict.items():
        if is_number(v):
            v = float(v)
            if v == -1:
                lhs.append(k)
            elif v == 1:
                rhs.append(k)
            elif v < 0 and v != -1 and v.is_integer():
                lhs.append(' '.join([str(-int(v)), k]))
            elif v > 0 and v != 1 and v.is_integer():
                rhs.append(' '.join([str(int(v)), k]))
            elif v < 0 and v != -1:
                lhs.append(' '.join([('{:.' + str(floatdecimal) + 'f}').format(-v), k]))
            elif v > 0 and v != 1:
                rhs.append(' '.join([('{:.' + str(floatdecimal) + 'f}').format(v), k]))
        else:
            if v[0] == '-':
                lhs.append(v[1:] + ' ' + k)
            else:
                rhs.append(v + ' ' + k)
                
    return ' '.join([ ' + '.join(lhs), arrow, ' + '.join(rhs)])

def extract_details_from_rxnid(rxn_id):
    idsplit = rxn_id.split('-')
    tag = idsplit[0]
    rxn_base_id = idsplit[1]
    enz_id = rxn_id[(len(tag) + len(rxn_base_id) + 2):]
    rxn_dir = rxn_base_id.split('_')[-1]
    rxn_base_id = rxn_base_id[:-len(rxn_dir)-1]
    return(tag,rxn_base_id,rxn_dir,enz_id)
