'''
Created on 2010-11-11

@author: Madlee
'''

from kuai.unit import default_units
from kuai.mol import Molecule

class Model:
    pass

def setup_model(mols, partable, par):
    if isinstance(mols, Molecule):
        mols = [(mols, 1)]
    elif len(mols) == 2 and isinstance(mols[0], Molecule):
        mols = [mols]
    result = Model()
    result.parameters = par
    result.mols = mols
    result.charges = []
    result.masses = []
    result.coords = []
    
    result.atomtype = []
    result.group_id = []
    
    masses_table = partable['atom/mass']
    type2id = partable['atom/type']
    
    amu, _, _ = default_units.parse(1, 'amu')
    angerstron, _, _ = default_units.parse(1, 'A')
    
    next_group_id = 0
    
    for m, n in mols:
        charges = []
        masses = []
        coords = []
        atomtype = []
        group_id = []
        for a in m.atoms:
            if hasattr(a, 'partial_charge'):
                charges.append(a.partial_charge)
            else:
                charges.append(0.0)
            if a.type in masses_table:
                masses.append(par[masses_table[a.type]])
            else:
                masses.append(a.weight * amu)
            if hasattr(a, 'group_id'):
                group_id.append(a.group_id)
                
            assert a.type in type2id
            atomtype.append(type2id[a.type])
            coords.append(a.coord.x * angerstron)
            coords.append(a.coord.y * angerstron)
            coords.append(a.coord.z * angerstron)
        
        result.masses += masses * n
        result.coords += coords * n
        result.charges += charges * n
        result.atomtype += atomtype * n
        
        group_id_map = {}
        for i in group_id:
            group_id_map[i] = 0
        gid_in_m = [i for i in group_id_map.iterkeys()]
        gid_in_m.sort()
        for i in range(len(gid_in_m)):
            group_id_map[gid_in_m[i]] = i
            
        for i in range(len(group_id)):
            group_id[i] = group_id_map[group_id[i]]
        
        for i in range(n):
            result.group_id += [gid + next_group_id for gid in group_id]
            next_group_id += len(gid_in_m)
            
    
    if len(result.charges) != len(result.masses):
        del result.charges
    if len(result.group_id) != len(result.masses):
        del result.group_id
    
    return result
