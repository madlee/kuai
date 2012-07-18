'''
Created on Jun 8, 2012

@author: fli
'''

from kuai.mol import Bond, Molecule
from kuai.mol.rings import RingSet

MAX_AROMATIC_RING = 10

AROMATIC_CANDIDATES = {'C', 'N', 'O', 'P', 'S', 'As', 'Se'}

class AromaticSet(object):
    """Contain all aromatic blocks of a molecule."""
    def __init__(self, mol, rings=None, hints=None):
        if hints:
            atoms = [i for i in mol.atoms if hints.is_aromatic(i)]
        else:
            atoms = [i for i in mol.atoms if i.symbol in AROMATIC_CANDIDATES]
        
        candidates = self.__get_aromatic_candidates(mol, atoms)
        for block, rings in candidates.iteritems():
            pass
        
    
    def __get_aromatic_candidates(self, mol, atoms):
        candidates = set(atoms)
         
        while candidates:
            discard = set()
            for atom_i in candidates :
                deg = mol.degree(atom_i)
                total_order = 0
                n_double = 0 
                for j in range(deg):
                    atom_j = mol.neighbor_atom(atom_i, j)
                    if atom_j in candidates:
                        bond_j = mol.neighbor_bond(atom_i, j)
                        if bond_j.order > Bond.DOUBLE_BOND:
                            discard.add(atom_i)
                            discard.add(atom_j)
                            total_order = 0
                            break
                            
                        total_order += bond_j.order
                        if bond_j.order > Bond.SINGLE_BOND:
                            n_double += 1
                
            if not discard:
                break
            else:
                candidates = candidates.difference(discard)
        
        mol = Molecule(atoms, None, mol)
        return RingSet(mol, MAX_AROMATIC_RING).blocks

