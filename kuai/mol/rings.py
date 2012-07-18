'''
Created on Jun 7, 2012

@author: fli
'''

from kuai.mol import Molecule
from kuai.mol.algo import MoleculeVisitor, bft, split

MAX_RING_SIZE = 1024

class RingFinder(MoleculeVisitor):
    def __init__(self, bond, max_ring_size):
        self.bond = bond
        self.trace = {bond.atom2:(None, 2)}
        self.max_ring_size = max_ring_size
    
    def start(self, mol, atom):
        if self.trace[atom][1] > self.max_ring_size:
            # Too big. abandoned.
            raise self
        
    def find(self, mol, target, source):
        assert target not in self.trace
        assert source in self.trace
        self.trace[target] = (source, self.trace[source][1]+1)
    
    def back(self, mol, target, source):
        if target == self.bond.atom1 and source != self.bond.atom2:
            self.trace[target] = (source, self.trace[source][1]+1)
            raise self
    
    def get_result(self):
        if self.bond.atom1 in self.trace:
            result = []
            o = self.bond.atom1
            while o != self.bond.atom2:
                result.append(o)
                o = self.trace[o][0]

            result.append(self.bond.atom2)
            return result
        else:
            return None
            
    @staticmethod
    def find_ring(mol, bond, max_ring_size):
        finder = RingFinder(bond, max_ring_size)
        try:
            flags = dict([(i, MoleculeVisitor.FLAG_WHITE,) for i in mol.atoms])
            flags[bond.atom1] = MoleculeVisitor.FLAG_BLACK
            bft(mol, finder, flags, [bond.atom2])
            return None
        except RingFinder:
            return finder.get_result()

class RingSet(object):
    """Contain all rings' info of a molecule."""
    def __init__(self, mol, max_ring_size=MAX_RING_SIZE, hints=None):
        self.max_ring_size = max_ring_size
        self.blocks = {}
        self.rings = []
        candidates = self.__get_ring_candidates(mol, hints)
        
        for i in candidates:
            if len(i.atoms) == len(i.bonds):
                self.blocks[i] = [i]
                self.rings.append(i)
            else:
                assert len(i.atoms) < len(i.bonds)
                self.__make_sssr(i, max_ring_size)
    
    
    def rank(self, v):
        """Return the smallest ring rank of the atom/bond V. Return None if v is not on ring"""
        if not hasattr(self, '__rank'):
            self.__rank = {}
            for i in self.rings:
                n = len(i.atoms)
                for atom in i.atoms:
                    if atom not in self.__rank[atom] or self.__rank[atom] > n:
                        self.__rank[atom] = n
                for bond in i.bonds:
                    if bond not in self.__rank[bond] or self.__rank[bond] > n:
                        self.__rank[atom] = n
        if v in self.__rank[v]:
            return self.__rank[v]
        else:
            return None
        
    def spiro_atoms(self):
        if not hasattr(self, '__spiro_atoms'):
            self.__spiro_atoms = []
            for block in self.blocks.iterkeys():
                for atom in block.atoms:
                    if block.degree(atom) >= 4:
                        self.__spiro_atoms.append(atom)
        return self.__spiro_atoms
    
    def fused_bonds(self):
        if not hasattr(self, '__fused_bonds'):
            self.__fused_bonds = []
            for block, rings in self.blocks.iteritems():
                for bond in block.bonds:
                    if block.degree(bond.atom1) >= 3 and block.degree(bond.atom2) >= 3:
                        n = 0
                        for ring in rings:
                            if bond in ring:
                                n += 1
                        if n >= 2:
                            self.__fused_bonds.append(bond)
        return self.__fused_bonds

    def __get_ring_candidates(self, mol, hints):
        flags = dict([(atom, 1) for atom in mol.atoms])
        if hints:
            for atom in mol.atoms:
                rank = hints.rank(atom)
                if not rank or rank > self.max_ring_size:
                    flags[atom] = 0                    
                
        while True:
            no_change = True
            for atom in mol.atoms:
                if flags[atom]:
                    deg = mol.degree(atom)
                    if deg > 1:
                        v = 0
                        for i in xrange(deg):
                            v += flags[mol.neighbor_atom(atom, i)]
                            if v > 1:
                                break
                        if v < 2:
                            flags[atom] = 0
                            no_change = False
                    else:
                        flags[atom] = 0
                        no_change =  False
            
            if no_change:
                break
        
        atoms = [k for k, v in flags if v]
        return split(Molecule(atoms, None, mol))
    
    def __make_sssr(self, part, max_ring_size):
        assert len(part.atoms) < len(part.bonds)
        rings = []

        visited_bonds = set()
        for bond in part.bonds:
            if bond not in visited_bonds: 
                ring = RingFinder.find_ring(part, bond, max_ring_size)
                if ring:
                    rings.append(ring)
                    for i in ring.bond:
                        visited_bonds.add(i)
        
        assert len(rings) == len(part.bonds) - len(part.atoms) + 1
        
        self.rings += rings
        if len(visited_bonds) == len(part.bonds):
            self.blocks[part] = rings
        else:
            self.blocks[part] = split(Molecule(None, list(visited_bonds), part))
