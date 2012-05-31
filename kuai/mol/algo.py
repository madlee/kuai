'''
Created on 2010-6-8

@author: Madlee
'''

from kuai import mol

def guess_bond_link(atoms, pbc=None, vdw_radius=None, tor = 0.2):
    if vdw_radius is None:
        vdw_radius = {}
    for i in atoms:
        if i.symbol not in vdw_radius:
            vdw_radius[i.symbol] = i.element.radius

    bonds = []
    natoms = len(atoms)
    for i, atomI in enumerate(atoms):
        for j in range(i+1, natoms):
            atomJ = atoms[j]
            r = atomI.coord - atomJ.coord
            if pbc:
                r = pbc.norm(r)
            if abs(r) < vdw_radius[atomI.symbol] + vdw_radius[atomJ.symbol] + tor:
                bonds.append(mol.Bond(atomI, atomJ))
    return bonds


class MoleculeVisitor(Exception):
    FLAG_WHITE, FLAG_GRAY, FLAG_BLACK = 0, 1, 2 
    
    def start(self, mol, atom):
        pass
    
    def finish(self, mol, atom):
        pass
    
    def find(self, mol, target, source):
        pass
    
    def back(self, mol, target, source):
        pass

def dft(mol, visitor, flags=None, stack=None):
    if flags is None:
        flags = dict([(i, MoleculeVisitor.FLAG_WHITE,) for i in mol.atoms])
    if stack is None:
        stack = [mol.atoms[0]]
    while len(stack) > 0:
        i = stack.pop()
        visitor.start(mol, i)
        flags[i] = MoleculeVisitor.FLAG_GRAY
        ndegree = mol.degree(i)
        for j in range(ndegree):
            atomJ = mol.neighbor_atom(i, j)
            if flags[atomJ] == MoleculeVisitor.FLAG_WHITE:
                visitor.find(mol, atomJ, i)
                stack.append(atomJ)
                dft(mol, visitor, flags, stack)
            else:
                visitor.back(mol, atomJ, i)
        flags[i] = MoleculeVisitor.FLAG_BLACK
        visitor.finish(mol, i)
        
def bft(mol, visitor, flags=None, queue=None):
    if flags is None:
        flags = dict([(i, MoleculeVisitor.FLAG_WHITE,) for i in mol.atoms])
    if queue is None:
        queue  = [mol.atoms[0]]
    for i in queue:
        visitor.start(mol, i)
        ndegree = mol.degree(i)
        for j in range(ndegree):
            atomJ = mol.neighbor_atom(i, j)
            if flags[atomJ] == MoleculeVisitor.FLAG_WHITE:
                visitor.find(mol, atomJ, i)
                flags[atomJ] = MoleculeVisitor.FLAG_GRAY
                queue.append(atomJ)
            else:
                visitor.back(mol, atomJ, i)
        flags[i] = MoleculeVisitor.FLAG_BLACK
        visitor.finish(mol, i)

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
    def find_ring(mol, bond):
        finder = RingFinder(bond)
        try:
            flags = dict([(i, MoleculeVisitor.FLAG_WHITE,) for i in mol.atoms])
            flags[bond.atom1] = MoleculeVisitor.FLAG_BLACK
            bft(mol, finder, flags, [bond.atom2])
            return None
        except RingFinder:
            return finder.get_result()

class ConnectedPartVisitor(MoleculeVisitor):
    def __init__(self, atom):
        self.atoms = []
        self.atoms.append(atom)
        
    def find(self, mol, target, source):
        self.atoms.append(target)


def split(v):
    flags = {}
    for i in v.atoms:
        flags[i] = MoleculeVisitor.FLAG_WHITE

    result = []
    for i in v.atoms:
        if flags[i] == MoleculeVisitor.FLAG_WHITE:
            visitor = ConnectedPartVisitor(i)
            dft(v, visitor, flags, [i])
            result.append(mol.submol(v, visitor.atoms, None))

    return result
        
