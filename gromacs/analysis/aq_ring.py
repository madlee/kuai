"""
Usage:
    aq-ring.py input.pdb
"""

from kuai.types.mytypes import Atom, Bond
from kuai.basicTools.xyz import XYZ
from kuai.basicTools.pbc import PBC, radians
from kuai.basicTools.top import angle

THRESHOLD = 3.5
MAX_ANGLE = radians(30)

MAX_RING_SIZE = 20

SINGLE_BOND = 10
HYDROGEN_BOND = 0


class Model:
    def __init__(self, pbc, atoms, bonds):
        self.pbc = pbc
        self.atoms = atoms[:]
        self.bonds = bonds[:]
        self.__adj = {}
        for i in self.atoms:
            self.__adj[i] = []
            
        for i in self.bonds:
            assert i.atom1 in self.__adj
            self.__adj[i.atom1].append((i.atom2, i))
            assert i.atom2 in self.__adj
            self.__adj[i.atom2].append((i.atom1, i))
            
    def degree(self, atom):
        return len(self.__adj[atom])
    
    def index(self, var):
        try:
            return self.atoms.index(var)
        except ValueError:
            return self.bonds.index(var)
        
    def neighbor_atom(self, atom, i):
        return self.__adj[atom][i][0]
    
    def neighbor_bond(self, atom, i):
        return self.__adj[atom][i][1]
    
    def get_bond(self, atom1, atom2):
        for j in range(self.degree(atom1)):
            if self.neighbor_atom(atom1, j) is atom2:
                return self.neighbor_bond(atom1, j)
        return None
    
    
def read_pdb(filename):
    """Read atoms from pdb file. Modify this to coop with your input file."""
    pbc = None
    atoms = []
    bonds = []
    with open(filename) as file:
        last_o = None
        for line in file:
            if line.startswith('ATOM'):
                element = line[13:17].strip()
                if element == 'OW':
                    last_o = Atom.parse(line)
                    atoms.append(last_o)
                elif element == 'HW1' or element == 'HW2':
                    h = Atom.parse(line)
                    atoms.append(h)
                    bonds.append(Bond(last_o, h, SINGLE_BOND))
                    
            elif line.startswith('CRYST1'):
                line = ' '.join(line.split()[1:7])
                pbc = PBC.parse(line)
                
    return Model(pbc, atoms, bonds)

def guess_hydrogen_bond(model):
    result = []
    
    o = [i for i in model.atoms if i.symbol == 'OW']  # All oxygen
    n = len(o)
    
    for i in xrange(n):
        oi = o[i]
        hi1 = model.neighbor_atom(oi, 0)
        delta_hi1 = model.pbc.norm(oi.coords - hi1.coords)
        hi2 = model.neighbor_atom(oi, 1)
        delta_hi2 = model.pbc.norm(oi.coords - hi2.coords)
        for j in xrange(i):
            oj = o[j]
            delta_oo = model.pbc.norm(oi.coords-oj.coords)
            if abs(delta_oo) <= THRESHOLD:
                hj1 = model.neighbor_atom(oj, 0)
                delta_hj1 = model.pbc.norm(hj1.coords - oj.coords)
                hj2 = model.neighbor_atom(oj, 1)
                delta_hj2 = model.pbc.norm(hj2.coords - oj.coords)
                
                if angle(delta_hi1, delta_oo) < MAX_ANGLE:
                    result.append(Bond(hi1, oj, HYDROGEN_BOND))
                    
                if angle(delta_hi2, delta_oo) < MAX_ANGLE:
                    result.append(Bond(hi2, oj, HYDROGEN_BOND))
                    
                if angle(delta_hj1, delta_oo) < MAX_ANGLE:
                    result.append(Bond(hj1, oi, HYDROGEN_BOND))
                    
                if angle(delta_hj2, delta_oo) < MAX_ANGLE:
                    result.append(Bond(hj2, oi, HYDROGEN_BOND))
                
    return result            


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
    def __init__(self, bond):
        self.bond = bond
        self.trace = {bond.atom2:(None, 2)}
    
    def start(self, mol, atom):
        if self.trace[atom][1] > MAX_RING_SIZE:
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


if __name__ == '__main__':
    import os
    os.chdir("/home/sungroup/Documents/chengtao/sim/all-atom-md/uni_2/spc")
    
    for i in range(40):
        model = read_pdb("all%d.pdb"%i)
        output = 'ring%02d.txt'%i
        o = open(output, 'w')
        hydrogen_bonds = guess_hydrogen_bond(model)
        
        model = Model(model.pbc, model.atoms, model.bonds+hydrogen_bonds)
        
        for bond in hydrogen_bonds:
            ring = RingFinder.find_ring(model, bond)
            if ring != None:
                o.write(' '.join([str(i) for i in ring])+'\n') 
        o.close()
            

        