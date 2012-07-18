"""This package include model of molecules"""
from kuai.mol.element import Element 
from kuai.xyz import XYZ
    
class Atom(object):
    def __init__(self, symbol='C', charge=0):
        try:
            self.__symbol = Element.get(symbol).symbol
        except KeyError:
            self.__symbol = symbol
            
        self.charge = charge

    @property
    def symbol(self):
        return self.__symbol
    
    @property
    def element(self):
        return Element.get(self.symbol)
    
    @property
    def number(self):
        return self.element.number
    
    @property
    def weight(self):
        return self.element.weight
        
class Bond(object):
    UNKNOWN_BOND    =  0
    SINGLE_BOND     = 10
    PARTIAL_BOND    = 15
    DOUBLE_BOND     = 20
    TRIPLE_BOND     = 30

    __symbol = {UNKNOWN_BOND:'?', SINGLE_BOND:'-', PARTIAL_BOND:'~', DOUBLE_BOND:'=', TRIPLE_BOND:'#'}

    def __init__(self, atom1, atom2, order=2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        
    @property
    def symbol(self):
        return Bond.__symbol[self.order]

class Molecule(object):
    def __init__(self, atoms=None, bonds=None, parent=None):
        if parent != None:
            if atoms == None and bonds == None:
                atoms = parent.atoms
                bonds = parent.bonds
            elif atoms == None:
                atomset = set()
                for i in bonds:
                    atomset.add(i.atom1)
                    atomset.add(i.atom2)
                atoms = list(atomset)
            elif bonds == None:
                atomset = set(atoms)
                bonds = []
                for i in parent.bonds:
                    if i.atom1 in atomset and i.atom2 in atomset:
                        bonds.append(i)

        self.atoms = atoms
        self.bonds = bonds
        
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
        
    def __contains__(self, v):
        return v in self.atoms or v in self.bonds
        
    def neighbor_atom(self, atom, i):
        return self.__adj[atom][i][0]
    
    def neighbor_bond(self, atom, i):
        return self.__adj[atom][i][1]
    
    def get_bond(self, atom1, atom2):
        for j in range(self.degree(atom1)):
            if self.neighbor_atom(atom1, j) is atom2:
                return self.neighbor_bond(atom1, j)
        return None
    
    def get_angle(self, atom1, atom3):
        for m in range(self.degree(atom1)):
            atomM = self.neighbor_atom(atom1, m)
            for n in range(self.degree(atom3)):
                atomN = self.neighbor_atom(atom3, n)
                if  atomM is atomN: 
                    return atom1, atomM, atom3, 
        return None

    def rank_of_atom(self, atom):
        rank = 0
        for j in range(self.degree(atom)):
            rank += self.neighbor_bond(atom, j).order
        return rank
    
    def implicit_H(self, atom):
        e = atom.element
        if e.no_H:
            return 0
        else:
            try:
                val = e.valence[atom.charge][0]
                rank = self.rank_of_atom(atom)
                nH = (val * Bond.SINGLE_BOND - rank) / Bond.SINGLE_BOND
                if nH < 0:
                    nH = 0
                return nH
            except IndexError:
                return 0
    
    @property
    def weight(self):
        result = 0.0
        nH = 0
        for i in self.atoms:
            result += i.weight
            nH += self.implicit_H(i)
        result += nH * Element.get(1).weight
        return result
    
    @property
    def formula(self):
        map = {}
        nH = 0
        charge = 0
        for i in self.atoms:
            if i.symbol in map:
                map[i.symbol] += 1
            else:
                map[i.symbol] = 1
            nH += self.implicit_H(i)
            charge += i.charge
        if nH > 0:
            if 'H' in map:
                map['H'] += nH
            else:
                map['H'] = nH
                
        result = ''
        if 'C' in map:
            result += 'C'
            if map['C'] > 1:
                result += str(map['C'])
            del map['C']
        if 'H' in map:
            result += 'H'
            if map['H'] > 1:
                result += str(map['H'])
            del map['H']
            
        keys = [i for i in map.iterkeys()]
        keys.sort()
        for i in keys:
            result += i
            if map[i] > 1:
                result += str(map[i])
        if charge > 0:
            result += '+'
            if charge > 1:
                result += str(charge)
        elif charge < 0:
            result += '-'
            if charge < -1:
                result += str(-charge)
        return result
