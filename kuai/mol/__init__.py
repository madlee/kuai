"""This package include model of molecules"""
from kuai.xyz import XYZ, angle
from kuai.mol.element import Element 
from math import cos, sin, acos, atan2, sqrt, pi, degrees, radians

class PBC:
    def __init__(self, a, b, c, alpha=pi/2, beta=pi/2, gamma=pi/2):
        cosA, cosB, cosG = cos(alpha), cos(beta), cos(gamma)
        sinG = sin(gamma)
        self.ax = a
        self.bx,  self.by = b*cosG, b*sinG
        self.cx, self.cy, = c*cosB, c*(cosA-cosB*cosG)/sinG 
        self.cz = sqrt(c*c-self.cx*self.cx - self.cy*self.cy)
    
    @property
    def a(self):
        return self.ax
    
    @property
    def b(self):
        return sqrt(self.bx*self.bx + self.by*self.by)
    
    @property
    def c(self):
        return sqrt(self.cx*self.cx + self.cy*self.cy + self.cz*self.cz)
    
    @property
    def va(self):
        return XYZ(self.ax, 0.0, 0.0)
    
    @property
    def vb(self):
        return XYZ(self.bx, self.by, 0.0)
    
    @property
    def vc(self):
        return XYZ(self.cx, self.cy, self.cz)
    
    @property
    def alpha(self):
        return angle(self.vb, self.vc)
        
    @property
    def beta(self):
        return acos(self.cx/self.c)
    
    @property
    def gamma(self):
        return atan2(self.by, self.bx)
    
    @property
    def volumn(self):
        return self.ax*self.by*self.cz
    
    def scale(self, v):
        try:
            self.x *= v.x
            self.y *= v.y
            self.z *= v.z
        except AttributeError:
            self.x *= v
            self.y *= v
            self.z *= v
            
    def norm(self, xyz):
        x = xyz.x
        y = xyz.y
        z = xyz.z
        while z < -self.cz/2:
            z += self.cz
            y += self.cy
            x += self.cx
        while z > self.cz/2:
            z -= self.cz
            y -= self.cy
            x -= self.cx
        while y < -self.by/2:
            y += self.by
            x += self.bx
        while y > self.by/2:
            y -= self.by
            x -= self.bx
        while x < -self.ax/2:
            x += self.ax
        while x > self.ax/2:
            x -= self.ax

        return XYZ(x, y, z)

    @staticmethod
    def parse(line):
        tokens = line.split()
        if len(tokens) >= 6:
            try:
                tokens = [float(i) for i in tokens[-6:]]
                tokens[3] = radians(tokens[3])
                tokens[4] = radians(tokens[4])
                tokens[5] = radians(tokens[5])
                return PBC(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5])
            except:
                pass
        return None
            
        
    def __str__(self):
        return "PBC: %12.4f %12.4f %12.4f %5.2f %5.2f %5.2f" % \
                (self.a, self.b, self.c, 
                 degrees(self.alpha), degrees(self.beta), degrees(self.gamma))
    
class Atom:
    def __init__(self, symbol='C'):
        try:
            self.__symbol = Element.get(symbol).symbol
        except KeyError:
            self.__symbol = symbol
        self.coord = XYZ(0.0, 0.0, 0.0)
        self.charge = 0
        self.type = '?'

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
        
class Bond:
    UNKNOWN_BOND = 0
    SINGLE_BOND = 2
    PARTIAL_BOND = 3
    DOUBLE_BOND = 4
    TRIPLE_BOND = 6
    
    __symbol = ['?', '?', '-', '~', '=', '?', '#']

    def __init__(self, atom1, atom2, order=2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        
    @property
    def symbol(self):
        return Bond.__symbol[self.order]

class Molecule:
    def __init__(self, atoms, bonds):
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
        rank = self.rank_of_atom(atom)
        nH = (atom.element.h_state * Bond.SINGLE_BOND - rank) / Bond.SINGLE_BOND
        if nH < 0:
            nH = 0
        return nH
    
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


def submol(mol, atoms, bonds):
    if atoms is None and bonds is None:
        atoms, bonds = mol.atoms, mol.bonds
    elif atoms is None:
        atomset = set()
        for i in bonds:
            atomset.add(i.atom1)
            atomset.add(i.atom2)
        atoms = [i for i in atomset]
    elif bonds is None:
        bonds = []
        for i in mol.bonds:
            if i.atom1 in atoms and i.atom2 in atoms:
                bonds.append(i)
    else:
        bonds2 = []
        for i in bonds:
            if i.atom1 in atoms and i.atom2 in atoms:
                bonds2.append(i)
        bonds = bonds2
    result = Molecule(atoms, bonds)
    return result
            