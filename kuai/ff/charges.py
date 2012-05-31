'''
Created on 2010-4-9

@author: madlee
'''
'''
Created on 2010-3-25

@author: Madlee
'''

from kuai.xyz import XYZ, torsion
from math import cos
from kuai.ff import FunctionAtom, FunctionBond, FunctionAngle, FunctionTorsion

class ChargeFunction:
    def calc_charge(self, mol, atoms, par):
        raise NotImplementedError() 
    
    def calc_dqdk(self, mol, atoms, par):
        raise NotImplementedError()
    
    def get_charges(self, mol, ff, charge0 = None):
        result = {}
        for i in mol.atoms:
            result[i] = 0.0
        
        lists = self.setup(mol, ff)
        for atoms, key, par in lists:
            q = self.calc_charge(mol, atoms, par)
            assert len(q) == len(atoms)
            for j in range(len(atoms)):
                result[atoms[j]] += q[j]
        if charge0 is None:
            charge0 = [0.0] * len(mol.atoms)
        else:
            assert len(charge0) == len(mol.atoms)
        for i in range(len(charge0)):
            charge0[i] += result[mol.atoms[i]]
        return charge0
    
    def get_dipole(self, mol, ff, dipole0 = None):
        q = self.get_charges(mol, ff)
        assert len(q) == len(mol.atoms)
        result = XYZ(0.0, 0.0, 0.0)
        for i in range(len(q)):
            result += mol.atoms[i].coord * q[i]
        if dipole0 is not None:
            result += dipole0
        return result
    
    def get_esp(self, mol, ff, coord):
        q = self.get_charges(mol, ff)
        assert len(q) == len(mol.atoms)
        result = 0
        for i in range(len(q)):
            r = abs(mol.atoms[i].coord - coord)
            result += q[i] / r
        return i
    
    def get_dqdk(self, mol, ff, par2index, n):
        result = {}
        for i in mol.atoms:
            result[i] = [0.0] * n

        lists = self.setup(mol, ff)
        for atoms, key, par in lists:
            assert len(atoms) == self.rank()
            rank = len(atoms)
            dqdk = self.calc_dqdk(mol, atoms, par)
            offset = par2index[key]
            assert len(dqdk) == rank
            for j in range(rank):
                assert len(dqdk[j]) == len(par)
                for k in range(len(par)):
                    result[atoms[j]][k+offset] += dqdk[j][k]
        return result
    
class FunctionQATC(FunctionAtom, ChargeFunction):
    def is_zero(self, par):
        return par[0] == 0
                
    def calc_charge(self, mol, atoms, par):
        return par[0],
    
    def calc_dqdk(self, mol, atoms, par):
        return (1.0,),
    
    def estimate(self, key):
        return [0.0];
    
     
class FunctionQBINC(FunctionBond, ChargeFunction):
    def is_zero(self, par):
        return par[0] == 0
                
    def calc_charge(self, mol, atoms, par):
        return par[0], -par[0],
    
    def calc_dqdk(self, mol, atoms, par):
        return (1.0,), (-1.0,),
    
    def estimate(self, key):
        return [0.0];

class FunctionQAINC(FunctionAngle, ChargeFunction):
    def is_zero(self, par):
        return par[0] == 0
                
    def calc_charge(self, mol, atoms, par):
        return par[0], 0.0, -par[0],
    
    def calc_dqdk(self, mol, atoms, par):
        return (1.0,), (0.0,), (-1.0,),
    
    def estimate(self, key):
        return [0.0];

class FunctionQTINC(FunctionTorsion, ChargeFunction):
    def is_zero(self, par):
        return par[0] == 0
                
    def calc_charge(self, mol, atoms, par):
        return par[0], 0.0, 0.0, -par[0],
    
    def calc_dqdk(self, mol, atoms, par):
        return (1.0,), (0.0,), (0.0,), (-1.0,),
    
    def estimate(self, key):
        return [0.0];


class FunctionQTCOS(FunctionTorsion, ChargeFunction):
    def is_zero(self, par):
        return par[0] == 0 and par[1] == 0 and par[2] == 0 and par[3] == 0    
                
    def calc_charge(self, mol, atoms, par):
        phi = torsion(atoms[0].coord, atoms[1].coord, atoms[2].coord, atoms[3].coord)
        delta = par[0] * cos(phi) + par[1] * cos(phi*2) + par[2] * cos(phi*3) + par[3]*cos(phi*4)
        return delta, 0.0, 0.0, -delta,
    
    def calc_dqdk(self, mol, atoms, par):
        phi = torsion(atoms[0].coord, atoms[1].coord, atoms[2].coord, atoms[3].coord)
        d1, d2, d3, d4 = cos(phi), cos(2*phi), cos(3*phi), cos(4*phi)
        return (d1, d2, d3, d4,), (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0,), (-d1, -d2, -d3, -d4),
    
    def estimate(self, key):
        return [0.0, 0.0, 0.0, 0.0];

class FunctionQCBT(FunctionTorsion, ChargeFunction):
    def is_zero(self, par):
        for i in par:
            if i != 0:
                return False
        return True    
                
    def calc_charge(self, mol, atoms, par):
        phi = torsion(atoms[0].coord, atoms[1].coord, atoms[2].coord, atoms[3].coord)
        delta1 = par[0] * cos(phi) + par[1] * cos(phi*2) + par[2] * cos(phi*3) + par[3]*cos(phi*4)
        delta2 = par[4] * cos(phi) + par[5] * cos(phi*2) + par[6] * cos(phi*3) + par[7]*cos(phi*4)
        return delta1, -delta1, -delta2, delta2,
    
    def calc_dqdk(self, mol, atoms, par):
        phi = torsion(atoms[0].coord, atoms[1].coord, atoms[2].coord, atoms[3].coord)
        d1, d2, d3, d4 = cos(phi), cos(2*phi), cos(3*phi), cos(4*phi)
        return  (d1, d2, d3, d4, 0.0, 0.0, 0.0, 0.0), \
                (-d1, -d2, -d3, -d4, 0.0, 0.0, 0.0, 0.0), \
                (0.0, 0.0, 0.0, 0.0, -d1, -d2, -d3, -d4), \
                (0.0, 0.0, 0.0, 0.0, d1, d2, d3, d4),
    
    def estimate(self, key):
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

