'''
Created on 2010-4-9

@author: madlee
'''

'''
Created on 2010-3-25

@author: Madlee
'''

from kuai.xyz import XYZ
from kuai.fit import FittingFunction
from kuai.ff.charges import ChargeFunction

class ChargeFittingData(FittingFunction):
    def __init__(self):
        self.charges = []
        self.dipoles = []
        self.esps = []
        self.func = {}
        self.fixed_parameters = {}
        
    def add_function(self, func, par=None):
        if isinstance(func, ChargeFunction):
            if par is None:
                par = {}
            self.func[func] = par
        else:
            raise TypeError("A ChargeFunction is required")
        
    def add_charge(self, mol, charges):
        assert (len(mol.atoms) == len(charges))
        self.charges.append((mol, charges))
    
    def add_dipole(self, mol, dipole):
        self.dipoles.append((mol, dipole, ))
    
    def add_esp(self, mol, esp, sigma):
        self.esps.append((mol, esp, ))
        
    def get_y0(self):
        result = []
        for mol, q in self.charges:
            result += q
        for mol, dipole in self.dipoles:
            result.append(dipole.x)
            result.append(dipole.y)
            result.append(dipole.z)
        for mol, esp in self.esps:
            for x, e in esp:
                result.append(e)
        return result

    def count_y(self):
        result = 0
        for mol, q in self.charges:
            result += len(q)
        for mol, dipole in self.dipoles:
            result += 3
        for mol, esp in self.esps:
            result += len(esp)
        return result
        
    def makeup_k(self):
        for func, par in self.func.iteritems():
            for mol, q in self.charges:
                func.makeup(mol, par)
            for mol, dipole in self.dipoles:
                func.makeup(mol, par)
            for mol, esp in self.esps:
                func.makeup(mol, par)
        
    def get_par2index(self, offset=0):
        result = {}
        for func, par in self.func.iteritems():
            par2index = {}
            for k, v in par.iteritems():
                par2index[k] = offset
                offset += len(v)
            result[func] = par2index
        return result
    
    def get_k0(self):
        result = []
        for func, par in self.func.iteritems():
            for k, v in par.iteritems():
                result += v
        return result

    def count_k(self):
        result = 0
        for func, par in self.func.iteritems():
            for k, v in par.iteritems():
                result += len(v)
        return result

    def calc_y(self):
        result = []
        for mol, q in self.charges:
            q = self.__get_q(mol)
            result += q
        for mol, dipole in self.dipoles:
            q = self.__get_q(mol)
            dipole = XYZ(0.0, 0.0, 0.0)
            for i in len(q):
                dipole += mol.atoms[i].coord * q[i]
            result.append(dipole.x)
            result.append(dipole.y)
            result.append(dipole.z)
        for mol, esp in self.esps:
            q = self.__get_q(mol)
            for x, e in esp:
                e = 0.0
                for i in len(q):
                    r = abs(x - mol.atoms[i].coord)
                    e += q[i]/r
                result.append(e)
        return result
    
    def calc_dydk(self, par2index):
        result = []
        nk = self.count_k()
        for mol, q in self.charges:
            dqdk = self.__get_dqdk(mol, par2index)
            result += dqdk
        for mol, dipole in self.dipoles:
            dqdk = self.__get_dqdk(mol, par2index)
            natoms = len(mol.atoms)
            assert len(dqdk) == natoms
            x, y, z = [0.0] * nk, [0.0] * nk, [0.0] * nk
            for i in range(natoms):
                for j in range(nk):
                    x[j] += dqdk[i][j] * mol.atoms[i].x
                    y[j] += dqdk[i][j] * mol.atoms[i].y
                    z[j] += dqdk[i][j] * mol.atoms[i].z
            result += [x, y, z]
        for mol, esp in self.esps:
            dqdk = self.__get_dqdk(mol, par2index)
            natoms = len(mol.atoms)
            assert len(dqdk) == natoms
            for x, e in esp:
                dedk = [0.0] * nk
                for i in range(natoms):
                    r = abs(mol.atoms[i].coord - x)
                    for j in range(nk):
                        dedk[j] += dqdk[i][j] / r 
                result.append(dedk)
        return result
    
    def get_fixed_index(self, par2index):
        result = []
        for func, items in self.fixed_parameters.iteritems():
            for key, indexes in items.iteritems():
                offset = par2index[func][key]
                for i in indexes:
                    result.append(offset+i)
        return result
    
    def reset_k(self, newk):
        offset = 0
        for func, par in self.func.iteritems():
            for k, v in par.iteritems():
                for i in range(len(v)):
                    v[i] = newk[offset]
                    offset += 1
                    
    def fix_function(self, func):
        self.fixed_parameters[func] = {}
        for k, par in self.func[func].iteritems():
            self.fixed_parameters[func][k] = [i for i in range(len(par))]
    
    def __get_q(self, mol):
        natoms = len(mol.atoms)
        q = [0.0] * natoms
        for func, par in self.func.iteritems():
            charges = func.get_charges(mol, par)
            assert len(charges) == natoms
            for j in range(natoms):
                q[j] += charges[j]
        return q
    
    def __get_dqdk(self, mol, par2index):
        nk = self.count_k()
        natoms = len(mol.atoms)
        result = []
        for i in range(natoms):
            result.append([0.0] * nk)
        for func, par in self.func.iteritems():
            dqdk = func.get_dqdk(mol, par, par2index[func], nk)
            assert len(dqdk) == natoms
            for i in range(natoms):
                atomI = mol.atoms[i]
                for j in range(nk):
                    result[i][j] += dqdk[atomI][j]
        return result
        
