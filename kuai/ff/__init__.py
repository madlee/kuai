from kuai.ff.inner_coords import *


class ForceFieldFunction:
    def list(self, mol):
        raise NotImplementedError()
    
    def uniform(self, key, index):
        raise NotImplementedError()
    
    def rank(self):
        raise NotImplementedError()
    
    def get_key(self, mol, atoms):
        raise NotImplementedError()
    
    def estimate(self, key):
        raise NotImplementedError()
    
    def list_keys(self, mol, list=None):
        if list is None:
            list = self.list(mol)
        n = self.rank()
        assert len(list) % n == 0
        
        result = set()
        i = 0
        while i < len(list):
            atoms = list[i:i+n]
            key = self.get_key(mol, atoms)
            result.add(key)
            i += n
        return [key for key in result]
    
    def setup(self, mol, ff, list=None):
        if list is None:
            list = self.list(mol)
        n = self.rank()
        assert len(list) % n == 0
        result = []
        i = 0
        while i < len(list):
            atoms = list[i:i+n]
            key = self.get_key(mol, atoms)
            if key in ff:
                result.append((atoms, key, ff[key], ))
            else:
                print "Warning: Missing parameters of " + key
            i += n
        return result
    
    def makeup(self, mol, par):
        keyset = self.list_keys(mol)
        for i in keyset:
            if i not in par:
                par[i] = self.estimate(i)
        return par

    def is_zero(self, par):
        return False
    
    @property
    def natoms(self):
        return self.rank()
    
class FunctionAtom(ForceFieldFunction):
    def list(self, mol):
        return [i for i in mol.atoms]
            
    def uniform(self, types, index):
        pass

    def rank(self):
        return 1
    
    def get_key(self, mol, atoms):
        assert len(atoms) == 1
        return atoms[0].type

class FunctionPair(ForceFieldFunction):
    def uniform(self, types, index):
        if types[0] > types[1]:
            types[0], types[1] = types[1], types[0]
            if index is not None:
                index[0], index[1] = index[1], index[0]
        return (types, index,)

    def rank(self):
        return 2
    
    def get_key(self, mol, atoms):
        assert len(atoms) == 2
        types = [atoms[0].type, atoms[1].type]
        types, atoms = self.uniform(types, atoms)
        return ' '.join(types)

class FunctionBond(FunctionPair):
    def list(self, mol):
        return list_bond_atoms(mol)
    
    
class FunctionAngle(ForceFieldFunction):
    def list(self, mol):
        return list_angle_atoms(mol)
            
    def uniform(self, types, index):
        if types[0] > types[2]:
            types[0], types[2] = types[2], types[0]
            if index is not None:
                index[0], index[2] = index[2], index[0]
        return (types, index)

    def rank(self):
        return 3
    
    def get_key(self, mol, atoms):
        assert len(atoms) == 3
        types = [atoms[0].type, atoms[1].type, atoms[2].type]
        types, atoms = self.uniform(types, atoms)
        return ' '.join(types)

class FunctionPair13(FunctionPair):
    def list(self, mol):
        return list_pair_13_atoms(mol)
    
class FunctionTorsion(ForceFieldFunction):
    def list(self, mol):
        return list_torsion_atoms(mol)
            
    def uniform(self, types, index):
        if types[1] > types[2] or types[1] == types[2] and types[0] > types[3]:
            types[0], types[1], types[2], types[3] = types[3], types[2], types[1], types[0]
            if index is not None:
                index[0], index[1], index[2], index[3] = index[3], index[2], index[1], index[0] 
        return (types, index)

    def rank(self):
        return 4
    
    def get_key(self, mol, atoms):
        assert len(atoms) == 4
        types = [atoms[0].type, atoms[1].type, atoms[2].type, atoms[3].type]
        types, atoms = self.uniform(types, atoms)
        return ' '.join(types)


class FunctionPair14(FunctionPair):
    def list(self, mol):
        return list_pair_14_atoms(mol)

class FunctionFork0(ForceFieldFunction):
    def list(self, mol):
        return list_fork_0_atoms(mol)
            
    def uniform(self, types, index):
        if types[1] > types[2]:
            types[1], types[2] = types[2], types[1]
            index[1], index[2] = index[2], index[1]
        if types[2] > types[3]:
            types[2], types[3] = types[3], types[2]
            index[2], index[3] = index[3], index[2]
        if types[1] > types[2]:
            types[1], types[2] = types[2], types[1]
            index[1], index[2] = index[2], index[1]
        return types, index
        

    def rank(self):
        return 4
    
    def get_key(self, mol, atoms):
        assert len(atoms) == 4
        types = [atoms[0].type, atoms[1].type, atoms[2].type, atoms[3].type]
        types, atoms = self.uniform(types, atoms)
        return ' '.join(types)


class FunctionFork2(FunctionFork0):
    def list(self, mol):
        return list_fork_2_atoms(mol)
            
    def uniform(self, types, index):
        types[0], types[2] = types[2], types[0]  
        index[0], index[2] = index[2], index[0]
        FunctionFork0.uniform(self, types, index)
        types[0], types[2] = types[2], types[0]  
        index[0], index[2] = index[2], index[0]
        return (types, index)

class FunctionNonbond(ForceFieldFunction):
    def list(self, mol):
        result = []
        n = len(mol.atoms)
        for i in range(n):
            for j in range(i+1, n):
                result.append(i)
                result.append(j)

    def uniform(self, types, index):
        if types[0] > types[1]:
            types[0], types[1] = types[1], types[0]
            if index is not None:
                index[0], index[1] = index[1], index[0]
        return (types, index,)

    def rank(self):
        return 2
    
    def get_key(self, mol, atoms):
        assert len(atoms) == 2
        types = [atoms[0].type, atoms[1].type]
        types, atoms = self.uniform(types, atoms)
        return ' '.join(types)

def pack_parameters(ff, result=None):
    if result is None:
        result = []
    par2index = {}
    for k, v in ff.iteritems():
        par2index[k] = len(result)
        result += v
    return par2index, result

def max_index(par2index):
    result = 0
    for v in par2index.itervalues():
        if v > result:
            result = v
    return result
