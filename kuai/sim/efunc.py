'''
Created on 2010-11-12

@author: Madlee
'''

from kuai.ff import FunctionBond, FunctionAngle
from kuai.ff import FunctionTorsion, FunctionFork2
from kuai.ff import FunctionPair13, FunctionPair14
from kuai.sim import ForceFieldEnergyProxy, NonbondEnergyProxy
from numpy import zeros, byte
from kuai.kuaisuan import screenPairList

def valecneFunctionCreator(name, model, key2index, filter):
    result = ForceFieldEnergyProxy()
    result.name = name

    type = name.split('/')[0]

    rank = 0
    if type == 'bond':
        rank = 2
        fffunc = FunctionBond()
    elif type == 'angle':
        rank = 3
        fffunc = FunctionAngle()
    elif type == 'torsion':
        rank = 4
        fffunc = FunctionTorsion()
    elif type == 'improper_torsion':
        fffunc = FunctionFork2()
    else:
        assert False
    
    lists = []
    start = 0
    for mol, n in model.mols:
        data = fffunc.setup(mol, key2index)
        for atom, _, _ in data:
            for i in range(len(atom)):
                atom[i] = mol.index(atom[i])
        natoms = len(mol.atoms)
        for _ in range(n):
            for atom, _, index in data:
                for j in atom:
                    lists.append(j+start)
                lists.append(index)
            start += natoms
    if lists:
        result.lists = lists    
        if rank != 0:
            screenPairList(filter, result, rank)
        return result
    else:
        return None
    


def pairCreator(name, model, key2index, filter):
    result = NonbondEnergyProxy()
    result.name = name
    type = name.split('/')[0]
    rank = 0
    if type == 'pair-14':
        fffunc = FunctionPair14()
        rank = 2
    elif type == 'pair-13':
        fffunc = FunctionPair13()
        rank = 2
    else:
        assert False
    
    lists = []
    start = 0
    for mol, n in model.mols:
        data = fffunc.setup(mol, key2index)
        for atom, _, _ in data:
            for i in range(len(atom)):
                atom[i] = mol.index(atom[i])
        natoms = len(mol.atoms)
        for _ in range(n):
            for atom, _, index in data:
                for j in atom:
                    lists.append(j+start)
                lists.append(index)
            start += natoms
    if lists:
        result.lists = lists
        if ':columnb-constant:' in key2index and hasattr(model, 'charges'):
            result.columnb_constant = key2index[':columnb-constant:']
    
        if rank != 0:
            screenPairList(filter, result, rank)
        return result
    else:
        return None

def pairNonbondCreator(name, model, key2index, filter):
    result = NonbondEnergyProxy()
    result.name = name
    result.filter = filter
    if ':columnb-constant:' in key2index and hasattr(model, 'charges'):
        result.columnb_constant = key2index[':columnb-constant:']
    result.par0 = key2index[':start:']
    return result
    
    
FUNCTION_CREATOR = {'bond/harmonic':valecneFunctionCreator,
                    'angle/harmonic':valecneFunctionCreator,
                    'torsion/cos_poly5':valecneFunctionCreator,
                    'improper_torsion/cos_poly5':valecneFunctionCreator,

                    'pair-14/lj12_6':pairCreator,
                    'pair-13/lj12_6':pairCreator,
                    
                    'pair-nonbond/lj12_6':pairNonbondCreator,
                   }

def make_efunc(model, ffindex):
    efunc = []
    natoms = len(model.coords) / 3
    nfilter = (natoms * (natoms - 1) / 2 + 7) / 8
    filter = zeros((nfilter), byte) 
    filter |= -1
    for k, v in ffindex.iteritems():
        if k in FUNCTION_CREATOR:
            funcI = FUNCTION_CREATOR[k](k, model, v, filter)
            if funcI:
                efunc.append(funcI)
    return efunc
