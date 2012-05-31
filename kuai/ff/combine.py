'''
Created on 2010-11-19

@author: Madlee
'''
from math import sqrt

def list_atomtypes(iter):
    types = []
    for i in iter:
        types += i.split()
    types = set(types)
    types =[i for i in types]
    types.sort()
    type2id = {}
    for i in range(len(types)):
        type2id[types[i]] = i
    return types, type2id

def combine_arithmatic(index, par, typeI, typeJ):
    assert typeJ <= typeI
    ii = index[typeI]
    r0i = par[ii]
    e0i = par[ii+1]

    jj = index[typeJ]
    r0j = par[jj]
    e0j = par[jj+1]
    return (r0i+r0j)/2, sqrt(e0i*e0j), 
    
def combine_geometry(index, par, typeI, typeJ):
    assert typeJ <= typeI
    ii = index[typeI]
    r0i = par[ii]
    e0i = par[ii+1]

    jj = index[typeJ]
    r0j = par[jj]
    e0j = par[jj+1]
    return sqrt(r0i*r0j), sqrt(e0i*e0j), 

def combine(index, par, function):
    types, type2id = list_atomtypes(index.iterkeys())
    newpar = []
    i0 = len(par)
    newindex = {}
    n = len(types)
    for i in range(n):
        typeI = types[i]
        for j in range(0, i+1):
            typeJ = types[j]
            key = ' '.join([typeJ, typeI])
            parIJ = function(index, par, typeI, typeJ)
            if key in index:
                newpar += par[index[key]:index[key]+len(parIJ)]
            else:
                newpar += function(index, par, typeI, typeJ)
            newindex[key] = i0
            i0 += len(parIJ)
    
    newindex[':start:'] = len(par)
    par += newpar
    newindex[':finish:'] = len(par)
    return type2id, newindex, 

def scale_vdw(index, par, scale_vdw, offsets):
    if scale_vdw == 1.0:
        return index
    else:
        i0 = index[':start:']
        i1 = index[':finish:']
        new_i0 = len(par) 
        new_par = par[i0:i1]
        newindex = {}
        for k, v in index.iteritems():
            if i0 <= v and v < i1:
                newindex[k] = v - i0
                try:
                    new_par[newindex[k] + offsets] *= scale_vdw
                except TypeError:
                    for j in offsets:
                        new_par[newindex[k] + j] *= scale_vdw
                newindex[k] += new_i0
        par += new_par
        newindex[':finish:'] = len(par)
        return newindex
        
    