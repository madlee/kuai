'''
Created on 2010-10-17

@author: madlee
'''

from math import sqrt
from scipy.linalg import eig

kcal_2_hz = 11791.0

def frequence(hess, mass):
    assert hess.shape == (len(mass)*3, len(mass)*3)
    natoms = len(mass)
    
    data = []
    for i in range(natoms):
        v1 = []
        v2 = []
        v3 = []
        for j in range(natoms):
            mw = kcal_2_hz/sqrt(mass[i]*mass[j])
            v1.append(hess[i*3+0, j*3+0]*mw)
            v1.append(hess[i*3+0, j*3+1]*mw)
            v1.append(hess[i*3+0, j*3+2]*mw)
            v2.append(hess[i*3+1, j*3+0]*mw)
            v2.append(hess[i*3+1, j*3+1]*mw)
            v2.append(hess[i*3+1, j*3+2]*mw)
            v3.append(hess[i*3+2, j*3+0]*mw)
            v3.append(hess[i*3+2, j*3+1]*mw)
            v3.append(hess[i*3+2, j*3+2]*mw)
        data.append(v1)
        data.append(v2)
        data.append(v3)
    la, _ = eig(data)
    result = []
    for i in la:
        i = i.real
        if i >= 0:
            result.append(sqrt(i))
        else:
            result.append(-sqrt(-i))
    result.sort(reverse=True)
    return result
