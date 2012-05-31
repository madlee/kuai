'''
Created on 2010-10-15

@author: Madlee
'''

from sympy import diff, N
from numpy import array, matrix

class EnergyFunction:
    def __init__(self):
        e = self.e
        dedx = []
        n = len(self.x)
        for i in self.x:
            dedx.append(diff(e, i))
        self.dedx = dedx
        
        d2edx2 = []
        for i in range(n):
            func = dedx[i]
            dfunc = []
            for j in range(n):
                dfunc.append(diff(func, self.x[j]))
            d2edx2.append(dfunc)
        self.d2edx2 = d2edx2
        
    def setup_data(self, k, x, list):
        natom = len(list)-1
        data = {}
        for i in range(natom):
            atom = list[i]
            data[self.x[i*3+0]] = x[atom*3+0]
            data[self.x[i*3+1]] = x[atom*3+1]
            data[self.x[i*3+2]] = x[atom*3+2]

        ik = list[natom]
        npar = len(self.k)
        for i in range(npar):
            data[self.k[i]] = k[ik+i]
        return data

    def gete(self, k, x, list):
        data = self.setup_data(k, x, list)
        result = N(self.e(data))
        return float(result)
    
    def getdedx(self, k, x, list):
        data = self.setup_data(k, x, list)
        result = []
        for i in self.dedx:
            vi = N(i(data))
            result.append(float(vi))
        return result
    
    def getd2edx2(self, k, x, list):
        data = self.setup_data(k, x, list)
        result = []
        for i in self.d2edx2:
            vi = []
            for j in i:
                vj = float(N(j(data)))
                vi.append(vj)
            result.append(vi)
        return result

    def energy(self, k, x):
        result = 0
        rank = self.rank()
        n = len(self.list)
        for i in range(0,n,rank+1):
            result += self.gete(k, x, self.list[i:i+rank+1])
        return result
             
    def gradian(self, par, x):
        result = [0.0] * 3 * self.natoms
        rank = self.rank()
        n = len(self.list)
        assert n % (rank+1) == 0
        for i in range(0,n,rank+1):
            templist = self.list[i:i+rank+1]
            v = self.getdedx(par, x, templist)
            for j in range(rank):
                atom = templist[j]
                result[atom*3+0] += v[j*3+0]
                result[atom*3+1] += v[j*3+1]
                result[atom*3+2] += v[j*3+2]
        return array(result)
             

    def hessian(self, par, x):
        result = [[0.0] * (3 * self.natoms) for i in range(3 * self.natoms)]
        rank = self.rank()
        n = len(self.list)
        assert n % (rank+1) == 0
        for i in range(0, n, rank+1):
            templist = self.list[i:i+rank+1]
            v = self.getd2edx2(par, x, templist)
            for j in range(rank):
                atomJ = templist[j]
                for k in range(rank):
                    atomK = templist[k]
                    result[atomJ*3+0][atomK*3+0] += v[j*3+0][k*3+0]
                    result[atomJ*3+0][atomK*3+1] += v[j*3+0][k*3+1]
                    result[atomJ*3+0][atomK*3+2] += v[j*3+0][k*3+2]
                    result[atomJ*3+1][atomK*3+0] += v[j*3+1][k*3+0]
                    result[atomJ*3+1][atomK*3+1] += v[j*3+1][k*3+1]
                    result[atomJ*3+1][atomK*3+2] += v[j*3+1][k*3+2]
                    result[atomJ*3+2][atomK*3+0] += v[j*3+2][k*3+0]
                    result[atomJ*3+2][atomK*3+1] += v[j*3+2][k*3+1]
                    result[atomJ*3+2][atomK*3+2] += v[j*3+2][k*3+2]
        return matrix(result)

class EnergySet:
    def __init__(self):
        self.funcs = []     # all (function, list) pair
        
    def energy(self, k, x):
        result = 0
        for foo in self.funcs:
            result += foo.energy(k, x) 
        return result
             
    def gradian(self, k, x):
        result = 0
        for foo in self.funcs:
            result += foo.gradian(k, x) 
        return result
             

    def hessian(self, k, x):
        result = 0
        for foo in self.funcs:
            result += foo.hessian(k, x)
        return result
