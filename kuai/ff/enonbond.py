'''
Created on 2010-10-16

@author: Madlee
'''

from kuai.ff import FunctionNonbond
from kuai.ff.efunc import EnergyFunction
from sympy import * 


class FunctionN12_6(FunctionNonbond, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')
        r0, e0 = symbols('r0 e0')
        self.x = x1, y1, z1, x2, y2, z2 
        self.k = r0, e0

        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        r2 = (r0*r0)/(dx*dx+dy*dy+dz*dz)
        r6 = r2 * r2 * r2
        e = e0 * r6 * (r6 - 2)
        self.e = e

        EnergyFunction.__init__(self)

class FunctionN_9_6(FunctionNonbond, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')
        r0, e0 = symbols('r0 e0')
        self.x = x1, y1, z1, x2, y2, z2 
        self.k = r0, e0

        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        r = sqrt(dx*dx+dy*dy+dz*dz)
        r = r0/r
        r3 = r * r * r
        e = e0 * r3 * r3 * (2 * r3 - 3)
        self.e = e

        EnergyFunction.__init__(self)

class FunctionCoulomb(EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')
        r0, e0 = symbols('r0 e0')
        self.x = x1, y1, z1, x2, y2, z2 
        self.k = r0, e0

        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        r = sqrt(dx*dx+dy*dy+dz*dz)
        r = r0/r
        r3 = r * r * r
        e = e0 * r3 * r3 * (2 * r3 - 3)
        self.e = e

        EnergyFunction.__init__(self)
