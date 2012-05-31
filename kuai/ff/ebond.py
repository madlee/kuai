'''
Created on 2010-10-15

@author: madlee
'''

from kuai.ff import FunctionBond
from kuai.ff.efunc import EnergyFunction
from sympy import * 


class FunctionBHARM(FunctionBond, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')
        e0, r0 = symbols('e0 r0')
        self.x = x1, y1, z1, x2, y2, z2 
        self.k = r0, e0

        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        r = sqrt(dx*dx+dy*dy+dz*dz)
        e = e0*(r-r0)*(r-r0)
        self.e = e

        EnergyFunction.__init__(self)

class FunctionBQUAR(FunctionBond, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')
        r0, k2, k3, k4 = symbols('r0, k2, k3, k4')
        self.x = x1, y1, z1, x2, y2, z2 
        self.k = r0, k2, k3, k4

        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        r = sqrt(dx*dx+dy*dy+dz*dz)
        dr = r - r0
        e = dr * dr *(k2 + dr*(k3 + k4*dr))
        self.e = e

        EnergyFunction.__init__(self)

