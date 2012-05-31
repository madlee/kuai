'''
Created on 2010-10-16

@author: Madlee
'''

from kuai.ff import FunctionTorsion
from kuai.ff.efunc import EnergyFunction
from sympy import * 

class FunctionTCOS5(FunctionTorsion, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2 = symbols('x1 y1 z1 x2 y2 z2')
        x3, y3, z3, x4, y4, z4 = symbols('x3 y3 z3 x4 y4 z4')
        e0, e1, e2, e3, e4, e5 = symbols('e0 e1, e2, e3, e4, e5')
        self.x = x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
        self.k = e0, e1, e2, e3, e4, e5

        vx1, vy1, vz1 = x1-x2, y1-y2, z1-z2
        vx2, vy2, vz2 = x4-x3, y4-y3, z4-z3
        
        ax, ay, az = x2-x3, y2-y3, z2-z3
        
        nx1, ny1, nz1 = vy1*az-vz1*ay, vz1*ax - vx1*az, vx1*ay - vy1*ax 
        nx2, ny2, nz2 = vy2*az-vz2*ay, vz2*ax - vx2*az, vx2*ay - vy2*ax

        r12 = sqrt((nx1*nx1+ny1*ny1+nz1*nz1)*(nx2*nx2+ny2*ny2+nz2*nz2))
        p = (nx1*nx2+ny1*ny2+nz1*nz2) / r12
        
        e = e0 + p*(e1 + p*(e2 + p*(e3 + p*(e4 +p*e5))))
        self.e = e

        EnergyFunction.__init__(self)    
