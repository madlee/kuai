'''
Created on 2010-10-15

@author: Madlee
'''
from kuai.ff import FunctionAngle
from kuai.ff.efunc import EnergyFunction
from sympy import * 

class FunctionAHARM(FunctionAngle, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2, x3, y3, z3 = symbols('x1 y1 z1 x2 y2 z2 x3 y3 z3')
        e0, theta0 = symbols('e0 theta0')
        self.x = x1, y1, z1, x2, y2, z2, x3, y3, z3 
        self.k = theta0, e0

        dx1, dy1, dz1 = x1-x2, y1-y2, z1-z2
        dx2, dy2, dz2 = x3-x2, y3-y2, z3-z2
        
        r12 = sqrt((dx1*dx1+dy1*dy1+dz1*dz1)*(dx2*dx2+dy2*dy2+dz2*dz2))
        theta = acos((dx1*dx2+dy1*dy2+dz1*dz2)/r12)

        e = e0*(theta-theta0)*(theta-theta0)
        self.e = e

        EnergyFunction.__init__(self)

class FunctionAQUAR(FunctionAngle, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2, x3, y3, z3 = symbols('x1 y1 z1 x2 y2 z2 x3 y3 z3')
        k2, k3, k4, theta0 = symbols('k2 k3 k4 theta0')
        self.x = x1, y1, z1, x2, y2, z2, x3, y3, z3 
        self.k = theta0, k2, k3, k4

        dx1, dy1, dz1 = x1-x2, y1-y2, z1-z2
        dx2, dy2, dz2 = x3-x2, y3-y2, z3-z2
        
        r12 = sqrt((dx1*dx1+dy1*dy1+dz1*dz1)*(dx2*dx2+dy2*dy2+dz2*dz2))
        theta = acos((dx1*dx2+dy1*dy2+dz1*dz2)/r12)
        dt = theta-theta0
        e = dt * dt *(k2 + dt*(k3 + k4*dt))
        self.e = e

        EnergyFunction.__init__(self)

class FunctionACOSH(FunctionAngle, EnergyFunction):
    def __init__(self):
        x1, y1, z1, x2, y2, z2, x3, y3, z3 = symbols('x1 y1 z1 x2 y2 z2 x3 y3 z3')
        e0, cos0 = symbols('e0 cos0')
        self.x = x1, y1, z1, x2, y2, z2, x3, y3, z3 
        self.k = cos0, e0

        dx1, dy1, dz1 = x1-x2, y1-y2, z1-z2
        dx2, dy2, dz2 = x3-x2, y3-y2, z3-z2
        
        r12 = sqrt((dx1*dx1+dy1*dy1+dz1*dz1)*(dx2*dx2+dy2*dy2+dz2*dz2))
        cost = (dx1*dx2+dy1*dy2+dz1*dz2)/r12
        e = e0 * (cost - cos0) * (cost - cos0)
        self.e = e

        EnergyFunction.__init__(self)

