'''
Created on 2010-4-9

@author: Madlee
'''

from math import sqrt, acos, pi

class XYZ:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __abs__(self):
        return sqrt(self.x*self.x + self.y*self.y+self.z*self.z)
    
    def __str__(self):
        return "(%12.4f, %12.4f, %12.4f)" % (self.x, self.y, self.z)
    
    def __add__(self, v2):
        return XYZ(self.x+v2.x, self.y+v2.y, self.z+v2.z)
    
    def __sub__(self, v2):
        return XYZ(self.x-v2.x, self.y-v2.y, self.z-v2.z)
    
    def __mul__(self, v2):
        try:
            x = self.y*v2.z-self.z*v2.y
            y = self.z*v2.x-self.x*v2.z
            z = self.x*v2.y-self.y*v2.x
            return XYZ(x, y, z)
        except AttributeError:
            return XYZ(self.x*v2, self.y*v2, self.z*v2)
    
    def __div__(self, v2):
        return XYZ(self.x/v2, self.y/v2, self.z/v2)
    
    def __neg__(self):
        return XYZ(-self.x, -self.y, -self.z)
    
    def __eq__(self, v2):
        return self.x == v2.x and self.y == v2.y and self.z == v2.z
    
    def __ne__(self, v2):
        return not (self == v2)
    
    def __repr__(self):
        return str(self)

def dot(v1, v2):
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z

def angle(v1, v2, v3=None):
    if v3 is not None:
        v1 = v1-v2
        v2 = v3-v2
    r = dot(v1, v2) / sqrt(dot(v1, v1) * dot(v2, v2))
    if r <= -1.0:
        return pi
    elif r >= 1:
        return 0
    else:
        return acos(r)
        
def torsion(v1, v2, v3, v4):
    axis = v3-v2
    v1 = (v1 - v2) * axis
    v4 = (v3 - v4) * axis
    return angle(v1, v4)

