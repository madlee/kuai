from xyz import *
from top import angle
from math import cos, acos, sin, degrees, atan2, radians, pi


class PBC:
    """ class for the  Protein Data Bank (PDB) format
    """
    def __init__(self, a, b, c, alpha=pi/2, beta=pi/2, gamma=pi/2):
        """@param a: box a
        @param b: box b
        @param c: box c
        @param alpha: alpha angle default is pi/2
        @param beta: alpha angle default is pi/2
        @param gamma: alpha angle default is pi/2
        """
        cosA, cosB, cosG = cos(alpha), cos(beta), cos(gamma)
        """@var cosA: value of cos(alpha)
        @var cosB: value of cos(beta)
        @var cosG: value of cos(gamma)
        """
        sinG = sin(gamma)
        """@var sinG: value of sin(gamma)"""
        self.ax = a
        """@var ax: xx tensor"""
        self.bx,  self.by = b*cosG, b*sinG
        """@var bx: yx tensor"""
        """@var bx: yy tensor"""
        self.cx, self.cy, = c*cosB, c*(cosA-cosB*cosG)/sinG
        """@var cx: zx tensor"""
        """@var cy: zy tensor"""
        self.cz = sqrt(c*c-self.cx*self.cx - self.cy*self.cy)
    
    @property
    def a(self):
        return self.ax
    
    @property
    def b(self):
        return sqrt(self.bx*self.bx + self.by*self.by)
    
    @property
    def c(self):
        return sqrt(self.cx*self.cx + self.cy*self.cy + self.cz*self.cz)
    
    @property
    def va(self):
        return XYZ(self.ax, 0.0, 0.0)
    
    @property
    def vb(self):
        return XYZ(self.bx, self.by, 0.0)
    
    @property
    def vc(self):
        return XYZ(self.cx, self.cy, self.cz)
    
    @property
    def alpha(self):
        return angle(self.vb, self.vc)
        
    @property
    def beta(self):
        return acos(self.cx/self.c)
    
    @property
    def gamma(self):
        return atan2(self.by, self.bx)
    
    @property
    def volumn(self):
        return self.ax*self.by*self.cz
    
    def scale(self, v):
        try:
            self.x *= v.x
            self.y *= v.y
            self.z *= v.z
        except AttributeError:
            self.x *= v
            self.y *= v
            self.z *= v
            
    def norm(self, xyz):
        x = xyz.x
        y = xyz.y
        z = xyz.z
        while z < -self.cz/2:
            z += self.cz
            y += self.cy
            x += self.cx
        while z > self.cz/2:
            z -= self.cz
            y -= self.cy
            x -= self.cx
        while y < -self.by/2:
            y += self.by
            x += self.bx
        while y > self.by/2:
            y -= self.by
            x -= self.bx
        while x < -self.ax/2:
            x += self.ax
        while x > self.ax/2:
            x -= self.ax

        return XYZ(x, y, z)

    @staticmethod
    def parse(line):
        tokens = line.split()
        if len(tokens) >= 6:
            try:
                tokens = [float(i) for i in tokens[-6:]]
                tokens[3] = radians(tokens[3])
                tokens[4] = radians(tokens[4])
                tokens[5] = radians(tokens[5])
                return PBC(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5])
            except:
                pass
        return None
            
        
    def __str__(self):
        return "PBC: %12.4f %12.4f %12.4f %5.2f %5.2f %5.2f" % \
                (self.a, self.b, self.c, 
                 degrees(self.alpha), degrees(self.beta), degrees(self.gamma))
