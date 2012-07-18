from math import cos, sin, acos, atan2, sqrt, pi, degrees, radians
from kuai.xyz import XYZ, angle
from kuai.mol.element import Element 
from kuai.kuaisuan import *
from kuai.sim.model import setup_model
from kuai.sim.efunc import make_efunc
from kuai.unit import default_units

from numpy import sum, dot

class PBC:
    def __init__(self, a, b, c, alpha=pi/2, beta=pi/2, gamma=pi/2):
        cosA, cosB, cosG = cos(alpha), cos(beta), cos(gamma)
        sinG = sin(gamma)
        self.ax = a
        self.bx,  self.by = b*cosG, b*sinG
        self.cx, self.cy, = c*cosB, c*(cosA-cosB*cosG)/sinG 
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


class BasicEnergyFunction:
    def start(self, model):
        pass
    
    def finish(self, model):
        pass
    
    def finalize(self, model):
        pass
    
    def getE(self, model):
        raise NotImplementedError()
    
    def getG(self, model):
        raise NotImplementedError()
    
    def getEG(self, model):
        return self.getE(model). self.getG(model)
    

class BasicSimulationControl:
    def start(self, model):
        pass
    
    def finish(self, model, nstep):
        pass

    def finalize(self, model):
        pass
    
    def act(self, model, step):
        raise NotImplementedError()
    
class SimpleReporter(BasicSimulationControl):
    def __init__(self):
        self.format = '%(step)8s: %(potential_e)12s %(kinetic_e)12s %(temperature)8s %(pressure)8s %(density)8s: %(total_e)12s'
        self.e_unit = 'kJ/mol'
        self.t_unit = 'K'
        self.p_unit = 'atm'
        self.d_unit = 'kg/l'
        
        self.ewidth = '%12.4f'
        self.twidth = '%8.2f'
        self.pwidth = '%8.2f'
        self.dwidth = '%8.3f'
        
    def start(self, model):
        self.e_factor, _, _ = default_units.format(1, self.e_unit)
        self.t_factor, _, _ = default_units.format(1, self.t_unit)
        self.p_factor, _, _ = default_units.format(1, self.p_unit)
        self.d_factor, _, _ = default_units.format(1, self.d_unit)
        
        ke2t = len(model.coords) 
        if hasattr(model, 'constrains'):
            ke2t -= model.constrains
        self.t_factor /= ke2t
        
        if hasattr(model, 'masses'):
            self.allmasses = sum(model.masses)
        else:
            self.allmasses = None
        
        data = {'step':'STEP', 'potential_e':'Potetial E',  'kinetic_e':'Kinetic E', 
                'temperature':'Temp', 'pressure':'Pres', 'density':'Dens',
                'total_e':'Total E'}
        line = self.format % data
        print '=' * len(line)
        print line
        data = {'step':'', 'potential_e':self.e_unit,  'kinetic_e':self.e_unit, 
                'temperature':self.t_unit, 'pressure':self.p_unit, 'density':self.d_unit,
                'total_e':self.e_unit}
        print self.format % data
        print '-' * len(line)
        print self.dump(model, 'START')
        
    
    def act(self, model, step):
        print self.dump(model, step)
        
    def dump(self, model, step):
        step = str(step)
        potential_e = sum(model.potential_e)
        total_e = potential_e
        potential_e = self.ewidth % (potential_e * self.e_factor) 
        
        if hasattr(model, 'kinetic_e'):
            kinetic_e = model.kinetic_e
            total_e += kinetic_e                    
            temperature = self.twidth % (kinetic_e * self.t_factor)
            kinetic_e = self.ewidth % (kinetic_e * self.e_factor)
        else:
            temperature = kinetic_e = 'N/A'

        total_e = self.ewidth % (total_e * self.e_factor)            
       
        if hasattr(model, 'pressure'):
            pressure = (model.pressure.x + model.pressure.y + model.pressure.z) / 3
            pressure = self.pwidth % (pressure * self.p_factor)
        else:
            pressure = 'N/A'
        
        if self.allmasses and hasattr(model, 'pbc'):
            density = self.allmasses / model.pbc.volumn 
            density = self.dwidth % (density * self.d_factor)
        else:
            density = 'N/A'

        line = self.format % vars()
        return line
        
    def finish(self, model, nstep):
        line = self.dump(model, 'FINISH')
        print '-'*len(line)
        print line
        print '='*len(line)

class BasicTemperatureControl(BasicSimulationControl):
    def start(self, model):
        self.vx = model.speeds[0::3]
        self.vy = model.speeds[1::3]
        self.vz = model.speeds[2::3]
        self.act(model, 0)
    
    def finish(self, model, nstep):
        pass

    def finalize(self, model):
        pass
    
    def act(self, model, step):
        ex = dot(self.vx * self.vx, model.masses)
        ey = dot(self.vy * self.vy, model.masses)
        ez = dot(self.vz * self.vz, model.masses)
        model.kinetic_e = (ex + ey + ez) / 2

class BasicPressureControl(BasicSimulationControl):
    def start(self, model):
        self.x = model.coords[0::3]
        self.y = model.coords[1::3]
        self.z = model.coords[2::3]
        self.fx = model.grediants[0::3]
        self.fy = model.grediants[1::3]
        self.fz = model.grediants[2::3]
        self.act(model, 0)
    
    def finish(self, model, nstep):
        pass

    def finalize(self, model):
        pass
    
    def act(self, model, step):
        vx = model.kinetic_e * 2 + dot(self.fx, self.x)
        vy = model.kinetic_e * 2 + dot(self.fy, self.y)
        vz = model.kinetic_e * 2 + dot(self.fz, self.z)
        factor = 1.0 / model.pbc.volumn
        vx *= factor
        vy *= factor
        vz *= factor
        
        model.pressure = XYZ(vx, vy, vz)
