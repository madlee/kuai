from kuai.kuaisuan import *
from kuai.sim.model import setup_model
from kuai.sim.efunc import make_efunc
from kuai.unit import default_units

from numpy import sum, dot


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
