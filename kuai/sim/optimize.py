'''
Created on 2010-11-23

@author: Madlee
'''

from kuai.kuaisuan import SimulationJob
from scipy import optimize

class __EnergyFunction:
    def __init__(self, model, funcs, job):
        self.model = model
        self.funcs = funcs
        self.job = job

    def __call__(self, x):
        self.model.coords = x
        e = self.job.getE(self.model, self.funcs)
        result = sum(e)
        return result

        
class __GradientFunction:
    def __init__(self, model, funcs, job):
        self.model = model
        self.funcs = funcs
        self.job = job
        
    def __call__(self, x):
        self.model.coords = x
        _, g = self.job.getEG(self.model, self.funcs)
        return g
    
class __HessianFunction:
    def __init__(self, model, funcs, job):
        self.model = model
        self.funcs = funcs
        self.job = job
        
    def __call__(self, x):
        self.model.coords = x
        _, _, h = self.job.getEGH(self.model, self.funcs)
        return h

class OptimizeJob(SimulationJob):
    def efunc(self, model, funcs):
        return __EnergyFunction(model, funcs, self)
    
    def gfunc(self, model, funcs):
        return __GradientFunction(model, funcs, self)
    
    def hfunc(self, model, funcs):
        return __HessianFunction(model, funcs, self)
        
    def fmin_simplex(self, model, funcs, *args, **kwargs):
        efunc = self.efunc(model, funcs)
        result = optimize.fmin(efunc, model.coords, *args, **kwargs)
        model.coords = result[0]
        return result
        
    def fmin_powell(self, model, funcs, *args, **kwargs):
        efunc = self.efunc(model, funcs)
        result = optimize.fmin_powell(efunc, model.coords, *args, **kwargs)
        model.coords = result[0]
        return result
        
    def fmin_cg(self, model, funcs, *args, **kwargs):
        efunc = self.efunc(model, funcs)
        gfunc = self.gfunc(model, funcs)
        result = optimize.fmin_cg(efunc, model.coords, gfunc, *args, **kwargs)
        model.coords = result[0]
        return result
        
    def fmin_bfgs(self, model, funcs, *args, **kwargs):
        efunc = self.efunc(model, funcs)
        gfunc = self.gfunc(model, funcs)
        result = optimize.fmin_bfgs(efunc, model.coords, gfunc, *args, **kwargs)
        model.coords = result[0]
        return result
    
    def fmin_ncg(self, model, funcs, *args, **kwargs):
        efunc = self.efunc(model, funcs)
        gfunc = self.gfunc(model, funcs)
        hfunc = self.hfunc(model, funcs)
        result = optimize.fmin_ncg(efunc, model.coords, gfunc, fhess = hfunc, *args, **kwargs)
        model.coords = result[0]
        return result
