'''
Created on 2011-1-23

@author: madlee
'''

from kuai.sim import SimulationJob
from kuai.constants import SPEED_OF_LIGHT, PI
from numpy import diag, reshape, matrix

def frequency(model, efuncs, job=None, masses = None):
    if job is None:
        job = SimulationJob()
    _, _, h = job.getEGH(model, efuncs)
    ndim = len(model.coords)
    h = matrix*(reshape(h, (ndim, ndim)))
    
    if masses is None:
        masses = []
        for i in model.masses:
            masses += [i] * 3
        masses = diag(masses)
        
    h = masses * h * masses
    freq = eigvals(h)
    freq *= 1/(2 * PI * SPEED_OF_LIGHT)
    return freq