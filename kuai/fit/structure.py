'''
Created on 2011-1-24

@author: madlee
'''

from math import sqrt
from numpy.random import rand
from kuai.xyz import XYZ, angle
from kuai.fit import ForceFieldFittingFunction
from kuai.ff.inner_coords import list_bond_atoms
from kuai.ff.inner_coords import list_angle_atoms
from kuai.sim.optimize import OptimizeJob
from kuai.sim.efunc import make_efunc
from kuai.sim.model import setup_model

class FittingStructureFunction(ForceFieldFittingFunction):
    def __init__(self):
        self.mols = []
        self.ny = 0
        self.optjob = OptimizeJob()
    
    def get_y0(self):
        return [0] * self.ny
    
    def calc_y(self):
        result = []
        for model, bonds, angles, efunc, y0 in self.mols:
            coords = model.coords
            model.coords = 2 * self.noise * rand(len(coords)) - self.noise
            model.coords += coords
            self.optjob.fmin_ncg(model, efunc)
            y1 = FittingStructureFunction.__calc_ab(model, bonds,angles)
            assert len(y0) == len(y1)
            for i in range(len(y0)):
                result.append((y0[i]-y1[i])/y0[i])
            model.coords = coords   # restore old coords
            
        raise NotImplementedError()

    @staticmethod
    def __calc_ab(model, bonds, angles):
        y0 = []
        for i1 in range(0, len(bonds), 2):
            i2 = i1 + 1
            x = model.coords[i1*3] - model.coords[i2*3]
            y = model.coords[i1*3+1] - model.coords[i2*3+1]
            z = model.coords[i1*3+2] - model.coords[i2*3+2]
            y0.append(sqrt(x*x+y*y+z*z))

        for i1 in range(0, len(angles), 3):
            i2 = i1 + 1
            i3 = i2 + 2
            v1 = XYZ(model.coords[i1*3], model.coords[i1*3+1], model.coords[i1*3+2])
            v2 = XYZ(model.coords[i2*3], model.coords[i2*3+1], model.coords[i2*3+2])
            v3 = XYZ(model.coords[i3*3], model.coords[i3*3+1], model.coords[i3*3+2]) 
            y0.append(angle(v1, v2, v3))
        
        return y0
        
    def add_mol(self, mol):
        bonds = list_bond_atoms(mol)
        bonds = [mol.index(i) for i in bonds]
        angles = list_angle_atoms(mol)
        bonds = [mol.index(i) for i in bonds]
        
        efunc = make_efunc(mol, self.index)
        model = setup_model(mol, self.index, self.par)
        
        y0 = FittingStructureFunction.__calc_ab(model, bonds, angles)
        self.mols.append((model, bonds, angles, efunc, y0))
        self.ny += len(bonds) / 2 + len(angles) / 3
