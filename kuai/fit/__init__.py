from numpy import array, matrix
from scipy.linalg import lstsq as linear_lstsq
from scipy.optimize import leastsq as lm_lstsq

class FittingFunction:
    def get_y0(self):
        raise NotImplementedError()
    
    def get_k0(self):
        raise NotImplementedError()
    
    def calc_y(self, k):
        raise NotImplementedError()
    
    def calc_dydk(self, k):
        return None
   
    def numerical_dydk(self, k, delta = None, rank = 1):
        if delta is None:
            delta = [0.1] * len(k)
        else:
            assert len(delta) == len(k)
            delta = delta[:]

        result = None
        for i in range(rank):
            resultI = []
            for j in range(len(delta)):
                k[j] += delta[j]
                y1 = self.calc_y(k)
                k[j] -= 2 * delta[j]
                y2 = self.calc_y(k)
                dYdKj = (array(y2) - array(y1)) / (2 * delta[j])
                resultI.append(dYdKj)
                k[j] += delta[j]
                delta[j] *= 0.5
            if result is None:
                result = array(resultI)
            else:
                resultI = array(resultI)*(4**i) - result
                result = resultI / (4**i - 1)
        return matrix(result).T
                
            
def linear_fit(func, full_output=True, cond=None):
    k0 = func.get_k0()
    y0 = func.get_y0()
    dydk = func.calc_dydk(k0)
    if dydk is None:
        dydk = func.numerical_dydk(k0)
    result = linear_lstsq(dydk, array(y0) - array(func.calc_y([0] * len(k0))), cond, True, True)
    if full_output:
        return result
    else:
        return result[0]
    
class __LevenbergMarquardtFittingFunction:
    def __init__(self, func):
        self.func = func
        self.y0 = array(func.get_y0())
        
    def __call__(self, k):
        y1 = self.func.calc_y(k)
        return array(y1) - self.y0
        
def lm_fit(func, **kwargs):
    k0 = func.get_k0()
    dydk = func.calc_dydk(k0)
    if dydk is None:
        dfunc = None
    else:
        dfunc = __LevenbergMarquardtFittingFunction(func)
    
    ffunc = __LevenbergMarquardtFittingFunction(func)
    
    return lm_lstsq(ffunc, k0, DFunc=dfunc, **kwargs)
    
    
    
class ForceFieldFittingFunction(FittingFunction):
    def __init__(self, index, parameters):
        self.index = index
        self.parameters = parameters
        self.k2par = []
        self.par2k = {}
        
    def calc_y(self, k):
        self.set_k(k)
        return self.get_y()
    
    def calc_dydk(self, k):
        self.set_k(k)
        return self.get_dydk()
   
    def get_y(self):
        raise NotImplementedError()
   
    def get_dydk(self):
        raise NotImplementedError()
    
    def set_k(self, k):
        assert len(k) == len(self.k2par)
        for i in len(k):
            try:
                self.parameters[self.k2par[i]] = k[i]
            except KeyError:
                for j in self.k2par[i]:
                    self.parameters[j] = k[i]
    
    def get_k0(self):
        result = []
        for i in self.k2par:
            try:
                result.append(self.parameters[i])
            except KeyError:
                result.append(self.parameters[i[0]])
        return result
    
