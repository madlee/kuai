'''
Created on 2010-10-29

@author: Madlee
'''
import re
import math
default_prefix = {'E':1e18, 'P':1e15, 'T':1e12, 'G':1e9, 'M':1e6, 
                  'k':1e3, 'd':0.1, 'c':0.01,'m':0.001, 'u':1e-6, 
                  'n':1e-9, 'p':1e-12, 'f':1e-15, 'a':1e-18}

class UnitConverter:
    __pattern = re.compile(r'\d*$')
    
    def __init__(self, basic_unit, prefix=default_prefix):
        self.data = {}
        for i in basic_unit.split():
            self.data[i] = (1, [i], [])
        self.prefix = prefix
        
    def add_unit(self, key, value, unit):
        v, u, o = self.uniform(unit)
        v *= value
        self.data[key] = v, u, o
    
    def uniform(self, unit):
        result_o = []
        result_u = []
        result_value = 1
        tokens = unit.split('/', 1)
        over = tokens[0]
        for i in over.split():
            for j in UnitConverter.__split(i):
                if j not in self.data:
                    prefix = j[0]
                    if prefix in self.prefix and j[1:] in self.data:
                        j = j[1:]
                        result_value *= self.prefix[prefix]
                    else:
                        raise ValueError("Undefined unit %s" % j)
                        
                factor, o, u = self.data[j]
                result_value *= factor
                result_o += o
                result_u += u
                
        if len(tokens) > 1:
            assert len(tokens) == 2
            under = tokens[1]
            under = under.replace('/', ' ')
            v, o, u = self.uniform(under)
            result_value /= v
            result_u += o
            result_o += u
        
        temp_u = []
        for i in result_u:
            try:
                result_o.remove(i)
            except ValueError:
                temp_u.append(i)

        result_u = temp_u
        result_o.sort()
        result_u.sort()
        return result_value, result_o, result_u
    
    @staticmethod
    def __split(unit):
        if unit == '1':
            return []
        else:
            power = UnitConverter.__pattern.search(unit).group(0)
            if len(power) == 0:
                return [unit]
            else:
                unit = unit[0:-len(power)]
                return [unit] * int(power)
        
        
    def parse(self, value, unit):
        v, o, u = self.uniform(unit)
        return value * v, o, u
        
    def format(self, value, unit):
        v, o, u = self.uniform(unit)
        return value/v, o, u
        
    def convert(self, value, source, target):
        value, o, u = self.parse(value, source)
        value, o2, u2 = self.format(value, target)
        if o != o2 or u != u2:
            raise ValueError("%s is different to %s", source, target)
        else:
            return value, o, u
"""
default_units = UnitConverter('g m s e rad')
default_units.add_unit('A', 1e-10, 'm')
default_units.add_unit('N', 1, 'kg m/s2')
default_units.add_unit('mol', 6.02214151e23, '')
default_units.add_unit('J', 1, 'N m')
default_units.add_unit('cal', 4.184, 'J')
default_units.add_unit('Pa', 1, 'N / m2')
default_units.add_unit('deg', math.pi / 180, 'rad')
default_units.add_unit('amu', 1, 'g/mol')
default_units.add_unit('Hz', 1, '1/s')
default_units.add_unit('C', 1 / 1.60217653e-19, 'e')
default_units.add_unit('l', 0.001, 'm3')
default_units.add_unit('K', 2*1.38065e-23, 'J')
default_units.add_unit('atm', 101325, 'Pa')
"""

default_units = UnitConverter('amu A fs e rad')
default_units.add_unit('g', 6.02214151e23, 'amu')
default_units.add_unit('m', 1e10, 'A')
default_units.add_unit('s', 1e15, 'fs')
default_units.add_unit('mol', 6.02214151e23, '')
default_units.add_unit('N', 1, 'kg m/s2')
default_units.add_unit('J', 1, 'N m')
default_units.add_unit('cal', 4.184, 'J')
default_units.add_unit('Pa', 1, 'N / m2')
default_units.add_unit('deg', math.pi / 180, 'rad')
default_units.add_unit('Hz', 1, '1/s')
default_units.add_unit('C', 1 / 1.60217653e-19, 'e')
default_units.add_unit('l', 0.001, 'm3')
default_units.add_unit('K', 2*1.38065e-23, 'J')
default_units.add_unit('atm', 101325, 'Pa')

default_units.add_unit('hartree', 2625.49962, 'kJ/mol')
default_units.add_unit('bohr', 5.2917720859e-11, 'm')
default_units.add_unit('debye', 1.0 / 299792458 * 1e-21, 'C m')
