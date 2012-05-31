'''
Created on 2010-10-29

@author: Madlee
'''

from kuai.unit import default_units
from kuai.ff.combine import *

def remove_star(par):
    tokens = [i.strip() for i in par.split(',')]
    for i in range(len(tokens)):
        if tokens[i].endswith('*'):
            tokens[i] = tokens[i][:-1]
    tokens = [float(i) for i in tokens]
    return tokens
   
def parse_line(line):
    tokens = [i.strip() for i in line.split(':')]
    if len(tokens) >= 3:
        tokens[1] = [i.strip() for i in tokens[1].split(',')]
        tokens[2] = remove_star(tokens[2])
        return tokens
    elif len(tokens) == 0:
        return None
    else:
        raise ValueError('Can not parse line %s' % line)

def read_ppf(file):
    line = file.next().rstrip()
    if line != '#DFF:PPF':
        raise RuntimeError("The input file is not a valid PPF file.")
    line = file.next().rstrip()
    if not line.startswith('#PROTOCOL'):
        raise RuntimeError("The input file is not a valid PPF file.")
    tokens = line.split('=')
    if len(tokens) != 2:
        raise RuntimeError("The input file is not a valid PPF file.")
    
    type = tokens[1].strip()
         
    COULOMB_CONST, _, _ = default_units.parse(8.987551787e9, 'N m2 / C2')
    index = {}
    parameters = [0.0, COULOMB_CONST]

    index['forcefield']={'type':type}

    for line in file:
        tokens = parse_line(line)
        if tokens:
            if tokens[0] == 'E0':
                # Useless term
                pass
                
            elif tokens[0] == 'ATYPE':
                if 'atom/mass' not in index:
                    index['atom/mass'] = {}
                key2par = index['atom/mass']

                key = ' '.join(tokens[1])
                tokens = tokens[2]
                assert len(tokens) == 2
                par = tokens[1]
                par, _, _ = default_units.parse(par, 'amu')
                key2par[key] = len(parameters)
                parameters.append(par)
                
            elif tokens[0] == 'N12_6':
                if 'nonbond/lj12_6' not in index:
                    index['nonbond/lj12_6'] = {}
                key2par = index['nonbond/lj12_6']

                key = ' '.join(tokens[1])
                tokens = tokens[2]
                assert len(tokens) == 2
                key2par[key] = len(parameters)

                par, _, _ = default_units.parse(tokens[0], 'A')
                parameters.append(par)
                par, _, _ = default_units.parse(tokens[1], 'kcal/mol')
                parameters.append(par)
                
            elif tokens[0] == 'ATC':
                if 'atom/charge' not in index:
                    index['atom/charge'] = {}
                key2par = index['atom/charge']
                key = ' '.join(tokens[1])
                tokens = tokens[2]
                assert len(tokens) == 1
                key2par[key] = len(parameters)
                
                par, _, _ = default_units.parse(tokens[0], 'e')
                parameters.append(par)
                
            elif tokens[0] == 'BINC':
                if 'bond/binc' not in index:
                    index['bond/binc'] = {}
                key2par = index['bond/binc']
                
                key = tokens[1]
                assert len(key) == 2
                if key[0] > key[1]:
                    key[0], key[1] = key[1], key[0]
                    reverse=True
                else:
                    reverse=False
                    
                key = ' '.join(key)

                tokens = tokens[2]
                assert len(tokens) == 1
                key2par[key] = len(parameters)
                par, _, _ = default_units.parse(tokens[0], 'e')
                if reverse:
                    parameters.append(-par)
                else:
                    parameters.append(par)

            elif tokens[0] == 'BHARM':
                if 'bond/harmonic' not in index:
                    index['bond/harmonic'] = {}
                key2par = index['bond/harmonic']
                
                key = tokens[1]           
                assert len(key) == 2
                if key[0] > key[1]:
                    key[0], key[1] = key[1], key[0]
                    reverse=True
                else:
                    reverse=False
                    
                key = ' '.join(key)

                tokens = tokens[2]
                assert len(tokens) == 2
                key2par[key] = len(parameters)
                
                par, _, _ = default_units.parse(tokens[0], 'A')
                parameters.append(par)
                par, _, _ = default_units.parse(tokens[1], 'kcal/mol/A2')
                parameters.append(par)

            elif tokens[0] == 'AHARM':
                if 'angle/harmonic' not in index:
                    index['angle/harmonic'] = {}
                key2par = index['angle/harmonic']
                
                key = tokens[1]
                assert len(key) == 3
                if key[0] > key[2]:
                    key[0], key[2] = key[2], key[0]
                    reverse=True
                else:
                    reverse=False
                    
                key = ' '.join([i.strip() for i in key])

                tokens = tokens[2]
                assert len(tokens) == 2
                key2par[key] = len(parameters)
                par, _, _ = default_units.parse(tokens[0], 'deg')
                parameters.append(par)
                par, _, _ = default_units.parse(tokens[1], 'kcal/mol/rad2')
                parameters.append(par)
                
            elif tokens[0] == 'TCOSP':
                if 'torsion/cos_poly5' not in index:
                    index['torsion/cos_poly5'] = {}
                key2par = index['torsion/cos_poly5']
                
                key = tokens[1]
                assert len(key) == 4
                if key[1] > key[2] or key[1] == key[2] and key[0] > key[3]:
                    key[0], key[3] = key[3], key[0]
                    key[1], key[2] = key[2], key[1]
                    reverse=True
                else:
                    reverse=False

                key = ' '.join([i.strip() for i in key])

                tokens = tokens[2]
                assert len(tokens) % 3 == 0
                key2par[key] = len(parameters)
                
                v = [0] * 5
                k = [0] * 6
                for i in range(0, len(tokens), 3):
                    n = int(tokens[i+2]+.1)
                    kI = tokens[i+1]
                    phi = tokens[i+0]
                    if phi > 150:
                        assert abs(phi-180) < 1
                        kI = -kI
                    else:
                        assert abs(phi) < 1
                    if n % 2 == 0:
                        kI = -kI
                    assert v[n] == 0
                    v[n] = 2 * kI
                
                k[0]=v[2]+(v[1]+v[3])/2
                k[1]=0.5*(3*v[3]-v[1])
                k[2]=-v[2]+4*v[4]
                k[3]=-2*v[3]
                k[4]=-4*v[4]
                k[5]=0
                for i in range(0, len(k)):
                    k[i], _, _ = default_units.parse(k[i], 'kcal/mol')
                parameters += k

            elif tokens[0] == 'IBCOS':
                if 'improper_torsion/cos_poly5' not in index:
                    index['improper_torsion/cos_poly5'] = {}
                key2par = index['improper_torsion/cos_poly5']
                
                key = tokens[1]
                assert len(key) == 4
                key = ' '.join(key)

                tokens = tokens[2]
                assert len(tokens) % 3 == 0
                key2par[key] = len(parameters)
                
                v = [0] * 6
                k = [0] * 6
                for i in range(0, len(tokens), 3):
                    n = int(tokens[i+2]+.1)
                    kI = tokens[i+1]
                    phi = tokens[i+0]
                    if phi > 150:
                        assert abs(phi-180) < 1
                        kI = -kI
                    else:
                        assert abs(phi) < 1
                    if n % 2 == 0:
                        kI = -kI
                    assert v[n] == 0
                    v[n] = 2 * kI
                
                k[0]=v[2]+(v[1]+v[3])/2
                k[1]=0.5*(3*v[3]-v[1])
                k[2]=-v[2]+4*v[4]
                k[3]=-2*v[3]
                k[4]=-4*v[4]
                k[5]=0
                for i in range(0, len(k)):
                    k[i], _, _ = default_units.parse(k[i], 'kcal/mol')
                parameters += k
                
            else:
                raise RuntimeError("The function %s is not supported yet." % tokens[0])

    if type == 'AMBER':
        index['forcefield']['scale14/vdw'] = 0.5
        index['forcefield']['scale14/columnb'] = 1 / 1.2
        atomtype, newindex = combine(index['nonbond/lj12_6'], parameters, combine_arithmatic)
        newindex[':columnb-constant:'] = 1
        index['pair-nonbond/lj12_6'] = newindex
        index['atom/type'] = atomtype
        
        index['pair-14/lj12_6'] = scale_vdw(newindex, parameters, 0.5, 1)
        columnb14 = parameters[newindex[':columnb-constant:']] / 1.2
        index['pair-14/lj12_6'][':columnb-constant:'] = len(parameters)
        parameters.append(columnb14)
    else:
        raise RuntimeError("The force field type %s is not supported yet." % type)
        
    return index, parameters

def read_ppffile(filename):
    with open(filename) as file:
        return read_ppf(file)
