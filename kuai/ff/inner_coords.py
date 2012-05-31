'''
Created on 2010-3-25

@author: Madlee
'''

def list_bond_atoms(mol):
    result = []
    for i in mol.bonds:
        result.append(i.atom1)
        result.append(i.atom2)
    return result

def list_angle_atoms(mol):
    result = []
    for i in mol.atoms:
        n = mol.degree(i)
        for j in range(1, n):
            for k in range(0, j):
                ja = mol.neighbor_atom(i, j)
                ka = mol.neighbor_atom(i, k)
                if ja is not i and ka is not i:
                    result.append(ja)
                    result.append(i)
                    result.append(ka)
    return result


def list_torsion_atoms(mol):
    result = []
    for i in mol.bonds:
        d1 = mol.degree(i.atom1)
        for j1 in range(d1):
            ja1 = mol.neighbor_atom(i.atom1, j1)
            if ja1 is not i.atom1 and ja1 is not i.atom2:
                d2 = mol.degree(i.atom2)
                for j2 in range(d2):
                    ja2 = mol.neighbor_atom(i.atom2, j2)
                    if ja2 is not i.atom1 and ja2 is not i.atom2 and ja1 is not ja2:
                        result.append(ja1)
                        result.append(i.atom1)
                        result.append(i.atom2)
                        result.append(ja2)
    return result


def list_fork_0_atoms(mol):
    result = []
    for i in mol.atoms:
        n = mol.degree(i)
        if n == 3:
            result.append(i)
            result.append(mol.neighbor_atom(i, 0))
            result.append(mol.neighbor_atom(i, 1))
            result.append(mol.neighbor_atom(i, 2))
    return result

def list_fork_1_atoms(mol):
    result = []
    for i in mol.atoms:
        n = mol.degree(i)
        if n == 3:
            result.append(mol.neighbor_atom(i, 0))
            result.append(i)
            result.append(mol.neighbor_atom(i, 1))
            result.append(mol.neighbor_atom(i, 2))
    return result

def list_fork_2_atoms(mol):
    result = []
    for i in mol.atoms:
        n = mol.degree(i)
        if n == 3:
            result.append(mol.neighbor_atom(i, 0))
            result.append(mol.neighbor_atom(i, 1))
            result.append(i)
            result.append(mol.neighbor_atom(i, 2))
    return result

def list_fork_3_atoms(mol):
    result = []
    for i in mol.atoms:
        n = mol.degree(i)
        if n == 3:
            result.append(mol.neighbor_atom(i, 0))
            result.append(mol.neighbor_atom(i, 1))
            result.append(mol.neighbor_atom(i, 2))
            result.append(i)
    return result

def list_pair_13_atoms(mol):
    angles = list_angle_atoms(mol)
    pairs = set()
    for i in range(0, len(angles), 3):
        atom1 = angles[i]
        atom3 = angles[i+2]
        if mol.get_bond(atom1, atom3) != None:
            continue
        if id(atom1) < id(atom3):
            atom1, atom3 = atom3, atom1
            pairs.add((atom1, atom3))
    result = []
    for atom1, atom4 in pairs:
        result.append(atom1)
        result.append(atom4)
    return result

def list_pair_14_atoms(mol):
    torsions = list_torsion_atoms(mol)
    pairs = set()
    for i in range(0, len(torsions), 4):
        atom1 = torsions[i]
        atom4 = torsions[i+3]
        if mol.get_bond(atom1, atom4) != None:
            continue
        if mol.get_angle(atom1, atom4) != None:
            continue
        if id(atom1) < id(atom4):
            atom1, atom4 = atom4, atom1
        pairI = (atom1, atom4)
        pairs.add(pairI)
    result = []
    for atom1, atom4 in pairs:
        result.append(atom1)
        result.append(atom4)
    return result




