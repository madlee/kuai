'''
Created on 2010-6-8

@author: Madlee
'''

from sqlite3 import connect
from os.path import isfile
from cStringIO import StringIO
from kuai import mol
from kuai.mol.algo import guess_bond_link, split
from kuai.mol.io import format_mol
from hashlib import md5
from time import clock
#from tinymol.tool import inchi

TRAJEC_DB_NAME = "trajec.db"

SQL_SET_CREATE_FRAMES_TABLE = ["""CREATE TABLE frames (
    frame INTEGER PRIMARY KEY,
    data TEXT
);
"""]

SQL_SET_CREATE_ANALYSIS_TABLE = [
"""CREATE TABLE moltypes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    natoms INTEGER,
    nbonds INTEGER,
    formula CHAR(50),
    inchi_key CHAR(50) UNIQUE,
    data TEXT DEFAULT ''
);
""",

"""CREATE TABLE molecules (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    moltype INTEGER REFERENCES moltypes (id),
    atom_hash CHAR(50),
    bond_hash CHAR(50),
    data TEXT DEFAULT '',
    UNIQUE (moltype, atom_hash, bond_hash)
);
""",

"""CREATE TABLE mol_in_frame (
    frame INTEGER REFERENCES frames (id),
    molecule INTEGER REFERENCES molecules (id),
    PRIMARY KEY(frame, molecule)
);
""",

"""CREATE TABLE atom_in_mol (
    atom INTEGER,
    molecule INTEGER REFERENCES molecules (id),
    PRIMARY KEY(atom, molecule)
);
""",
"""CREATE TABLE bond_in_mol (
    atom1 INTEGER,
    atom2 INTEGER,
    molecule INTEGER REFERENCES molecules (id),
    PRIMARY KEY(atom1, atom2, molecule)
);
""",
"""CREATE VIEW mol_info_view AS 
    SELECT molecules.id as id, formula, inchi_key, moltype,
        moltypes.data AS inchi_code, atom_hash, bond_hash 
    FROM molecules INNER JOIN moltypes ON molecules.moltype = moltypes.id
""",
"""CREATE VIEW frame_statics AS 
    SELECT frame, moltype, formula, inchi_key, inchi_code, count(*) as hits
    FROM mol_info_view INNER JOIN mol_in_frame ON mol_info_view.id = mol_in_frame.molecule 
    GROUP BY frame, moltype
    ORDER BY frame, moltype""",
]

SQL_CREATE_TABLE_MOL_LIFE = """CREATE TABLE mol_life (
    molecule INTEGER REFERENCES molecules (id),
    birth INTEGER REFERENCES frames (frame),
    death INTEGER,
    PRIMARY KEY(molecule, birth)
);
"""
SQL_CREATE_TABLE_BOND_LIFE = """CREATE TABLE bond_life (
    atom1 INTEGER,
    atom2 INTEGER,
    birth INTEGER REFERENCES frames (frame),
    death INTEGER,
    PRIMARY KEY(atom1, atom2, birth)
);
"""

SQL_INSERT_FRAME = """INSERT INTO frames (frame, data) VALUES (?, ?)"""

SQL_SELECT_ALL_FRAME = """SELECT frame, data FROM frames;"""

SQL_FIND_MOLTYPE = """SELECT id FROM moltypes WHERE inchi_key = ?;"""

SQL_INSERT_MOLTYPE = """INSERT INTO moltypes (natoms, nbonds, formula, inchi_key, data) VALUES (?, ?, ?, ?, ?);"""

SQL_FIND_MOLECULE = """SELECT id FROM molecules WHERE moltype = ? AND atom_hash = ? AND bond_hash = ?;"""

SQL_INSERT_MOLECULE = """INSERT INTO molecules (moltype, atom_hash, bond_hash) VALUES (?, ?, ?)"""

SQL_INSERT_MOL_IN_FRAME = """INSERT INTO mol_in_frame (frame, molecule) VALUES (?, ?);"""

SQL_INSERT_ATOM_IN_MOL = """INSERT INTO atom_in_mol (atom, molecule) VALUES (?, ?);"""
 
SQL_INSERT_BOND_IN_MOL = """INSERT INTO bond_in_mol (atom1, atom2, molecule) VALUES (?, ?, ?);"""

SQL_TRACE_MOL_IN_FRAME = """SELECT frame
    FROM mol_in_frame 
    WHERE molecule = ? 
    ORDER BY frame"""

SQL_INSERT_MOL_LIFE = """INSERT INTO mol_life (molecule, birth, death) VALUES (?, ?, ?)"""

SQL_TRACE_ATOM_IN_FRAME = """SELECT frame, molecule 
    FROM mol_in_frame 
    WHERE molecule IN (SELECT molecule FROM atom_in_mol WHERE atom = ?)
    ORDER BY frame"""
    
SQL_TRACE_BOND_IN_FRAME = """SELECT frame, molecule
    FROM mol_in_frame 
    WHERE molecule IN (SELECT molecule FROM bond_in_mol WHERE atom1 = ? AND atom2 = ?)
    ORDER BY frame"""

SQL_SELECT_ALL_MOLID = """SELECT id FROM molecules ORDER BY id"""

SQL_SELECT_TOP2_FRAME = """SELECT frame  
    FROM frames
    ORDER BY frame LIMIT 2"""
    
SQL_PEEK_MOL = """SELECT * FROM mol_info_view WHERE id = ?"""
    
SQL_SELECT_ALL_BOND = """SELECT DISTINCT atom1, atom2 FROM bond_in_mol"""

SQL_INSERT_BOND_LIFE = """INSERT INTO bond_life (atom1, atom2, birth, death) VALUES (?, ?, ?, ?)"""

SQL_PEEK_MOL_ATOM_AT_FRAME = """SELECT id, formula, inchi_key, moltype, inchi_code
    FROM mol_info_view 
    WHERE id in (
        SELECT DISTINCT mol_in_frame.molecule 
            FROM mol_in_frame INNER JOIN atom_in_mol ON mol_in_frame.molecule = atom_in_mol.molecule 
            WHERE frame = ? AND atom in (
                SELECT atom FROM atom_in_mol WHERE molecule=?)) 
        ORDER BY id"""
        
SQL_SELECT_BIRTH_OF_TYPE = """SELECT id, birth 
    FROM mol_info_view INNER JOIN mol_life ON mol_info_view.id == mol_life.molecule 
    WHERE moltype = ?"""

def create_db(filename=TRAJEC_DB_NAME):
    if isfile(filename):
        raise RuntimeError("%s is existed in current folder." % filename)

    conn = connect(filename)
    cursor = conn.cursor()

    data = globals()
    for i in SQL_SET_CREATE_FRAMES_TABLE:
        cursor.execute(i % data)
    conn.commit()
    del cursor
    return conn

def open_db(filename=TRAJEC_DB_NAME):
    conn = connect(filename)
    cursor = conn.cursor()
    return conn, cursor

def save_frame(cursor, frame, data):
    cursor.execute(SQL_INSERT_FRAME, (frame, data))

def parse_mol(data, vdw_radius=None, tor=0.2):
    file = StringIO(data)
    i = file.next()
    natom = int(i.strip())
    pbc = mol.PBC.parse(file.next())
    atoms = []
    for i in file:
        tokens = i.split()
        assert len(tokens) >= 4
        a = mol.Atom(tokens[0])
        a.coord.x = float(tokens[1])
        a.coord.y = float(tokens[2])
        a.coord.z = float(tokens[3])
        atoms.append(a)

    assert len(atoms) == natom
    bonds = guess_bond_link(atoms, pbc, vdw_radius, tor)
    result = mol.Molecule(atoms, bonds)
    if pbc:
        result.pbc = pbc
    return result

    
def formula(mol):
    count = {}
    for i in mol.atoms:
        if i.symbol in count:
            count[i.symbol] += 1
        else:
            count[i.symbol] = 1
    result = []
    if 'C' in count:
        result.append('C')
        if count['C'] > 1:
            result.append(str(count['C']))
        del count['C']
    if 'H' in count:
        result.append('H')
        if count['H'] > 1:
            result.append(str(count['H']))
        del count['H']
    keys = [i for i in count.iterkeys()]
    keys.sort()
    for i in keys:
        result.append(i)
        if count[i] > 1:
            result.append(str(count[i]))
    return ''.join(result)

def hash_atom(atoms):
    atoms.sort()
    atomset = [str(i) for i in atoms]
    return md5(",".join(atomset)).hexdigest()

def hash_bond(bonds):
    bonds.sort()
    atomset = [str(i) + '-' + str(j) for i , j in bonds]
    return md5(",".join(atomset)).hexdigest()

def save_mol_type(cursor, natoms, nbonds, formula, key, code):
    try:
        cursor.execute(SQL_INSERT_MOLTYPE, (natoms, nbonds, formula, key, code))
    except:
        pass
    cursor.execute(SQL_FIND_MOLTYPE, (key, ))
    result = cursor.fetchone()
    return result[0]

def save_molecule(cursor, moltype, atoms, bonds):
    hasha = hash_atom(atoms)
    hashb = hash_bond(bonds)
    try:
        cursor.execute(SQL_INSERT_MOLECULE, (moltype, hasha, hashb))
        inserted = True
    except:
        inserted = False
    cursor.execute(SQL_FIND_MOLECULE, (moltype, hasha, hashb))
    result = cursor.fetchone()[0]
    if inserted:
        for i in atoms:
            cursor.execute(SQL_INSERT_ATOM_IN_MOL, (i+1, result))
        for i in bonds:
            cursor.execute(SQL_INSERT_BOND_IN_MOL, (i[0]+1, i[1]+1, result))
    return result

def save_mol_set(cursor, frame, data, vdw_radius=None, tor=0.2):
    mol = parse_mol(data, vdw_radius, tor)
    molset = split(mol)
    for molI in molset:
        f = formula(molI)
        s = format_mol(molI, 'mol')
        code, _, k = inchi(s)
        moltype = save_mol_type(cursor, len(molI.atoms), len(molI.bonds), f, k, code)
        atomid = [mol.atoms.index(i) for i in molI.atoms]

        bondid = []
        for bond in molI.bonds:
            atom1 = mol.atoms.index(bond.atom1)
            atom2 = mol.atoms.index(bond.atom2)
            if atom1 > atom2:
                atom1, atom2 = atom2, atom1
            bondid.append((atom1, atom2))
        
        molid = save_molecule(cursor, moltype, atomid, bondid)
        cursor.execute(SQL_INSERT_MOL_IN_FRAME, (frame, molid))
        
def trace_mol(cursor, molid, gap):
    """generator of (birth, death, life_time)"""
    cursor.execute(SQL_TRACE_MOL_IN_FRAME, (molid,))
    data = cursor.fetchall()
    i = 0
    start = finish = data[i][0]
    hasNext = True
    while True:
        hasNext = False
        while i < len(data)-1:
            i += 1
            finish += gap
            if finish != data[i][0]:
                hasNext = True
                break
        yield (start, finish, finish-start)
        if hasNext:
            start = finish = data[i][0]
        else:
            break
    # raise StopIteration()

def trace_atom(cursor, atomid):
    """generator (molecule, birth, death, life_time)"""
    cursor.execute(SQL_TRACE_ATOM_IN_FRAME, (atomid,))
    data = cursor.fetchall()
    i = 0
    start = data[i][0]
    molid = data[i][1]
    hasNext = True
    while hasNext:
        hasNext = False
        while i < len(data)-1:
            i += 1
            if molid != data[i][1]:
                hasNext = True
                break
        yield molid, start, data[i][0], data[i][0] - start
        if hasNext:
            start = data[i][0]
            molid = data[i][1]
        else:
            break


def trace_bond(cursor, atom1, atom2, gap):
    """generator (molecule, birth, death, life_time) for bond (atom1 atom2)"""
    if atom1 > atom2:
        atom1, atom2 = atom2, atom1
    cursor.execute(SQL_TRACE_BOND_IN_FRAME, (atom1, atom2))
    data = cursor.fetchall()
    i = 0
    start = finish = data[i][0]
    molid = data[i][1]
    hasNext = True
    while hasNext:
        hasNext = False
        while i < len(data)-1:
            i += 1
            finish += gap
            if finish != data[i][0]:
                hasNext = True
                break
        yield molid, start, finish, finish - start
        if hasNext:
            start = finish = data[i][0]
            molid = data[i][1]
        else:
            break

def get_time_gap(cursor):
    cursor.execute(SQL_SELECT_TOP2_FRAME)
    rows = cursor.fetchall()
    return rows[1][0] - rows[0][0]
         
    
def parse_trajec_db(cursor, vdw_radius=None, tor=0.2):
    for i in SQL_SET_CREATE_ANALYSIS_TABLE:
        cursor.execute(i)

    cursor.execute(SQL_SELECT_ALL_FRAME)
    allframe = cursor.fetchall()
    for frame, data in allframe:
        data = data.encode('ascii')
        t1 = clock()
        save_mol_set(cursor, frame, data, vdw_radius, tor)
        print frame, clock() - t1

def save_mol_life_time(cursor, gap=None):
    cursor.execute(SQL_CREATE_TABLE_MOL_LIFE)
    if gap is None:
        gap = get_time_gap(cursor)
    cursor.execute(SQL_SELECT_ALL_MOLID)
    allmolid = cursor.fetchall()
    for id, in allmolid:
        t1 = clock()
        for i, j, _ in trace_mol(cursor, id, gap):
            # print id, i, j,
            cursor.execute(SQL_INSERT_MOL_LIFE, (id, i, j,))
        print "Save mol", id, 'life time cost', clock() - t1, '"'


def save_bond_life_time(cursor, gap=None):
    cursor.execute(SQL_CREATE_TABLE_BOND_LIFE)
    if gap is None:
        gap = get_time_gap(cursor)
    cursor.execute(SQL_SELECT_ALL_BOND)
    allbonds = cursor.fetchall()
    for atom1, atom2 in allbonds:
        assert atom1 < atom2
        t1 = clock()
        for molid, birth, death, _ in trace_bond(cursor, atom1, atom2, gap):
            print molid, atom1, atom2, birth, death
            cursor.execute(SQL_INSERT_BOND_LIFE, (atom1, atom2, birth, death,))
        print "Save Bond", atom1, atom2, 'life time cost', clock() - t1, '"'

def peek_mol(cursor, id):
    cursor.execute(SQL_PEEK_MOL, (id,))
    return cursor.fetchone()

def peek_mol_atom_in_frame(cursor, molid, frame):
    cursor.execute(SQL_PEEK_MOL_ATOM_AT_FRAME, molid, frame)
    return cursor.fetchall()

def count_moltype_precedes(cursor, moltype, gap=None):
    # TODO:
    if gap is None:
        gap = get_time_gap(cursor)
    result = {}
    cursor.execute(SQL_SELECT_BIRTH_OF_TYPE, (moltype,))
    for id, birth in cursor:
        row = peek_mol_atom_in_frame(cursor, id, birth-gap)
        
        for _, formula, inchi_key, moltype, inchi_code in row:
            pass
    return result
            
    


