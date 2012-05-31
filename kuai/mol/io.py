from kuai.iotool import Reader, ReaderManager, Writer, WriterManager
import kuai.mol as mol

class MoleculeReader(Reader):
    """The base class of all molecule reader"""
    pass

class MdlMolReader(MoleculeReader):
    """The reader to MDL Mol Files"""
    
    @staticmethod
    def parse_atom(line):
        x = float(line[0:10].strip())
        y = float(line[10:20].strip())
        z = float(line[20:30].strip())
        symbol = line[31:34].strip()
        if symbol == 'D':
            result = mol.Atom('H')
            result.coord = mol.XYZ(x, y, z)
            result.iso = 1
        elif symbol == 'T':
            result = mol.Atom('H')
            result.coord = mol.XYZ(x, y, z)
            result.iso = 2
        else:
            result = mol.Atom(symbol)
            result.coord = mol.XYZ(x, y, z)
        if len(line) > 36:
            iso = int(line[34:36].strip())
            if iso != 0:
                result.iso = iso
        result.charge = 0 
        if len(line) > 39:
            chg = int(line[36:39].strip())
            if chg != 0:
                result.charge = 4 - chg
        return result
    
    @staticmethod
    def parse_bond(line, atomlist):
        atom1 = atomlist[int(line[0:3])-1]
        atom2 = atomlist[int(line[3:6])-1]
        order = int(line[6:9])
        if order == 1:
            order = mol.Bond.SINGLE_BOND
        elif order == 2:
            order = mol.Bond.DOUBLE_BOND
        elif order == 3:
            order = mol.Bond.TRIPLE_BOND
        elif order == 4:
            order = mol.Bond.PARTIAL_BOND
        else:
            order = mol.Bond.UNKNOWN_BOND
        return mol.Bond(atom1, atom2, order)
    
    @staticmethod
    def parse_prop(mol, prop):
        return mol
    
    def read(self, file):
        """read molecule from a file like object. """
        name = file.next().rstrip()
        tag = file.next().rstrip()
        file.next()
        counts_line = file.next()
        try:
            natoms = int(counts_line[0:3].strip())
            nbonds = int(counts_line[3:6].strip())
            if len(counts_line) > 15:
                chiral_flag = int(counts_line[12:15].strip())
            else:
                chiral_flag = 0
            
            if len(counts_line) > 33:
                nproplines = counts_line[30:33].strip()
                if nproplines == '999':
                    nproplines = None
                else:
                    nproplines = int(nproplines)
            else:
                nproplines = 0
        except:
            raise RuntimeError("")

        atoms = []
        for i in range(natoms):
            atoms.append(MdlMolReader.parse_atom(file.next()))
        bonds = []
        for i in range(nbonds):
            bonds.append(MdlMolReader.parse_bond(file.next(), atoms))
        result = mol.Molecule(atoms, bonds)
        if name != '':
            result.name = name
        if tag != '':
            result.tag = tag
            
        if chiral_flag != 0:
            result.chiral_flag = chiral_flag
        
        prop = []
        if nproplines is None:
            for i in file:
                i = i.strip()
                if i == 'M  END':
                    break
                else:
                    prop.append(i)
        else:
            for i in range(nproplines):
                prop.append(file.next().strip())
        result = MdlMolReader.parse_prop(result, prop)        
        return result
    
        
class MdlSdfReader(MdlMolReader):
    """The reader to MDL SDF Files"""
    def read(self, file):
        result = MdlMolReader.read(self, file)
        result.data = {}
        tag = None
        data = []
        for i in file:
            if i.startswith('>  <') or i.startswith('$$$$'):
                if tag is not None:
                    result.data[tag] = ''.join(data).rstrip()
                data = []
                if i.startswith('$$$$'):
                    break
                else:
                    tag = i[4:i.find('>', 4)]
            else:
                i = i.rstrip()
                if len(i) == 81 and i[80] == '+':
                    i = i[:-1]
                    data.append(i)
                else:
                    data.append(i + '\n')
            
        return result
    
class DffMsdReader(MoleculeReader):
    def read(self, file):
        line = file.next()
        assert line.startswith('#DFF:MSD')
        line = file.next()
        assert line.startswith('#Model Structure Data File')

        line = file.next()
        pbc = None
        if line.startswith('PBC:'):
            pbc = mol.PBC.parse(line)
            line = file.next()
        natoms = int(line.strip())
        atoms = []
        for _ in range(natoms):
            atom = DffMsdReader.parse_atom(file.next())
            atoms.append(atom)
        nbonds = int(file.next().strip())
        bonds = []
        for _ in range(nbonds):
            bond = DffMsdReader.parse_bond(file.next(), atoms)
            bonds.append(bond)
        result = mol.Molecule(atoms, bonds)
        if pbc is not None:
            result.pbc = pbc
        
        for line in file:
            if line.startswith('#Formal Charge'):
                DffMsdReader.parse_formal_charge(file, result)
            elif line.startswith('M  END'):
                break
        return result
    
    @staticmethod
    def parse_atom(line):
        tokens = line.split()
        number = int(tokens[2])
        result = mol.Atom(number)
        result.name = tokens[1]
        result.type = tokens[3]
        if tokens[4] != '?':
            result.utype = tokens[4]
        result.partial_charge = float(tokens[5])
        result.coord.x = float(tokens[6])
        result.coord.y = float(tokens[7])
        result.coord.z = float(tokens[8])
        if len(tokens) > 9:
            result.mol_id = int(tokens[9])
        if len(tokens) > 10:
            if tokens[10] != 'UNK':
                result.residue = tokens[10]
        if len(tokens) > 11:
            group_id = int(tokens[11])
            result.group_id = group_id
        return result
    
    @staticmethod
    def parse_bond(line, atoms):
        tokens = line.split()
        a1 = int(tokens[0]) - 1
        a2 = int(tokens[1]) - 1
        order = int(tokens[2])
        if order == 1:
            order = mol.Bond.SINGLE_BOND
        elif order == 2:
            order = mol.Bond.DOUBLE_BOND
        elif order == 3:
            order = mol.Bond.TRIPLE_BOND
        elif order == -2:
            order = mol.Bond.PARTIAL_BOND
        else:
            order = mol.Bond.UNKNOWN_BOND
        return mol.Bond(atoms[a1], atoms[a2], order)
        
    @staticmethod
    def parse_formal_charge(file, mol):
        ncharged = int(file.next().strip())
        for _ in range(ncharged):
            tokens = file.next().split()
            index = int(tokens[0])-1
            charge = int(tokens[1])
            mol.atoms[index].charge = charge
        
    
MOL_READERS = ReaderManager()

MOL_READERS.readers['mdl'] = MOL_READERS.readers['mol'] = MdlMolReader()
MOL_READERS.readers['sdf'] = MdlSdfReader()
MOL_READERS.readers['msd'] = DffMsdReader()


def read_mol(file, type):
    return MOL_READERS.read(file, type)

def parse_mol(file, type):
    return MOL_READERS.parse(file, type)
    
def read_molfile(filename, type=None):
    return MOL_READERS.read_file(filename, type)


class MoleculeWriter(Writer):
    """The base class of all molecule writer"""
    pass


class MdlMolWriter(MoleculeWriter):
    def write(self, fileobj, var):
        fileobj.write('\n    Tiny Mol\n\n')
        fileobj.write('%3d%3d  0  0  0  0  0  0  0  0999 V2000\n' % (len(var.atoms), len(var.bonds)))
        for i in var.atoms:
            fileobj.write('%10.4f%10.4f%10.4f %-2s  0  0  0  0  0  0  0  0  0  0  0  0\n'
                          % (i.coord.x, i.coord.y, i.coord.z, i.symbol))
        for i in var.bonds:
            if i.order == mol.Bond.SINGLE_BOND:
                order = 1
            elif i.order == mol.Bond.DOUBLE_BOND:
                order = 2
            elif i.order == mol.Bond.TRIPLE_BOND:
                order = 3
            elif i.order == mol.Bond.PARTIAL_BOND:
                order = 4
            atom1 = var.index(i.atom1) + 1
            atom2 = var.index(i.atom2) + 1
            fileobj.write('%3d%3d%3d  0  0  0  0\n'
                          % (atom1, atom2, order))
        fileobj.write('M  END\n')
            
MOL_WRITERS = WriterManager()
MOL_WRITERS.writers['mdl'] = MOL_WRITERS.writers['mol'] = MdlMolWriter()

def write_mol(file, mol, type):
    return MOL_WRITERS.write(file, mol, type)

def format_mol(mol, type):
    return MOL_WRITERS.format(mol, type)
    
def write_molfile(filename, mol, type=None):
    return MOL_WRITERS.write_file(filename, mol, type)
