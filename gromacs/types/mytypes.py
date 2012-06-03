from kuai.basicTools.xyz import XYZ

class Atom:
    def __init__(self, id=None, symbol=None, coords=None):
        self.id = id
        self.symbol = symbol
        self.coords = coords
        
    @staticmethod
    def parse(line):
        name = line[4:11].strip()
        element = line[13:17].strip()
        coords = line[30:54].strip().split()
        assert len(coords) == 3
        x = float(coords[0])
        y = float(coords[1])
        z = float(coords[2])
        
        return Atom(name, element, XYZ(x, y, z))
        
    def __str__(self):
        return self.id
    
class Bond:
    def __init__(self, a1, a2, order):
        self.atom1 = a1
        self.atom2 = a2
        self.order = order
