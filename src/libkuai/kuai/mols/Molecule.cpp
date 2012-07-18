#include <kuai/tools/error.h>
#include <kuai/mols/Molecule.h>


namespace kuai {

	Molecule::Molecule(AtomArray& atoms, BondArray& bonds, bool swap_array) 
		: _atoms(new AtomArray), _bonds(new BondArray)
	{
		if (swap_array) {
			_atoms->swap(atoms);
			_bonds->swap(bonds);
		}
		else {
			*_atoms = atoms;
			*_bonds = bonds;
		}
		_setup_0();
	}

	Molecule::Molecule(AtomArrayPtr atoms, BondArrayPtr bonds) 
		: _atoms(atoms), _bonds(bonds)
	{
		_setup_0();
	}

	
	Molecule::Molecule(SharedPtr<Molecule>& pmol, const Array<Atom*>& atoms, const Array<Bond*>& bonds) 
		: _atoms(pmol->_atoms), _bonds(pmol->_bonds), _chiral_buffer(pmol->_chiral_buffer)
	{
		_setup_ct(*pmol, atoms, bonds);
		_setup_chirals(pmol);
	}
	Molecule::Molecule(SharedPtr<Molecule>& pmol, const Array<Atom*>& atoms) 
		: _atoms(pmol->_atoms), _bonds(pmol->_bonds), _chiral_buffer(pmol->_chiral_buffer)
	{
		_setup(pmol, atoms);
		_setup_chirals(pmol);
	}

	Molecule::Molecule(SharedPtr<Molecule>& pmol, const Array<Bond*>& bonds) 
		: _atoms(pmol->_atoms), _bonds(pmol->_bonds), _chiral_buffer(pmol->_chiral_buffer)
	{
		_setup(pmol, bonds);
		_setup_chirals(pmol);
	}

	Molecule::~Molecule() {
		// Do nothing.
	}

	SharedPtr<Molecule> Molecule::clone() {
		AtomArray atoms;
		BondArray bonds;
		ChiralCenterArray chirals;

		Index nAtoms = count_atoms(), nBonds = count_bonds(), nChirals = count_chirals();
		atoms.reserve(nAtoms);
		bonds.reserve(nBonds);
		chirals.reserve(nChirals);

		for (Index i = 0; i < nAtoms; ++i) {
			atoms.push_back(*get_atom(i));
		} 
		for (Index i = 0; i < nBonds; ++i) {
			Bond* bondI = get_bond(i);
			bonds.push_back(*bondI);
			bonds.back().atom1 = index(get_atom(bondI, 0));
			bonds.back().atom2 = index(get_atom(bondI, 1));
		}

		for (Index i = 0; i < nChirals; ++i) {
			ChiralCenter* center = get_chiral_center(i);
			chirals.push_back(*center);

			if (Atom* atom0 = get_atom(center)) {
				chirals.back().center_atom = index(atom0);
			}
			else {
				chirals.back().center_atom = INVALID_INDEX;
			}

			for (Index j = 0; j < 4; ++j) {
				if (Atom* atomJ = get_atom(center, j)) {
					chirals.back().neighbor[j] = index(atomJ);
				}
				else {
					chirals.back().neighbor[j] = INVALID_INDEX;
				}
			}
		}

		SharedPtr<Molecule> result(new Molecule(atoms, bonds, chirals, true));
		result->name = name;
		result->text_info = text_info;
		return result;
	}

	void Molecule::_setup_ct(const Molecule& mol, const Array<Atom*>& atoms, const Array<Bond*>& bonds) {
		Index buffersize = atoms.size() + 3*bonds.size() + 2;
		Array<Index> buffer; buffer.reserve(buffersize);

		buffer.push_back(atoms.size());
		buffer.push_back(bonds.size());

		Array<Index> atom_index(mol._atoms->size(), INVALID_INDEX);
		for (Array<Atom*>::const_iterator
			it = atoms.begin(); it != atoms.end(); ++it)
		{
			ptrdiff_t delta = (*it)-&(mol._atoms->front());
			assert (0 <= delta && delta < mol._atoms->size());
			assert (atom_index[delta] == INVALID_INDEX);
			atom_index[delta] = it-atoms.begin();
			buffer.push_back(Index(delta));
		}

		for (Array<Bond*>::const_iterator
			it = bonds.begin(); it != bonds.end(); ++it)
		{
			ptrdiff_t delta = (*it)-&(mol._bonds->front());
			assert (0 <= delta && delta < mol._bonds->size());
			buffer.push_back(Index(delta));

			assert (atom_index[(**it).atom1] != INVALID_INDEX);
			buffer.push_back(atom_index[(**it).atom1]);

			assert (atom_index[(**it).atom2] != INVALID_INDEX);
			buffer.push_back(atom_index[(**it).atom2]);
		}

		_ct.setup(&buffer.front(), mol._atoms->size(), mol._bonds->size());
	}

	void Molecule::_setup_0() 
	{
		Index buffersize = _atoms->size() + 3*_bonds->size() + 2;
		Array<Index> buffer;
		buffer.reserve(buffersize);
		buffer.push_back(_atoms->size());
		buffer.push_back(_bonds->size());
		for (Index i = 0; i < _atoms->size(); ++i) {
			buffer.push_back(i);
		}
		for (Index i = 0; i < _bonds->size(); ++i) {
			buffer.push_back(i);
			const Bond& bondI = (*_bonds)[i];
			
			assert (bondI.atom1 < _atoms->size());
			assert (bondI.atom2 < _atoms->size());

			buffer.push_back(bondI.atom1);
			buffer.push_back(bondI.atom2);
		}

		_ct.setup(&buffer.front(), _atoms->size(), _bonds->size());
	}

	void Molecule::_setup(SharedPtr<Molecule>& pmol, const Array<Atom*>& atoms)
	{
		BitSet bitmasks(pmol->_atoms->size());
		for (Array<Atom*>::const_iterator
			it = atoms.begin(); it != atoms.end(); ++it)
		{
			ptrdiff_t delta = *it-&pmol->_atoms->front();
			assert (0 <= delta && delta < bitmasks.size());
			bitmasks.set(delta);
		}
		
		Index nBonds = pmol->count_bonds();
		Array<Bond*> bonds; bonds.reserve(nBonds);
		for (Index i = 0; i < nBonds; ++i) {
			Bond* bondI = pmol->get_bond(i);
			if (bitmasks.test(bondI->atom1) && bitmasks.test(bondI->atom2)) {
				bonds.push_back(bondI);
			}
		}
		_setup_ct(*pmol, atoms, bonds);
	}

	void Molecule::_setup(SharedPtr<Molecule>& pmol, const Array<Bond*>& bonds) 
	{ 
		Array<Atom*> atoms; atoms.reserve(bonds.size()*2);

		for (Array<Bond*>::const_iterator
			it = bonds.begin(); it != bonds.end(); ++it)
		{
			Bond* bondI = (*it);
			atoms.push_back(pmol->get_atom(bondI, 0));
			atoms.push_back(pmol->get_atom(bondI, 1));
		}

		std::sort(atoms.begin(), atoms.end());
		Array<Atom*>::iterator it = std::unique(atoms.begin(), atoms.end());
		atoms.erase(it, atoms.end());

		_setup_ct(*pmol, atoms, bonds);
	}

	namespace {
		Index count_implicit_H(Molecule* mol) {
			Index result = 0;
			Index nAtoms = mol->count_atoms();
			Array<Index> orders(nAtoms, 0);

			Index nBonds = mol->count_bonds();
			for (Index i = 0; i < nBonds; ++i) {
				BondPtr bondI = mol->get_bond(i);
				AtomPtr atom1 = mol->get_atom(bondI, 0);
				AtomPtr atom2 = mol->get_atom(bondI, 1);
				orders[mol->index(atom1)] += bondI->order;
				orders[mol->index(atom2)] += bondI->order;
			}

			for (Index i = 0; i < nAtoms; ++i) {
				AtomPtr atomI = mol->get_atom(i);
				if (atomI->num_iso_H[0] == INVALID_INDEX) {
					const Element* e = atomI->element();
					if (!e->do_not_add_H && Element::MIN_ATOM_CHARGE <= atomI->charge && atomI->charge <= Element::MAX_ATOM_CHARGE) {
						Integer v = Integer(e->valence[atomI->charge - Element::MIN_ATOM_CHARGE][0]);
						v -= orders[i] / SINGLE_BOND;
						if (v > 0) {
							result += v;
						}
					}
				}
				else {
					result += atomI->num_iso_H[0];
				}
			}

			return result;
		}
	}

	String Molecule::formula() {
		Index nC = 0, nH = 0;
		ChargeType nCharge = 0;
		
		Index nAtoms = count_atoms();
		OrderedMap<String, Index> elements;
		for (Index i = 0; i < nAtoms; ++i) {
			AtomPtr atomI = get_atom(i);
			if (atomI->symbol() == "H") {
				nH += 1;
			}
			else if (atomI->symbol() == "C") {
				nC += 1;
			}
			else {
				elements[atomI->symbol()] += 1;
			}
			nCharge += atomI->charge; 
		}
		nH += count_implicit_H(this);
		std::ostringstream result;
		if (nC > 0) {
			result << "C";
			if (nC > 1) {
				result << nC;
			}
		}
		if (nH > 0) {
			result << "H";
			if (nH > 1) {
				result << nH;
			}
		}

		for (OrderedMap<String, Index>::const_iterator
			it = elements.begin(); it != elements.end(); ++it)
		{
			result << it->first;
			if (it->second > 1) {
				result << it->second;
			}
		}

		if (nCharge > 0) {
			result << "+";
			if (nCharge > 1) {
				result << nCharge;
			}
		}
		else if (nCharge < 0) {
			result << "-";
			if (nCharge < 1) {
				result << -nCharge;
			}
		}

		return result.str();
	}
	RealNumber Molecule::weight() {
		static const Element* H = Element::get(1);

		Index nAtoms = count_atoms();
		RealNumber result = 0;
		for (Index i = 0; i < nAtoms; ++i) {
			AtomPtr atomI = get_atom(i);
			result += atomI->weight();
		}
		Index nH = count_implicit_H(this);
		result += nH * H->weight;
		return result;
	}


	Molecule::Molecule(AtomArray& atoms, BondArray& bonds, ChiralCenterArray& chirals, bool swap_array) 
		: _atoms(new AtomArray), _bonds(new BondArray), _chiral_buffer(new ChiralCenterArray)
	{
		if (swap_array) {
			_atoms->swap(atoms);
			_bonds->swap(bonds);
			_chiral_buffer->swap(chirals);
		}
		else {
			*_atoms = atoms;
			*_bonds = bonds;
			*_chiral_buffer = chirals;
		}
		_setup_0();
		_chiral_centers.reserve(chirals.size());
		for (Index i = 0; i < chirals.size(); ++i) {
			const ChiralCenter& c = chirals[i];
			assert (c.center_atom == INVALID_INDEX || c.center_atom < atoms.size());
			assert (c.neighbor[0] == INVALID_INDEX || c.neighbor[0] < atoms.size());
			assert (c.neighbor[1] == INVALID_INDEX || c.neighbor[1] < atoms.size());
			assert (c.neighbor[2] == INVALID_INDEX || c.neighbor[2] < atoms.size());
			assert (c.neighbor[3] == INVALID_INDEX || c.neighbor[3] < atoms.size());
			_chiral_centers.push_back(i);
		}
	}
	Molecule::Molecule(SharedPtr<Molecule>& pmol, const std::vector<Atom*>& atoms, const std::vector<Bond*>& bonds, const std::vector<ChiralCenter*>& chirals)
		: _atoms(pmol->_atoms), _bonds(pmol->_bonds), _chiral_buffer(pmol->_chiral_buffer)
	{
		_setup_ct(*pmol, atoms, bonds);
		_setup_chirals(chirals);
	}
		
	Molecule::Molecule(SharedPtr<Molecule>& pmol, const std::vector<Atom*>& atoms, const std::vector<ChiralCenter*>& chirals) 
		: _atoms(pmol->_atoms), _bonds(pmol->_bonds), _chiral_buffer(pmol->_chiral_buffer)
	{
		_setup(pmol, atoms);
		_setup_chirals(chirals);
	}
		
	Molecule::Molecule(SharedPtr<Molecule>& pmol, const std::vector<Bond*>& bonds, const std::vector<ChiralCenter*>& chirals) 
		: _atoms(pmol->_atoms), _bonds(pmol->_bonds), _chiral_buffer(pmol->_chiral_buffer)
	{
		_setup(pmol, bonds);
		_setup_chirals(chirals);
	}

	ChiralCenterPtr Molecule::get_chiral_center(AtomPtr atom) {
		assert (atom != NULL);

		Index n = count_chirals();
		for (Index i = 0;  i < n; ++i) {
			ChiralCenterPtr center = get_chiral_center(i);
			if (get_atom(center) == atom) {
				return center;
			}
		}
		return NULL;
	};
	ChiralCenterPtr Molecule::get_chiral_center(BondPtr bond) {
		assert (bond != NULL);
		assert (index(bond) != INVALID_INDEX);

		Index n = count_chirals();
		for (Index i = 0;  i < n; ++i) {
			ChiralCenterPtr center = get_chiral_center(i);
			if (get_bond(center) == bond) {
				return center;
			}
		}
		return NULL;
	}

	AtomPtr Molecule::get_atom(ChiralCenterPtr center) {
		assert (center != NULL);
		Index v = center->center_atom;
		if (v != INVALID_INDEX && v < _atoms->size()) {
			AtomPtr result = &(*_atoms)[v];
			if (index(result) != INVALID_INDEX) {
				return result;
			}
		}
		return NULL;
	}
	AtomPtr Molecule::get_atom(ChiralCenterPtr center, Index i) {
		assert (center != NULL);
		assert (i < 4);
		Index v = center->neighbor[i];
		if (v != INVALID_INDEX && v < _atoms->size()) {
			AtomPtr result = &(*_atoms)[v];
			if (index(result) != INVALID_INDEX) {
				return result;
			}
		}
		return NULL;
	}
	BondPtr Molecule::get_bond(ChiralCenterPtr center) {
		assert (center != NULL);
		if (center->type == ChiralCenter::STEREO_TYPE_DOUBLEBOND) {
			assert (center->neighbor[1] < _atoms->size() && center->neighbor[2] < _atoms->size());
			AtomPtr a1 = &(*_atoms)[center->neighbor[1]];
			AtomPtr a2 = &(*_atoms)[center->neighbor[2]];
			return get_bond(a1, a2);
		}
		return NULL;
	}
	BondPtr Molecule::get_bond(ChiralCenterPtr center, Index i) {
		assert (center != NULL);
		switch (center->type) {
		case ChiralCenter::STEREO_TYPE_DOUBLEBOND: 
			assert (i < 3);
			assert (center->neighbor[i] < _atoms->size() && center->neighbor[i+1] < _atoms->size());
			{
				AtomPtr a1 = &(*_atoms)[center->neighbor[i]];
				AtomPtr a2 = &(*_atoms)[center->neighbor[i+1]];
				return get_bond(a1, a2);
			}

		case ChiralCenter::STEREO_TYPE_TETRAHEDRAL:
			assert (i < 4);
			assert (center->center_atom < _atoms->size() && center->neighbor[i] < _atoms->size());
			{
				AtomPtr a1 = &(*_atoms)[center->center_atom];
				AtomPtr a2 = &(*_atoms)[center->neighbor[i]];
				return get_bond(a1, a2);
			}

		default:
			throw error("Undefined Stereo type %d.", center->type);
			break;
		}
		
		return NULL;
	}

	Index Molecule::index(const ChiralCenterPtr center) const {
		if (_chiral_buffer.get()) {
			ptrdiff_t delta = center - &_chiral_buffer->front();
			Index result = std::find(_chiral_centers.begin(), _chiral_centers.end(), Index(delta)) - _chiral_centers.begin();
			if (result < _chiral_centers.size()) {
				return result;
			}
		}
		return INVALID_INDEX;
	}

	void Molecule::_setup_chirals(SharedPtr<Molecule>& pmol) {
		if (Index n = pmol->count_chirals()) {
			std::vector<ChiralCenter*> chirals;
			chirals.reserve(n);
			for (Index i = 0; i < n; ++i) {
				ChiralCenter* center = pmol->get_chiral_center(i);
				chirals.push_back(center);
			}
			_setup_chirals(chirals);
		}
	}

	void Molecule::_setup_chirals(const std::vector<ChiralCenter*>& chirals) {
		assert (_chiral_buffer.get() != NULL);
		for (std::vector<ChiralCenter*>::const_iterator 
			it = chirals.begin(); it != chirals.end(); ++it)
		{
			ChiralCenter* center = *it;

			switch (center->type) {
			case ChiralCenter::STEREO_TYPE_DOUBLEBOND: 
				{
					assert (center->neighbor[1] != INVALID_INDEX);
					assert (center->neighbor[2] != INVALID_INDEX);
					AtomPtr a1 = &(*_atoms)[center->neighbor[1]];
					AtomPtr a2 = &(*_atoms)[center->neighbor[2]];
					if (index(a1) == INVALID_INDEX || index(a2) == INVALID_INDEX || get_bond(a1, a2) == NULL) {
						continue;
					}
					
					if (center->neighbor[0] != INVALID_INDEX) {
						AtomPtr a0 = &(*_atoms)[center->neighbor[0]];
						if (index(a0) == INVALID_INDEX || get_bond(a0, a1) == NULL) {
							continue;
						}
					}

					if (center->neighbor[3] != INVALID_INDEX) {
						AtomPtr a3 = &(*_atoms)[center->neighbor[3]];
						if (index(a3) == INVALID_INDEX || get_bond(a2, a3) == NULL) {
							continue;
						}
					}
				}
				break;

			case ChiralCenter::STEREO_TYPE_TETRAHEDRAL:
				{
					assert (center->center_atom != INVALID_INDEX);
					AtomPtr pCenter = &(*_atoms)[center->center_atom];
					if (index(pCenter) == INVALID_INDEX || degree(pCenter) < 3) {
						continue;
					}
					Index j = 0;
					for (j = 0; j < 4; ++j) {
						if (center->neighbor[j] != INVALID_INDEX) {
							AtomPtr pA = &(*_atoms)[center->neighbor[j]];
							BondPtr pB = get_bond(pA, pCenter);
							if (index(pA) == INVALID_INDEX || pB == NULL) {
								j = 0;
								break;
							}
						}
					}
					if (j != 4) {
						continue;
					}
				}
				break;

			default:
				throw error("Undefined Stereo type.", center->type);
				break;
			}

			Index delta = Index(center - &_chiral_buffer->front());
			assert (delta < _chiral_buffer->size());
			_chiral_centers.push_back(delta);
		}
	}

	MoleculePtr drop_hydrogen(MoleculePtr& mol) {
		Index nAtoms = mol->count_atoms();
		Array<AtomPtr> atoms; atoms.reserve(nAtoms);
		for (Index i = 0; i < nAtoms; ++i)  {
			AtomPtr atomI = mol->get_atom(i);
			if (atomI->symbol() != "H") {
				atoms.push_back(atomI);
			}
		}
		if (atoms.size() == nAtoms) {
			return mol;
		}
		else {
			return MoleculePtr(new Molecule(mol, atoms));
		}
	}



}
