#include <kuai/tools/common.h>
#include <kuai/mols/hybird.h>


namespace kuai {

	static const Index TAUTOMERIC_ATOMS[] = {7, 8, 15, 16, 33, 34};

	HYBIRD_STATE get_hybird_state(MoleculePtr& mol, AtomPtr atom) {
		const Element* e = atom->element();
		if (e == NULL) {
			return HYBIRD_UNKNOWN;
		}
		
		Index deg = mol->degree(atom);
		switch (deg) {
		case 7:
			return HYBIRD_SP3D3;

		case 6:
			return HYBIRD_SP3D2;

		case 5:
			return HYBIRD_SP3D;

		case 4:
			return HYBIRD_SP3;

		default: 
			break;
		}
		
		Index total_order = 0;
		Index unsat = 0;
		for (Index i = 0; i < deg; ++i) {
			BondPtr bond_i = mol->get_neighbor_bond(atom, i);
			total_order += bond_i->order;
			unsat += bond_i->order - SINGLE_BOND;
		}

		switch (unsat) {
		case 0:
			// Check implicit hydrogen
			break;

		case SINGLE_BOND*1:
			return HYBIRD_SP2;

		case SINGLE_BOND*2:
			return HYBIRD_SP;

		default:
			return HYBIRD_UNKNOWN;
		}

		Index implicit_h = 0;
		if (!e->do_not_add_H) {
			Index offset = atom->charge - Element::MIN_ATOM_CHARGE;
			if (offset < 5) {
				implicit_h = (e->valence[offset][0] * SINGLE_BOND - total_order) / total_order;
			}
		}
		deg += implicit_h;

		switch (deg) {
		case 0:
			return HYBIRD_ION;

		case 1:
			return HYBIRD_S;

		case 2:
		case 3:
			// Check taut atoms
			break;
		
		case 4:
			return HYBIRD_SP3;

		default:
			return HYBIRD_UNKNOWN;
		}

		static const Index* const PEND = TAUTOMERIC_ATOMS+ARRAY_LENGTH(TAUTOMERIC_ATOMS);
		const Index* const p = std::find(TAUTOMERIC_ATOMS, PEND, atom->number());

		if (p != PEND) {
			deg -= implicit_h;
			for (Index i = 0; i < deg; ++i) {
				AtomPtr atom_i = mol->get_neighbor_atom(atom, i);
				Index deg_j = mol->degree(atom_i);
				for (Index j = 0; j < deg_j; ++j) {
					BondPtr bond_j = mol->get_neighbor_bond(atom_i, j);
					if (bond_j->order == DOUBLE_BOND) {
						return HYBIRD_SP2;
					}
				}
			}		
		}
		
		return HYBIRD_SP3;
	}

	void get_hybird_states(MoleculePtr& mol, Array<HYBIRD_STATE>& result) {
		Index natoms = mol->count_atoms();
		result.clear();
		result.reserve(natoms);
		for (Index i = 0; i < natoms; ++i) {
			AtomPtr atom_i = mol->get_atom(i);
			result.push_back(get_hybird_state(mol, atom_i));
		}
	}

}