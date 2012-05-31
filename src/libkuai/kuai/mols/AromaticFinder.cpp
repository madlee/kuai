#include <boost/tuple/tuple.hpp>
#include <kuai/mols/AromaticFinder.h>
#include <kuai/mols/RingFinder.h>

namespace kuai {

	namespace {
		void set_aromatic_candidates_flag(kuai::MoleculePtr& mol, std::vector<bool>& flags) {
			for (Index i = 0; i < flags.size(); ++i) {
				if (flags[i]) {
					AtomPtr atomI = mol->get_atom(i);
					Index deg = mol->degree(atomI);
					if (deg < 2) {
						flags[i] = false;
					}
					else {
						switch (atomI->number()) {
						case 6:
							if (atomI->charge == 0) {
								Index double_bond = 0;
								Index partial_bond = 0;
								for (Index j = 0; j < deg; ++j) {
									BondPtr bondJ = mol->get_neighbor_bond(atomI, j);
									if (bondJ->order == TRIPLE_BOND) {
										double_bond = 0;
										break;
									}
									else if (bondJ->order == DOUBLE_BOND) {
										double_bond += 1;
									}
									else if (bondJ->order == PARTIAL_BOND) {
										partial_bond += 1;
									}
								}
								flags[i] = (double_bond == 1 || partial_bond >= 2);
							}
							else {
								flags[i] = true;
							}
							break;

						case 7:	// N
						case 15: // P
						case 33: // As
							flags[i] = true;
							break;

						case 8:		// O
						case 16:	// S
						case 34:	// S
							flags[i] = true;
							break;

						default:
							flags[i] = false;
							break;
						}
					}
				}
			}
		}
	}

	AromaticFinder::AromaticFinder(MoleculePtr mol) 
		: _mol(mol)
	{
		Array<bool> flags(mol->count_atoms(), true);
		set_aromatic_candidates_flag(mol, flags);
		_setup(mol, flags);
	}

	AromaticFinder::AromaticFinder(MoleculePtr mol, Array<bool>& flags)
		: _mol(mol)
	{
		set_aromatic_candidates_flag(mol, flags);
		_setup(mol, flags);
	}

	void AromaticFinder::_setup(MoleculePtr mol, Array<bool>& hints) {
		RingFinder rings(mol, hints);
		RingFinder::RingBlockIterator it, itend;
		for (boost::tie(it, itend) = rings.iter_blocks(); it != itend; ++it) {
						
		}
	}
}

