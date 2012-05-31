#include <kuai/mols/Molecule.h>
#include <kuai/mols/RingFinder.h>


#ifndef _KUAI_MOLS_AROMATIC_FINDER_H_
#define _KUAI_MOLS_AROMATIC_FINDER_H_

namespace kuai {

	class AromaticFinder {
	public:
		explicit AromaticFinder(MoleculePtr mol);
		explicit AromaticFinder(MoleculePtr mol, Array<bool>& flags);
		explicit AromaticFinder(MoleculePtr mol, AromaticFinder& superset);

		bool is_aromatic(AtomPtr atom) const {
			return _aromatic.get() != NULL && _aromatic->index(atom) != INVALID_INDEX;
		}
		bool is_aromatic(BondPtr bond) const {
			return _aromatic.get() != NULL && _aromatic->index(bond) != INVALID_INDEX;
		}

	private:
		MoleculePtr _mol;
		MoleculePtr _aromatic;
		SharedPtr<RingFinder> _rings;

	private:
		void _setup(MoleculePtr mol, Array<bool>& hints);
	};

}
#endif
