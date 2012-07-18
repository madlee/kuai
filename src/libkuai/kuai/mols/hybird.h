#include <kuai/typedef.h>
#include <kuai/mols/Molecule.h>

#ifndef _KUAI_MOLS_HYBIRD_H_
#define _KUAI_MOLS_HYBIRD_H_

namespace kuai {

	enum HYBIRD_STATE {
		HYBIRD_UNKNOWN = -1,
		HYBIRD_NONE,
		HYBIRD_ION,
		HYBIRD_S,
		HYBIRD_SP,
		HYBIRD_SP2,
		HYBIRD_SP3,
		HYBIRD_SP3D, 
		HYBIRD_SP3D2, 
		HYBIRD_SP3D3
	};

	HYBIRD_STATE get_hybird_state(MoleculePtr& mol, AtomPtr atom);
	void get_hybird_states(MoleculePtr& mol, Array<HYBIRD_STATE>& result);

}

#endif
