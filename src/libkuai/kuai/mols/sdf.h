#include <kuai/mols/io.h>


#ifndef _KUAI_MOLS_INCHI_H_
#define _KUAI_MOLS_INCHI_H_

namespace kuai {

	class SdfReader 
		: public BasicMoleculeReader
	{
	public:
		virtual MoleculePtr read(std::istream& stream) const;
	};

	class SdfWriter
		: public BasicMoleculeWriter 
	{
	public:
		virtual bool write(std::ostream& stream, const Molecule& mol) const;
	};

}

#endif
