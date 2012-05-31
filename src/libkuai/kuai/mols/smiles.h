#include <kuai/mols/io.h>


#ifndef _KUAI_MOLS_SMILES_H_
#define _KUAI_MOLS_SMILES_H_

namespace kuai {

	/** parse molecule according to the input string */
	MoleculePtr parse_smiles(const String& v);
	String smiles(MoleculePtr mol, bool aromatic=false);
	String unique_smiles(MoleculePtr mol, bool aromatic=false);

	class SmilesReader
		: public BasicMoleculeReader
	{
	public:
		virtual bool read(std::istream& stream, Mutex& mutex, const FileName* filename, size_t count, Reference data);
	};

}

#endif 
