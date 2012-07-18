#include <inchi_api.h>
#include <kuai/typedef.h>
#include <kuai/mols/Molecule.h>


#ifndef _KUAI_MOLS_INCHI_H_
#define _KUAI_MOLS_INCHI_H_

namespace kuai {

	class InchiOutput
		: public inchi_Output, Noncopyable
	{
	public:
		InchiOutput() {
			memset(this, 0, sizeof(inchi_Output));
		}
		~InchiOutput() {
			FreeStdINCHI(this);
		}
	};

	class InchiOutputStruct
		: public inchi_OutputStruct, Noncopyable
	{
	public:
		InchiOutputStruct() {
			memset(this, 0, sizeof(inchi_OutputStruct));
		
		}
		~InchiOutputStruct() {
			FreeStructFromStdINCHI(this);
		}

		MoleculePtr get_result() const;
	};

	MoleculePtr parse_inchi(String v, String options, InchiOutputStruct& output);
	/** parse molecule according to the input string */
	inline MoleculePtr parse_inchi(const String& v, const String& options="") {
		InchiOutputStruct output;
		return parse_inchi(v, options, output);
	}
	
	void inchi(MoleculePtr& mol, const String& options, InchiOutput& output);
	inline String inchi(MoleculePtr& mol, const String& options="") {
		InchiOutput output;
		inchi(mol, options, output);
		return String(output.szInChI);
	}
	String inchi_key(const String& inchi_code);
	String inchi_key(MoleculePtr& mol, const String& options="");

	inline std::pair<String, String> inchi_and_key(MoleculePtr& mol, const String& options = "") {
		String inchi_code = inchi(mol, options);
		return std::make_pair(inchi_code, inchi_key(inchi_code));
	}

}

#endif 
