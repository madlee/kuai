#include <inchi_api.h>
#include <kuai/typedef.h>
#include <kuai/mols/Molecule.h>
#include <kuai/mols/io.h>


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
	
	void inchi(MoleculePtr mol, const String& options, InchiOutput& output);
	inline String inchi(MoleculePtr mol, const String& options="") {
		InchiOutput output;
		inchi(mol, options, output);
		return String(output.szInChI);
	}
	String inchi_key(const String& inchi_code);

	inline std::pair<String, String> inchi_and_key(MoleculePtr mol, const String& options = "") {
		String inchi_code = inchi(mol, options);
		return std::make_pair(inchi_code, inchi_key(inchi_code));
	}

	class InchiReader 
		: public BasicMoleculeReader
	{
	public:
		virtual bool read(std::istream& stream, Mutex& mutex, const FileName* filename, size_t count, Reference data);

		void set_option(const String& v) {
			_option = v;
		}

		const String& get_option() const {
			return _option;
		}

	private:
		String _next_name;
		String _next_inchi;
		String _option;

	private:
		bool _to_next_record(std::istream& stream, std::ostream& os);
		
	};

}

#endif 
