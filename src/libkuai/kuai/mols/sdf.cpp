#include <cctype>
#include <boost/tuple/tuple.hpp>
#include <boost/regex.hpp>
#include <kuai/tools/common.h>
#include <kuai/mols/sdf.h>
#include <kuai/mols/algo.h>
#include <kuai/mols/AromaticFinder.h>


namespace kuai {

	MoleculePtr SdfReader::read(std::istream& stream) const {
		return MoleculePtr();
	};

	bool SdfWriter::write(std::ostream& stream, const Molecule& mol) const {
		return true;
	};
}
