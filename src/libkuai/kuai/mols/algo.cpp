#include <kuai/mols/algo.h>

namespace kuai { 
	namespace {
		class ConnectPartFinder
			: public BasicMolecularVisitor
		{
		public:
			ConnectPartFinder(MoleculePtr mol) { 
				atoms.reserve(mol->count_atoms());
			}

			void start(AtomPtr v) {
				atoms.push_back(v);
			}

			MoleculePtr get_result(MoleculePtr mol) {
				return MoleculePtr(new Molecule(mol, atoms));
			}

			void clear() {
				atoms.clear();
			}

		private:
			std::vector<Atom*> atoms;
		};
	}


	void split(MoleculePtr mol, std::vector<MoleculePtr>& result) {
		Index nAtoms = mol->count_atoms();
		std::vector<VisitedFlag> flags(nAtoms, WHITE_FLAG);

		result.clear();

		ConnectPartFinder finder(mol);
		for (Index i = 0; i < nAtoms; ++i) {
			finder.clear();
			if (flags[i] == WHITE_FLAG) {
				dft(mol, finder, &flags[0], mol->get_atom(i));
				result.push_back(finder.get_result(mol));
			}					
		}
	}

	

}
