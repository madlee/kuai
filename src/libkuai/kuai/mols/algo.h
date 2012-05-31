// BasicMolecularVisitor and broad 1st travel and depth 1st travel are defined in this file.

#include <queue>

#include <kuai/typedef.h>
#include <kuai/mols/Molecule.h>


#ifndef _KUAI_MOLS_ALGO_H_
#define _KUAI_MOLS_ALGO_H_

namespace kuai {

	enum VisitedFlag
	{
		WHITE_FLAG,
		GRAY_FLAG,
		BLACK_FLAG
	};


	class BasicMolecularVisitor {
	public:
		typedef Molecule			MolType;
		typedef SharedPtr<Molecule>	MolPtr;
		typedef Atom				AtomType;
		typedef Atom*				AtomPtr;
		typedef Bond				BondType;
		typedef Bond*				BondPtr;

	public:
		void start(AtomPtr) 
		{ }
		void finish(AtomPtr) 
		{ }
		void find(AtomPtr target, AtomPtr source) 
		{ }
		void back(AtomPtr target, AtomPtr source) 
		{ }
	};

	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr mol, VisitorType& visitor, VisitedFlag flags[], typename VisitorType::AtomPtr v1)  {
		visitor.start(v1);
		flags[mol->index(v1)] = GRAY_FLAG;
		Index nDegree = mol->degree(v1);
		for (Index i = 0; i < nDegree; ++i) { 
			typename VisitorType::AtomPtr neighbor = mol->get_neighbor_atom(v1, i);
			Index iNeighbor = mol->index(neighbor);
			if (flags[iNeighbor] == WHITE_FLAG) { 
				visitor.find(neighbor, v1);
				dft(mol, visitor, flags, neighbor);
			}
			else { 
				visitor.back(neighbor, v1);
			}
		}
		flags[mol->index(v1)] = BLACK_FLAG;
		visitor.finish(v1);
	}
	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr mol, VisitorType& visitor, VisitedFlag flags[]) {
		Index nAtoms = mol->count_atoms();
		for (Index i = 0; i < nAtoms; ++i) { 
			if (flags[i] == WHITE_FLAG) { 
				dft(mol, visitor, flags, mol->get_atom(i));
			}
		}
	}
	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr mol, VisitorType& visitor)  { 
		Array<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
		dft(mol, visitor, &flags[0]);
	};
	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr mol, VisitorType& visitor, typename VisitorType::AtomPtr v1)  { 
		Array<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
		dft(mol, visitor, &flags[0], v1);
	}


	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr mol, VisitorType& visitor, VisitedFlag flags[], std::queue<typename VisitorType::AtomPtr>& queue) {
		for (; !queue.empty(); queue.pop()) {
			typename VisitorType::AtomPtr a = queue.front();
			Index i = mol->index(a);
			assert (flags[i] == GRAY_FLAG);
			visitor.start(a);
			Index nDegree = mol->degree(a);
			for (Index j = 0; j < nDegree; ++j) { 
				typename VisitorType::AtomPtr neighbor = mol->get_neighbor_atom(a, j);
				Index iNeighbor = mol->index(neighbor);
				if (flags[iNeighbor] == WHITE_FLAG) { 
					visitor.find(neighbor, a);
					queue.push(neighbor);
					flags[iNeighbor] = GRAY_FLAG;
				}
				else { 
					visitor.back(neighbor, a);
				}
			}
			visitor.finish(a);
			flags[i] = BLACK_FLAG;
		}
	}

	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr mol, VisitorType& visitor, VisitedFlag flags[], typename VisitorType::AtomPtr v1)  { 
		std::queue<typename VisitorType::AtomPtr> queue; 
		queue.push(v1);
		flags[mol->index(v1)] = GRAY_FLAG;
		bft(mol, visitor, flags, queue);
	}
	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr mol, VisitorType& visitor, typename VisitorType::AtomPtr v1) { 
		std::vector<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
		bft(mol, visitor, &flags[0], v1);
	}
	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr mol, VisitorType& visitor, VisitedFlag flags[]) { 
        Index nAtoms = mol->count_atoms();
		for (Index i = 0; i < nAtoms; ++i) { 
			if (flags[i] == WHITE_FLAG) { 
				bft(mol, visitor, flags, mol->get_atom(i));
			}
		}
	}
	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr mol, VisitorType& visitor)  { 
		std::vector<VisitedFlag> flags(mol->countAtoms(), WHITE_FLAG);
		bft(mol, visitor, &flags[0]);
	}

	void split(MoleculePtr mol, std::vector<MoleculePtr>& result);	

	inline std::vector<MoleculePtr> split(MoleculePtr mol) {
		std::vector<MoleculePtr> result;
		split(mol, result);
		return result;
	}

}

#endif
