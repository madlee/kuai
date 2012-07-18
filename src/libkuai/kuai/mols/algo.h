// BasicMolecularVisitor and broad 1st travel and depth 1st travel are defined in this file.

#include <queue>

#include <kuai/typedef.h>
#include <kuai/tools/common.h>
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

	class TraceBackVisitor
		: public BasicMolecularVisitor
	{
	public:
		void find(AtomPtr target, AtomPtr source) 
		{ 
			_sources[target] = source;
		}

		AtomPtr get_source(AtomPtr p) {
			KuaiMap<AtomPtr, AtomPtr>::const_iterator it = _sources.find(p);
			if (it == _sources.end()) {
				return NULL;
			}
			else {
				return it->second;
			}
		}

	protected:
		KuaiMap<AtomPtr, AtomPtr> _sources;
	};

	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr& mol, VisitorType& visitor, VisitedFlag flags[], typename VisitorType::AtomPtr v1, typename VisitorType::AtomPtr from = NULL)  {
		visitor.start(v1);
		flags[mol->index(v1)] = GRAY_FLAG;
		Index nDegree = mol->degree(v1);
		for (Index i = 0; i < nDegree; ++i) { 
			typename VisitorType::AtomPtr neighbor = mol->get_neighbor_atom(v1, i);
			if (neighbor != from) {
				Index iNeighbor = mol->index(neighbor);
				if (flags[iNeighbor] == WHITE_FLAG) { 
					visitor.find(neighbor, v1);
					flags[iNeighbor] = GRAY_FLAG;
					dft(mol, visitor, flags, neighbor, v1);
				}
				else { 
					visitor.back(neighbor, v1);
				}
			}
		}
		flags[mol->index(v1)] = BLACK_FLAG;
		visitor.finish(v1);
	}

	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr& mol, VisitorType& visitor, VisitedFlag flags[]) {
		Index nAtoms = mol->count_atoms();
		for (Index i = 0; i < nAtoms; ++i) { 
			if (flags[i] == WHITE_FLAG) { 
				dft(mol, visitor, flags, mol->get_atom(i));
			}
		}
	}
	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr& mol, VisitorType& visitor)  { 
		Array<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
		dft(mol, visitor, &flags[0]);
	};
	template<typename VisitorType>
	void dft(typename VisitorType::MolPtr& mol, VisitorType& visitor, typename VisitorType::AtomPtr v1)  { 
		Array<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
		dft(mol, visitor, &flags[0], v1);
	}


	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr& mol, VisitorType& visitor, VisitedFlag flags[], std::queue<typename VisitorType::AtomPtr>& queue) {
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
	void bft(typename VisitorType::MolPtr& mol, VisitorType& visitor, VisitedFlag flags[], typename VisitorType::AtomPtr v1)  { 
		std::queue<typename VisitorType::AtomPtr> queue; 
		queue.push(v1);
		flags[mol->index(v1)] = GRAY_FLAG;
		bft(mol, visitor, flags, queue);
	}
	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr& mol, VisitorType& visitor, typename VisitorType::AtomPtr v1) { 
		std::vector<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
		bft(mol, visitor, &flags[0], v1);
	}
	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr& mol, VisitorType& visitor, VisitedFlag flags[]) { 
        Index nAtoms = mol->count_atoms();
		for (Index i = 0; i < nAtoms; ++i) { 
			if (flags[i] == WHITE_FLAG) { 
				bft(mol, visitor, flags, mol->get_atom(i));
			}
		}
	}
	template<typename VisitorType>
	void bft(typename VisitorType::MolPtr& mol, VisitorType& visitor)  { 
		std::vector<VisitedFlag> flags(mol->countAtoms(), WHITE_FLAG);
		bft(mol, visitor, &flags[0]);
	}

	void split(MoleculePtr& mol, std::vector<MoleculePtr>& result);	

	inline std::vector<MoleculePtr> split(MoleculePtr& mol) {
		std::vector<MoleculePtr> result;
		split(mol, result);
		return result;
	}


	template<typename IteratorType>
	MoleculePtr compine_fragments(MoleculePtr& parent, IteratorType it1, IteratorType it2) {
		Array<AtomPtr> atoms;
		Array<BondPtr> bonds;
		Index total_atoms = 0, total_bonds = 0;
		for (;it1 != it2; ++it1) {
			total_atoms += (*it1)->count_atoms();
			total_bonds += (*it1)->count_bonds();
		}
		atoms.reserve(total_atoms);
		bonds.reserve(total_bonds);
		for (;it1 != it2; ++it1) {
			Index natoms = (*it1)->count_atoms();
			for (Index j = 0; j < natoms; ++j) {
				atoms.push_back((*it1)->get_atom(j));
			}
			Index nbonds = (*it1)->count_bonds();
			for (Index j = 0; j < nbonds; ++j) {
				bonds.push_back((*it1)->get_bond(j));
			}
		}
		unique(atoms);
		unique(bonds);

		return MoleculePtr(new Molecule(parent, atoms, bonds));
	}

}

#endif
