#include <kuai/typedef.h>

#ifndef _KUAI_MOLS_CONNECT_TABLE_H_
#define _KUAI_MOLS_CONNECT_TABLE_H_

namespace kuai {

	///> A connect table
	/***************************************************************************/
	/* This class represent a Connect Table. It maintains the connection 
	/* infomation in an array of index. So we can share the atoms and bonds 
	/* data independently to the connection data. It is great if we want to
	/* retrieve and manipulate some fragments of a molecule. The only thing 
	/* you need to do is re-construct a connect table and then use it. There
	/* is no need to copy the atoms and bonds array again. It saved many time
	/* and mem, I think so.
	/* This class is used by Molecule and should not be used directly in most 
	/* situation.
	/***************************************************************************/
	class ConnectTable : Noncopyable
	{
	public:
		explicit ConnectTable() 
		{ }

		void setup(const Index v[], Index max_atoms = INVALID_INDEX, Index max_bonds = INVALID_INDEX);

	public:
		Index count_atoms() const {
			return _atoms[-2];
		}
		Index count_bonds() const {
			return _atoms[-1];
		}
		Index get_atom(Index index) const {
			assert (index < count_atoms());
			return _atoms[index];
		}
		Index get_bond(Index index) const {
			assert (index < count_bonds());
			return _bonds[3 * index];
		}
		Index get_atom_index(const Index bond, Index index) const {
			assert (index < 2);
			Index bID = index_of_bond(bond);
			return _bonds[3*bID+1+index];
		}
		Index get_atom(const Index bond, Index index) const {
			return get_atom(get_atom_index(bond, index));
		}

		Index degree(Index a) const {
			Index id = index_of_atom(a);
			return _degrees[id];
		}
		Index get_neighbor_atom_index(Index a, Index j) const {
			Index id = index_of_atom(a);
			id = _offsets[id];
			return _adjs[id+j*2];
		}
		Index get_neighbor_atom(Index a, Index j) const {
			return get_atom(get_neighbor_atom_index(a, j));
		}
		Index get_neighbor_bond_index(Index a, Index j) const {
			Index id = index_of_atom(a);
			id = _offsets[id];
			return _adjs[id+j*2+1];
		}
		Index get_neighbor_bond(Index a, Index j) const {
			return get_bond(get_neighbor_bond_index(a, j));
		}
		Index get_bond(Index a1, Index a2) const {
			Index d = degree(a1);
			for (Index i = 0; i < d; ++i) {
				if (get_neighbor_atom(a1, i) == a2) {
					return get_neighbor_bond(a1, i);
				}
			}
			return INVALID_INDEX;
		};
		Index index_of_atom(Index a) const {
			return _ratoms[a];
		};
		Index index_of_bond(Index b) const {
			return _rbonds[b];
		};

		bool empty() const {
			return count_atoms() == 0;
		};

	private:
		std::vector<Index> _vArray;	// Array of index.
		Index* _atoms;
		Index* _bonds;
		Index* _degrees;
		Index* _offsets;
		Index* _adjs;
		Index* _ratoms;
		Index* _rbonds;
	};

}

#endif
