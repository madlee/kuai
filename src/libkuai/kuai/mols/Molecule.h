#include <boost/tuple/tuple.hpp>

#include <kuai/typedef.h>
#include <kuai/mols/Atom.h>
#include <kuai/mols/Bond.h>
#include <kuai/mols/ChiralCenter.h>
#include <kuai/mols/ConnectTable.h>


#ifndef _KUAI_MOLS_MOLECULE_H_
#define _KUAI_MOLS_MOLECULE_H_

namespace kuai {

	class Molecule 
		: public Noncopyable 
	{
	///////////////////////////////////////////////////////////////////
	// Basic functions
	public:
		explicit Molecule(AtomArray& atoms, BondArray& bonds, bool swap_array=false);
		explicit Molecule(AtomArrayPtr atoms, BondArrayPtr bonds);
		explicit Molecule(SharedPtr<Molecule> pmol, const std::vector<Atom*>& atoms, const std::vector<Bond*>& bonds);
		explicit Molecule(SharedPtr<Molecule> pmol, const std::vector<Atom*>& atoms);
		explicit Molecule(SharedPtr<Molecule> pmol, const std::vector<Bond*>& bonds);

		virtual ~Molecule();

		SharedPtr<Molecule> clone();	

	public:
		String name;						//> Name of the molecule.
		HashMap<String, String> text_info;	//> Text information readed from input.

		Index count_atoms() const {
			return _ct.count_atoms();
		}
		Index count_bonds() const {
			return _ct.count_bonds();
		}

		Atom* get_atom(Index i) {
			assert (i < count_atoms());
			return &(*_atoms)[_ct.get_atom(i)];
		}

		Bond* get_bond(Index i) {
			assert (i < count_bonds());
			return &(*_bonds)[_ct.get_bond(i)];
		}

		Atom* get_atom(const Bond* const bond, Index i) {
			Index bondID = bond - &(_bonds->front());
			assert (bondID < _bonds->size());
			Index atomID = _ct.get_atom(bondID, i);
			assert (atomID < _atoms->size());
			return &(*_atoms)[atomID];
		}
		Index degree(const Atom* const atom) const {
			Index id = atom-&(_atoms->front());
			assert (id < _atoms->size());
			return _ct.degree(id);
		}
		
		Atom* get_neighbor_atom(const Atom* const a, Index j) {
			Index id = a - &(_atoms->front());
			assert (id < _atoms->size());
			return &(*_atoms)[ _ct.get_neighbor_atom(id, j)];
		}
		Bond* get_neighbor_bond(const Atom* const a, Index j) {
			Index id = a-&(_atoms->front());
			assert (id < _atoms->size());
			return &(*_bonds)[ _ct.get_neighbor_bond(id, j)];
		}
		
		Bond* get_bond(const Atom* const atom1, const Atom* const atom2) {
			Index a1 = atom1-&(_atoms->front()), a2 = atom2-&(_atoms->front());
			assert (a1 < _atoms->size() && a2 < _atoms->size());
			Index result = _ct.get_bond(a1, a2);
			if (result != INVALID_INDEX) {
				assert (result < _bonds->size());
				return &(*_bonds)[result];
			}
			else {
				return NULL;
			}
		};
		Index index(const Atom* const a) const {
			Index delta = a - &(_atoms->front());
			if (delta < _atoms->size()) {
				return _ct.index_of_atom(delta);
			}
			else {
				return INVALID_INDEX;
			}
		};
		Index index(const Bond* const b) const {
			Index delta = b - &(_bonds->front());
			if (delta < _bonds->size()) {
				return _ct.index_of_bond(delta);
			}
			else {
				return INVALID_INDEX;
			}
		};

		bool empty() const {
			return count_atoms() == 0;
		};

		String formula() ;
		RealNumber weight() ;

	private:
		AtomArrayPtr _atoms;			//> Atom array. May be shared with other molecules.
		BondArrayPtr _bonds;			//> Bond array. May be shared with other molecules.
		ConnectTable _ct;				//> Connect table.

	private:
		void _setup_0();
		void _setup_ct(const Molecule& mol, const Array<Atom*>& atoms, const Array<Bond*>& bonds);
		void _setup(SharedPtr<Molecule> pmol, const Array<Atom*>& atoms);
		void _setup(SharedPtr<Molecule> pmol, const Array<Bond*>& bonds);

	// Basic functions
	///////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////
	// Chiral related functions
	public:

		explicit Molecule(AtomArray& atoms, BondArray& bonds, ChiralCenterArray& chirals, bool swap_array=false);
		explicit Molecule(SharedPtr<Molecule> pmol, const std::vector<Atom*>& atoms, const std::vector<Bond*>& bonds, const std::vector<ChiralCenter*>& chirals);
		explicit Molecule(SharedPtr<Molecule> pmol, const std::vector<Atom*>& atoms, const std::vector<ChiralCenter*>& chirals);
		explicit Molecule(SharedPtr<Molecule> pmol, const std::vector<Bond*>& bonds, const std::vector<ChiralCenter*>& chirals);

		Index count_chirals() const {
			return _chiral_centers.size();
		}

		ChiralCenterPtr get_chiral_center(Index i) {
			assert (_chiral_buffer.get() != NULL);
			assert (i < _chiral_centers.size());
			assert (_chiral_centers[i] < _chiral_buffer->size());
			return &((*_chiral_buffer)[_chiral_centers[i]]);
		}
		ChiralCenterPtr get_chiral_center(AtomPtr i);
		ChiralCenterPtr get_chiral_center(BondPtr i);

		AtomPtr get_atom(ChiralCenterPtr center);
		AtomPtr get_atom(ChiralCenterPtr center, Index i);
		BondPtr get_bond(ChiralCenterPtr center);
		BondPtr get_bond(ChiralCenterPtr center, Index i);

		Index index(const ChiralCenterPtr center) const;

	private:
		ChiralCenterArrayPtr _chiral_buffer;	//> Buffer to chiral center
		IndexArray _chiral_centers;				//> Chiral Center

	private:
		void _setup_chirals(SharedPtr<Molecule>& pmol);
		void _setup_chirals(const std::vector<ChiralCenter*>& chirals);


	// Chiral related functions
	///////////////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////////////
	// Molecular modifications. 
	// Maybe time consuming and maybe cause break of saved pointers.
	public:
		AtomPtr add_atom(AtomPtr atom);
		BondPtr add_bond(BondPtr bond);
		BondPtr add_bond(AtomPtr atom1, AtomPtr atom2);

		void add_atoms(String symbol, Index n = 1);
		void add_bonds(const Array<std::pair<Atom*, Atom*> >& atom_pairs, BondOrder order);

		AtomPtr drop_atom(AtomPtr atom);
		BondPtr drop_bond(BondPtr bond);


		void drop_atoms(const Array<Atom*>& atoms);
		void drop_bonds(const Array<Bond*>& bonds);

		Index drop_hydrogens();
		Index add_hydrogens(bool generate_coords=false);

		void clean();
		void clean_chirals();

	// Molecular modifications 
	///////////////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////////////
	// Buffer access. 
	// Should not be used unless you know what you are doing.
	
	public:
		boost::tuple<AtomArrayPtr, BondArrayPtr, ChiralCenterArrayPtr> buffer() {
			return boost::make_tuple(_atoms, _bonds, _chiral_buffer);
		}

		boost::tuple<Index, Index, Index> buffer_size() {
			return boost::make_tuple(_atoms->size(), _bonds->size(), _chiral_buffer.get()?_chiral_buffer->size():0);
		}

	// Buffer access. 
	///////////////////////////////////////////////////////////////////
	
	
	};

	typedef SharedPtr<Molecule> MoleculePtr;
	typedef std::vector<Molecule> MoleculeArray;
	typedef SharedPtr<MoleculeArray> MoleculeArrayPtr;

	MoleculePtr drop_hydrogen(MoleculePtr mol);
	MoleculePtr add_hydrogen(MoleculePtr mol, bool generate_coords=false);

}

#endif
