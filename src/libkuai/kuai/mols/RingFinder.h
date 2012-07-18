#include <boost/optional.hpp>

#include <kuai/mols/Molecule.h>




#ifndef _KUAI_MOLS_RING_FINDER_H_
#define _KUAI_MOLS_RING_FINDER_H_

namespace kuai {

	class RingFinder {
	public:
		typedef Array<MoleculePtr>::iterator Iterator;

	public:
		explicit RingFinder(MoleculePtr mol, Index max_ring_rank = INVALID_INDEX);
		explicit RingFinder(MoleculePtr mol, BitSet& hints, Index max_ring_rank = INVALID_INDEX);
		explicit RingFinder(MoleculePtr mol, RingFinder& superset, Index max_ring_rank = INVALID_INDEX);

		enum {
			NOT_ON_RING = 0,
			UNCERTAIN = 1,
			TOO_BIG = 2
		};

		bool is_on_ring(BondPtr bond) const {
			return rank(bond) >= 3;
		}
		Index rank(BondPtr bond) const {
			assert (_mol->index(bond) != INVALID_INDEX);
			return _bond_ranks[_mol->index(bond)];
		}

		Index rank(AtomPtr atom);

		bool is_on_ring(AtomPtr atom) {
			return rank(atom) >= 3;
		}

		bool empty() const {
			return _blocks.empty();
		}

		static MoleculePtr find_min_ring(MoleculePtr mol, BondPtr bond, Index max_ring_rank);

		Index count_rings() const {
			return _rings.size();
		};
		Index count_rings(MoleculePtr block) const {
			KuaiMap<MoleculePtr, std::pair<Iterator, Iterator> >::const_iterator it = _block_2_rings.find(block);
			if (it != _block_2_rings.end()) {
				return it->second.second - it->second.first;
			}
			else {
				return 0;
			}
		};
		Index count_blocks() const {
			return _blocks.size();
		};

		std::pair<Iterator, Iterator>  iter_blocks() {
			return std::make_pair(_blocks.begin(), _blocks.end());
		}
		std::pair<Iterator, Iterator>  iter_rings() {
			return std::make_pair(_rings.begin(), _rings.end());
		}

		std::pair<Iterator, Iterator> iter_rings(MoleculePtr block) {
			KuaiMap<MoleculePtr, std::pair<Iterator, Iterator> >::iterator it = _block_2_rings.find(block);
			if (it == _block_2_rings.end()) {
				return std::make_pair(_rings.end(), _rings.end());
			}
			else {
				return it->second;
			}
		}

		std::pair<Iterator, Iterator> iter_aromatic_blocks() {
			if (!_aromatic_tested) {
				test_aromatic();
			}
			return std::make_pair(_aromatic_blocks.begin(), _aromatic_blocks.end());
		}

		void test_aromatic(RingFinder* hints = NULL);

		bool is_aromatic(const AtomPtr a) {
			test_aromatic();
			return _aromatic_atoms.test(_mol->index(a));
		};
		bool is_aromatic(const BondPtr a) {
			test_aromatic();
			return _aromatic_bonds.test(_mol->index(a));
		}

		bool is_aromatic(MoleculePtr& mol) {
			test_aromatic();
			Array<MoleculePtr>::const_iterator it = std::find(_aromatic_blocks.begin(), _aromatic_blocks.end(), mol);
			if (it != _rings.end()) {
				return true;
			}
			else {
				for (Index i = 0; i < mol->count_atoms(); ++i) {
					AtomPtr atom_i = mol->get_atom(i);
					if (!is_aromatic(atom_i)) {
						return false;
					}
				}
				for (Index i = 0; i < mol->count_bonds(); ++i) {
					BondPtr bond_i = mol->get_bond(i);
					if (!is_aromatic(bond_i)) {
						return false;
					}
				}
				return true;
			}
		}

		/*
		Array<BondPtr> get_fused_bonds();
		Array<AtomPtr> get_spiro_atoms();
		*/


	private:
		MoleculePtr	_mol;				// Original molecule
		Array<Index> _bond_ranks;		// rank of bonds. rank is the smallest ring size
		Array<Index> _atom_ranks;		// rank of atoms. rank is the smallest ring size

		Array<MoleculePtr> _blocks, _rings;
		KuaiMap<MoleculePtr, std::pair<Iterator, Iterator> > _block_2_rings;

		bool _aromatic_tested;
		Array<MoleculePtr> _aromatic_blocks;
		BitSet _aromatic_atoms, _aromatic_bonds;

		/*
		boost::optional<Array<Atom> > _spiro_atoms;
		boost::optional<Array<Bond> > _fused_bonds;
		*/

	private:
		void _setup(MoleculePtr& mol, Index max_ring_size, BitSet& hints);
		bool _is_aromatic_block(MoleculePtr& ring, RingFinder* hints);
		MoleculePtr	_get_aromatic_candidates(MoleculePtr& block, RingFinder* hints);
		ChargeType _count_pi_electrons(MoleculePtr& ring);
	};

}

#endif
