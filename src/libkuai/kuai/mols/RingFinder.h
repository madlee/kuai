#include <kuai/mols/Molecule.h>


#ifndef _KUAI_MOLS_RING_FINDER_H_
#define _KUAI_MOLS_RING_FINDER_H_

namespace kuai {

	class RingFinder {
	public:
		class RingBlockIterator 
			: public std::forward_iterator_tag
		{ 
		public:
			RingBlockIterator()
			{ }
		private:
			RingBlockIterator(HashMap<MoleculePtr, Array<MoleculePtr> >::iterator it)
				: _it(it)
			{ };

		public:
			RingBlockIterator& operator++() {
				++_it;
				return *this;
			}
			RingBlockIterator operator++(int) {
				RingBlockIterator result(*this);
				++(*this);
				return result;
			}

			bool operator==(const RingBlockIterator& v0) const {
				return _it == v0._it;
			}
			bool operator!=(const RingBlockIterator& v0) const {
				return !(*this == v0);
			}

			MoleculePtr operator*() {
				return _it->first;
			}

		private:
			HashMap<MoleculePtr, Array<MoleculePtr> >::iterator _it;
			friend class RingFinder;
		}; 

		class SingleRingIterator 
			: public std::forward_iterator_tag
		{ 
		public:
			SingleRingIterator()
			{ }
		private:
			SingleRingIterator(HashMap<MoleculePtr, Array<MoleculePtr> >::iterator it, HashMap<MoleculePtr, Array<MoleculePtr> >::iterator itend)
				: _it1(it), _itend(itend)
			{ 
				if (_it1 != _itend) {
					_it2 = _it1->second.begin();
				}
			};

			SingleRingIterator(HashMap<MoleculePtr, Array<MoleculePtr> >::iterator it1, Array<MoleculePtr>::iterator it2, HashMap<MoleculePtr, Array<MoleculePtr> >::iterator itend)
				: _it1(it1), _it2(it2), _itend(itend)
			{ };

		public:
			SingleRingIterator& operator++() {
				++_it2;
				while (_it1 != _itend && _it2 == _it1->second.end()) {
					++_it1;
					_it2 = _it1->second.begin();
				}
				return *this;
			}
			SingleRingIterator operator++(int) {
				SingleRingIterator result(*this);
				++(*this);
				return result;
			}

			bool operator==(const SingleRingIterator& v0) const {
				if (_it1 == _itend) {
					return v0._it1 == v0._itend;
				}
				else {
					return _it1 == v0._it1 && _it2 == v0._it2;
				}
			}
			bool operator!=(const SingleRingIterator& v0) const {
				return !(*this == v0);
			}

			MoleculePtr operator*() {
				return *_it2;
			}

		private:
			HashMap<MoleculePtr, Array<MoleculePtr> >::iterator _it1, _itend;
			Array<MoleculePtr>::iterator _it2;
			friend class RingFinder;
		}; 

	public:
		explicit RingFinder(MoleculePtr mol, bool check_rank_for_each_bond=false);
		explicit RingFinder(MoleculePtr mol, Array<bool>& hints, bool check_rank_for_each_bond=false);
		explicit RingFinder(MoleculePtr mol, RingFinder& superset, bool check_rank_for_each_bond=false);

		bool is_on_chain(BondPtr bond) const {
			return rank_of_bond(bond) == INVALID_INDEX;
		}
		bool is_on_ring(BondPtr bond) const {
			return !is_on_chain(bond);
		}
		Index rank_of_bond(BondPtr bond) const {
			assert (_mol->index(bond) != INVALID_INDEX);
			return _bond_ranks[_mol->index(bond)];
		}

		Index rank_of_atom(AtomPtr atom);

		bool is_on_chain(AtomPtr atom) {
			return rank_of_atom(atom) == INVALID_INDEX;
		}
		bool is_on_ring(AtomPtr atom) {
			return !is_on_chain(atom);
		}

		bool empty() const {
			return _blocks.empty();
		}

		static MoleculePtr find_min_ring(MoleculePtr mol, BondPtr bond);

		Index count_rings() const {
			Index result = 0;
			for (HashMap<MoleculePtr, Array<MoleculePtr> >::const_iterator
				it = _blocks.begin(); it != _blocks.end(); ++it)
			{
				result += it->second.size();
			}
			return result;
		};
		Index count_rings(MoleculePtr block) const {
			HashMap<MoleculePtr, Array<MoleculePtr> >::const_iterator
				it = _blocks.find(block);
			assert (it != _blocks.end());
			assert (it->second.size() > 0);
			return it->second.size();
		};
		Index count_blocks() const {
			return _blocks.size();
		};

		std::pair<RingBlockIterator, RingBlockIterator>  iter_blocks() {
			return std::pair<RingBlockIterator, RingBlockIterator>(RingBlockIterator(_blocks.begin()), RingBlockIterator(_blocks.end()));
		}
		std::pair<SingleRingIterator, SingleRingIterator>  iter_rings() {
			return std::pair<SingleRingIterator, SingleRingIterator>(
				SingleRingIterator(_blocks.begin(), _blocks.end()), 
				SingleRingIterator(_blocks.end(), _blocks.end())
			);
		}

		std::pair<SingleRingIterator, SingleRingIterator> iter_rings(MoleculePtr block) {
			assert (_blocks.find(block) != _blocks.end());
			HashMap<MoleculePtr, Array<MoleculePtr> >::iterator it1 = _blocks.find(block);
			HashMap<MoleculePtr, Array<MoleculePtr> >::iterator itend = boost::next(it1);

			return std::pair<SingleRingIterator, SingleRingIterator>(
				SingleRingIterator(it1, it1->second.begin(), itend), 
				SingleRingIterator(itend, itend)
			);
		}


		Array<BondPtr> get_fused_bonds();
		Array<AtomPtr> get_spiro_atoms();


	private:
		MoleculePtr	_mol;				// Original molecule
		Array<Index> _bond_ranks;		// rank of bonds. rank is the smallest ring size
		Array<Index> _atom_ranks;		// rank of atoms. rank is the smallest ring size

		HashMap<MoleculePtr, Array<MoleculePtr> > _blocks;

	private:
		void _setup(MoleculePtr mol, Array<bool>& hints, bool check_rank_for_each_bond);
	};

}

#endif
