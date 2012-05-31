#include <boost/tuple/tuple.hpp>

#include <kuai/mols/RingFinder.h>
#include <kuai/mols/algo.h>

namespace kuai {

	namespace {
		kuai::MoleculePtr get_ring_candidates(kuai::MoleculePtr& mol, std::vector<bool>& flags) {
			Index nAtoms = mol->count_atoms();

			bool changed;
			do {
				changed = false;
				for (Index i = 0; i < nAtoms; ++i) {
					if (flags[i]) {
						Atom* atomI = mol->get_atom(i);
						Index degree = mol->degree(atomI);
						if (degree < 2) {
							changed = true;
							flags[i] = false;
						}
						else {
							Index count = 0;
							for (Index j = 0; j < degree; ++j) {
								Atom* atomJ = mol->get_neighbor_atom(atomI, j);
								if (flags[mol->index(atomJ)]) {
									count += 1;
								}
							}
							if (count < 2) {
								changed = true;
								flags[i] = false;
							}
						}
					}
				}		
			} while(changed);

			
			std::vector<Atom*> atoms; atoms.reserve(nAtoms);
			for (Index i = 0; i < nAtoms; ++i) {
				if (flags[i]) {
					atoms.push_back(mol->get_atom(i));
				}
			}
			if (atoms.empty()) {
				return kuai::MoleculePtr();
			}
			else {
				return kuai::MoleculePtr(new Molecule(mol, atoms));
			}
		}

		class MinRingVisitor
			: public BasicMolecularVisitor
		{
		public:
			static MoleculePtr find_min_ring(MoleculePtr& mol, BondPtr bond);

		private:
			explicit MinRingVisitor(MoleculePtr& mol, BondPtr bond) 
				:	_mol(mol), _bond(bond), _source(mol->count_atoms(), INVALID_INDEX),
					_atom1(mol->get_atom(bond, 0)), _atom2(mol->get_atom(bond, 1))
			{ }

		public:
			void back(AtomPtr target, AtomPtr source) {
				if (target == _atom1){
					if (source != _atom2) {
						_source[_mol->index(target)] = _mol->index(source);
						throw *this;	// We find it.
					}
				}
			}
			void find(Atom* target, Atom* source) { 
				_source[_mol->index(target)] = _mol->index(source);
			}
		
		private:
			MoleculePtr _mol;
			AtomPtr _atom1, _atom2;
			BondPtr _bond;
			std::vector<Index> _source;
		};

		MoleculePtr MinRingVisitor::find_min_ring(MoleculePtr& mol, BondPtr bond) {
			MinRingVisitor visitor(mol, bond);
			std::vector<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
			try {
				flags[mol->index(visitor._atom1)] = BLACK_FLAG;
				bft (mol, visitor, &flags.front(), visitor._atom2);
				return MoleculePtr();	// Return a null pointer
			}
			catch (MinRingVisitor&) {
				std::vector<AtomPtr> atoms; atoms.reserve(mol->count_atoms());
				Index iend = mol->index(visitor._atom2);
				for (Index i = mol->index(visitor._atom1); i != iend; i = visitor._source[i]) {
					atoms.push_back(mol->get_atom(i));
				}
				atoms.push_back(visitor._atom2);
				assert (atoms.size() >= 3);

				std::vector<BondPtr> bonds; bonds.reserve(atoms.size());
				for (Index i = 1; i< atoms.size(); ++i) {
					BondPtr bondI = mol->get_bond(atoms[i-1], atoms[i]);
					assert (bondI != NULL);
					bonds.push_back(bondI);
				}
				bonds.push_back(bond);
				MoleculePtr result = MoleculePtr(new Molecule(mol, atoms, bonds));
				assert (result->count_atoms() == result->count_bonds());
				return result;
			}
		}
	}

	RingFinder::RingFinder(kuai::MoleculePtr mol, bool check_rank_for_each_bond) 
		: _bond_ranks(mol->count_bonds(), INVALID_INDEX), _mol(mol)
	{
		Array<bool> flags(mol->count_atoms(), true);
		_setup(mol, flags, check_rank_for_each_bond);
	}

	RingFinder::RingFinder(kuai::MoleculePtr mol, Array<bool>& hints, bool check_rank_for_each_bond) 
		: _bond_ranks(mol->count_bonds(), INVALID_INDEX), _mol(mol)
	{
		_setup(mol, hints, check_rank_for_each_bond);
	}

	RingFinder::RingFinder(MoleculePtr mol, RingFinder& superset, bool check_rank_for_each_bond)
		: _bond_ranks(mol->count_bonds(), INVALID_INDEX), _mol(mol)
	{
		Index nAtoms = mol->count_atoms();
		Array<bool> flags(nAtoms, true);
		for (Index i = 0; i < nAtoms; ++i) {
			AtomPtr atomI = mol->get_atom(i);
			flags[i] = (superset.rank_of_atom(atomI) != INVALID_INDEX);
		}
		_setup(mol, flags, check_rank_for_each_bond);
	}


	void RingFinder::_setup(MoleculePtr mol, Array<bool>& flags, bool check_rank_for_each_bond) {
		kuai::MoleculePtr candidates = get_ring_candidates(mol, flags);
		if (candidates.get()) {
			std::vector<kuai::MoleculePtr> blocks;
			split(candidates, blocks); 
			for (std::vector<kuai::MoleculePtr>::iterator
				it = blocks.begin(); it != blocks.end(); ++it)
			{
				if ((**it).count_atoms() == (**it).count_bonds()) {
					// Single ring, put then into _rings and _blocks.
					Index nRank = (**it).count_atoms();
					_blocks[*it].push_back(*it);
					
					for (Index j = 0; j < nRank; ++j) {
						BondPtr bondJ = (**it).get_bond(j);
						_bond_ranks[mol->index(bondJ)] = nRank;
					}
				}
				else {
					Index nBonds = (**it).count_bonds();
					Array<MoleculePtr> rings;
					if (check_rank_for_each_bond) {
						rings.reserve(nBonds);
					}
					else {
						rings.reserve(nBonds-(**it).count_bonds()+1);
					}
					for (Index j = 0; j < nBonds; ++j) {
						BondPtr bondJ = (**it).get_bond(j);
						if (check_rank_for_each_bond || _bond_ranks[mol->index(bondJ)] == INVALID_INDEX) {
							MoleculePtr p_ring = find_min_ring(*it, bondJ);
							if (p_ring.get()) {
								// If we find a ring on the bondJ
								rings.push_back(p_ring);
								Index rank = p_ring->count_bonds();
								if (check_rank_for_each_bond) {
									// Save the rank for this bond only
									Index iJ = mol->index(bondJ);
									if (rank < _bond_ranks[iJ]) {
										_bond_ranks[iJ] = rank;
									}
								}
								else {
									// Save the rank for all bonds on the ring
									for (Index k = 0; k < rank; ++k) {
										BondPtr bondK = p_ring->get_bond(k);
										Index iK = mol->index(bondK);
										if (rank < _bond_ranks[iK]) {
											_bond_ranks[iK] = rank;
										}
									}
								}
							}
						}
					}
					std::vector<BondPtr> ring_bonds; ring_bonds.reserve(nBonds);
					for (Index j = 0; j < nBonds; ++j) {
						BondPtr bondJ = (**it).get_bond(j);
						Index iJ = mol->index(bondJ);
						if (_bond_ranks[iJ] != INVALID_INDEX) {
							ring_bonds.push_back(bondJ);
						}
					}
					if (ring_bonds.size() == nBonds) {
						_blocks[*it].swap(rings);
					}
					else {
						std::vector<MoleculePtr> blockI; blockI.reserve((**it).count_bonds()-(**it).count_atoms()+1);
						split(MoleculePtr(new Molecule(*it, ring_bonds)), blockI);
						for (std::vector<MoleculePtr>::const_iterator
							jt = blockI.begin(); jt != blockI.end(); ++jt)
						{
							Array<MoleculePtr>& v = _blocks[*jt];
							for (Array<MoleculePtr>::const_iterator
								kt = rings.begin(); kt != rings.end(); ++kt)
							{
								assert ((**kt).count_bonds() >= 3);
								BondPtr ptr = (**kt).get_bond(Index(0));
								if ((**jt).index(ptr) != INVALID_INDEX) {
									v.push_back(*kt);									
								}
							}
						}
					}
				}
			}
		}
	}

	MoleculePtr RingFinder::find_min_ring(MoleculePtr mol, Bond* bond) {
		return MinRingVisitor::find_min_ring(mol, bond);
	}

	Index RingFinder::rank_of_atom(AtomPtr atom) {
		if (_atom_ranks.empty()) {
			_atom_ranks.resize(_mol->count_atoms(), INVALID_INDEX);
			for (Index i = 0; i < _atom_ranks.size(); ++i) {
				AtomPtr atomI = _mol->get_atom(i);
				Index deg = _mol->degree(atomI);
				for (Index j = 0; j < deg; ++j) {
					BondPtr bondJ = _mol->get_neighbor_bond(atomI, j);
					Index rankB = rank_of_bond(bondJ);
					if (rankB < _atom_ranks[i]) {
						_atom_ranks[i] = rankB;
					}
				}
			}
		}
		assert (_mol->index(atom) != INVALID_INDEX);
		return _atom_ranks[_mol->index(atom)];
	}


	
	Array<BondPtr> RingFinder::get_fused_bonds() {
		Array<BondPtr> result;
		RingBlockIterator it, itend;
		for (boost::tie(it, itend) = iter_blocks(); it != itend; ++it) {
			MoleculePtr molI = *it;
			Index nAtoms = molI->count_atoms(), nBonds = molI->count_bonds();
			if (nAtoms < nBonds) {
				Array<Index> rings(nBonds, 0);
				SingleRingIterator jt, jtend;
				for (boost::tie(jt, jtend) = iter_rings(molI); jt != jtend; ++jt) {
					MoleculePtr molJ = *jt;
					Index nBondJ = molJ->count_bonds();
					for (Index k = 0; k < nBondJ; ++k) {
						BondPtr bondK = molJ->get_bond(k);
						rings[molI->index(bondK)] += 1;
					}
				}
				for (Index j = 0; j < rings.size(); ++j) {
					assert (rings[j] > 0);
					if (rings[j] > 1) {
						result.push_back(molI->get_bond(j));
					}
				}
			}
		}
		return result;
	}

	Array<AtomPtr> RingFinder::get_spiro_atoms() {
		Array<AtomPtr> result;
		RingBlockIterator it, itend;
		for (boost::tie(it, itend) = iter_blocks(); it != itend; ++it) {
			MoleculePtr molI = *it;
			Index nAtoms = molI->count_atoms(), nBonds = molI->count_bonds();
			if (nAtoms < nBonds) {
				for (Index j = 0; j < nAtoms; ++j) {
					AtomPtr atomJ = molI->get_atom(j);
					if (molI->degree(atomJ) >= 4) {
						result.push_back(atomJ);
					}
				}
			}
		}
		return result;
	}

	
}
