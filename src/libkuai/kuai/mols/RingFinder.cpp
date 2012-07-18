#include <boost/tuple/tuple.hpp>

#include <kuai/mols/RingFinder.h>
#include <kuai/mols/algo.h>
#include <kuai/mols/hybird.h>

namespace kuai {

	namespace {
		static const Index AROMATIC_ATOMS[] = {6, 7, 8, 15, 16, 33, 34};

		kuai::MoleculePtr get_ring_candidates(MoleculePtr& mol, BitSet& flags) {
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
							flags.reset(i);
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
								flags.reset(i);
							}
						}
					}
				}		
			} while(changed);

			
			std::vector<Atom*> atoms; atoms.reserve(flags.count());
			for (Index i = 0; i < nAtoms; ++i) {
				if (flags.test(i)) {
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
			static MoleculePtr find_min_ring(MoleculePtr& mol, BondPtr bond, Index max_ring_size);

		private:
			typedef KuaiMap<AtomPtr, std::pair<AtomPtr, Index> > SourceMap;
			explicit MinRingVisitor(MoleculePtr& mol, BondPtr bond, Index max_ring_size) 
				:	_mol(mol), _bond(bond), 
					_atom1(mol->get_atom(bond, 0)), 
					_atom2(mol->get_atom(bond, 1)),
					_max_ring_size(max_ring_size)
			{ 
				_source.insert(std::make_pair(_atom2, std::pair<AtomPtr, Index>(NULL, 2)));
			}

		public:
			void start(AtomPtr atom) {
				assert (_source.find(atom) != _source.end());

				SourceMap::const_iterator it = _source.find(atom);

				if (it->second.second > _max_ring_size) {
					throw this;
				}
			}

			void back(AtomPtr target, AtomPtr source) {
				if (target == _atom1 && source != _atom2) {
					assert (_source.find(source) != _source.end());
					Index depth = _source.find(source)->second.second + 1;
					_source.insert(std::make_pair(target, std::pair<AtomPtr, Index>(source, depth)));
					throw this;	// We find it.
				}
			}
			void find(AtomPtr target, AtomPtr source) { 
				assert (_source.find(target) == _source.end());
				assert (_source.find(source) != _source.end());

				Index depth = _source.find(source)->second.second + 1;
				_source.insert(std::make_pair(target, std::pair<AtomPtr, Index>(source, depth)));
			}

			MoleculePtr get_result() {
				SourceMap::const_iterator it = _source.find(_atom1);
				if (it == _source.end()) {
					return MoleculePtr();
				}
				else {
					std::vector<Atom*> atoms;
					atoms.reserve(_source.size());
					for (Atom* v = _atom1; v != _atom2; ) {
						atoms.push_back(v);
						v = _source.find(v)->second.first;
					}
					atoms.push_back(_atom2);
					std::vector<Bond*> bonds;
					bonds.reserve(atoms.size());
					bonds.push_back(_bond);
					std::vector<Atom*>::const_iterator it = atoms.begin();
					for(;;) {
						Atom* p1 = *it;
						++it;
						if (it == atoms.end()) {
							break;
						}
						else {
							Atom* p2 = *it;
							Bond* bond = _mol->get_bond(p1, p2);
							assert (bond != NULL);
							bonds.push_back(bond);
						}
					}
					return MoleculePtr(new Molecule(_mol, atoms, bonds));
				}
			}
		
		private:
			MoleculePtr _mol;
			AtomPtr _atom1, _atom2;
			BondPtr _bond;
			SourceMap _source;
			Index _max_ring_size;
		};

		MoleculePtr MinRingVisitor::find_min_ring(MoleculePtr& mol, BondPtr bond, Index max_ring_size) {
			MinRingVisitor visitor(mol, bond, max_ring_size);
			std::vector<VisitedFlag> flags(mol->count_atoms(), WHITE_FLAG);
			try {
				flags[mol->index(visitor._atom1)] = BLACK_FLAG;
				bft (mol, visitor, &flags.front(), visitor._atom2);
				return MoleculePtr();	// Return a null pointer
			}
			catch (MinRingVisitor*) {
				return visitor.get_result();
			}
		}
	}

	RingFinder::RingFinder(kuai::MoleculePtr mol, Index max_ring_rank) 
		: _bond_ranks(mol->count_bonds(), NOT_ON_RING), _mol(mol), _aromatic_tested(false)
	{
		BitSet flags(mol->count_atoms());
		flags.set();
		_setup(mol, max_ring_rank, flags);
	}

	RingFinder::RingFinder(kuai::MoleculePtr mol, BitSet& hints, Index max_ring_rank)
		: _bond_ranks(mol->count_bonds(), NOT_ON_RING), _mol(mol), _aromatic_tested(false)
	{
		_setup(mol, max_ring_rank, hints);
	}

	RingFinder::RingFinder(MoleculePtr mol, RingFinder& superset, Index max_ring_rank)
		: _bond_ranks(mol->count_bonds(), NOT_ON_RING), _mol(mol), _aromatic_tested(false)
	{
		Index nAtoms = mol->count_atoms();
		BitSet flags(nAtoms);
		for (Index i = 0; i < nAtoms; ++i) {
			AtomPtr atomI = mol->get_atom(i);
			if (superset.rank(atomI) >= 3 && superset.rank(atomI) <= max_ring_rank) {
				flags[i] = true;
			}
		}
		_setup(mol, max_ring_rank, flags);
	}


	void RingFinder::_setup(MoleculePtr& mol, Index max_ring_size, BitSet& flags) {
		MoleculePtr candidates = get_ring_candidates(mol, flags);
		if (candidates.get()) {
			std::vector<Index> rings_in_blocks;

			std::vector<kuai::MoleculePtr> blocks;
			split(candidates, blocks); 
			for (std::vector<kuai::MoleculePtr>::iterator
				it = blocks.begin(); it != blocks.end(); ++it)
			{
				Index nBonds = (**it).count_bonds();
				if ((**it).count_atoms() == nBonds) {
					if (nBonds <= max_ring_size) {
						// Single ring, put then into _rings and _blocks.
						_blocks.push_back(*it);
						_rings.push_back(*it);
						rings_in_blocks.push_back(1);
						for (Index j = 0; j < nBonds; ++j) {
							BondPtr bondJ = (**it).get_bond(j);
							_bond_ranks[mol->index(bondJ)] = nBonds;
						}
					}
					else {
						for (Index j = 0; j < nBonds; ++j) {
							BondPtr bondJ = (**it).get_bond(j);
							_bond_ranks[mol->index(bondJ)] = RingFinder::TOO_BIG;
						}
					}
				}
				else {
					bool has_chain = false;
					Index nrings = nBonds-(**it).count_atoms()+1;
					Array<MoleculePtr> rings; rings.reserve(nrings);
					for (Index j = 0; j < nBonds; ++j) {
						BondPtr bondJ = (**it).get_bond(j);
						if (_bond_ranks[mol->index(bondJ)] == NOT_ON_RING) {
							MoleculePtr p_ring = find_min_ring(*it, bondJ, max_ring_size);
							if (p_ring.get()) {
								// If we find a ring on the bondJ
								rings.push_back(p_ring);
								Index rank = p_ring->count_bonds();
								
								// Save the rank for all bonds on the ring
								for (Index k = 0; k < rank; ++k) {
									BondPtr bondK = p_ring->get_bond(k);
									Index iK = mol->index(bondK);
									if (_bond_ranks[iK] < 3 || rank < _bond_ranks[iK]) {
										_bond_ranks[iK] = rank;
									}
								}
							}
							else {
								has_chain = true;
							}
						}
					}
					if (has_chain) {
						std::vector<BondPtr> ring_bonds; ring_bonds.reserve(nBonds);
						for (Index j = 0; j < nBonds; ++j) {
							BondPtr bondJ = (**it).get_bond(j);
							Index iJ = mol->index(bondJ);
							if (_bond_ranks[iJ] >= 3) {
								ring_bonds.push_back(bondJ);
							}
						}
					
						std::vector<MoleculePtr> blockI; blockI.reserve(nrings);
						split(MoleculePtr(new Molecule(*it, ring_bonds)), blockI);
						for (std::vector<MoleculePtr>::const_iterator
							jt = blockI.begin(); jt != blockI.end(); ++jt)
						{
							Index n = 0;
							_blocks.push_back(*jt);
							for (Array<MoleculePtr>::const_iterator
								kt = rings.begin(); kt != rings.end(); ++kt)
							{
								assert ((**kt).count_bonds() >= 3);
								BondPtr ptr = (**kt).get_bond(Index(0));
								if ((**jt).index(ptr) != INVALID_INDEX) {
									_rings.push_back(*kt);
									++n;
								}
							}
							rings_in_blocks.push_back(n);
						}
					}
					else {
						_blocks.push_back(*it);
						_rings.insert(_rings.end(), rings.begin(), rings.end());
						rings_in_blocks.push_back(_rings.size());
						assert (_rings.size() == (*it)->count_bonds() - (*it)->count_atoms() + 1);
					}
				}
			}
			
			Index n = 0;
			assert (_blocks.size() == rings_in_blocks.size());
			for (Index i = 0; i < _blocks.size(); ++i) {
				_block_2_rings.insert(std::make_pair(_blocks[i], std::make_pair(_rings.begin()+n, _rings.begin()+n+rings_in_blocks[i])));
				n += rings_in_blocks[i];
			}
			assert (n == _rings.size());
		}
	}

	MoleculePtr RingFinder::find_min_ring(MoleculePtr mol, Bond* bond, Index min_ring_rank) {
		return MinRingVisitor::find_min_ring(mol, bond, min_ring_rank);
	}

	Index RingFinder::rank(AtomPtr atom) {
		assert (_mol->index(atom) != INVALID_INDEX);
		if (_atom_ranks.empty()) {
			_atom_ranks.resize(_mol->count_atoms(), NOT_ON_RING);
			for (Index i = 0; i < _atom_ranks.size(); ++i) {
				AtomPtr atomI = _mol->get_atom(i);
				Index deg = _mol->degree(atomI);
				for (Index j = 0; j < deg; ++j) {
					BondPtr bondJ = _mol->get_neighbor_bond(atomI, j);
					Index rankB = rank(bondJ);
					if (_atom_ranks[i] == NOT_ON_RING || rankB < _atom_ranks[i]) {
						_atom_ranks[i] = rankB;
					}
				}
			}
		}
		
		return _atom_ranks[_mol->index(atom)];
	}

	void RingFinder::test_aromatic(RingFinder* hints) {
		if (!_aromatic_tested) {
			Iterator it1, it2;
			for (boost::tie(it1, it2) = iter_blocks(); it1 != it2; ++it1) {
				MoleculePtr block_i = *it1;
				if (_is_aromatic_block(block_i, hints)) {
					_aromatic_blocks.push_back(block_i);
				}
				else if (block_i->count_atoms() <  block_i->count_bonds()) {
					MoleculePtr candidate = _get_aromatic_candidates(block_i, hints);
					RingFinder finder(candidate);
					Iterator jt1, jt2; 
					for (boost::tie(jt1, jt2) = finder.iter_blocks(); jt1 != jt2; ++jt1) {
						MoleculePtr block_j = *jt1;
						if (_is_aromatic_block(block_j, NULL)) {
							_aromatic_blocks.push_back(block_j);
						}
						else if (block_j->count_atoms() < block_j->count_bonds()) {
							Iterator kt1, kt2;
							std::vector<MoleculePtr> rings;
							for (boost::tie(kt1, kt2) = finder.iter_rings(block_j); kt1 != kt2; ++kt1) {
								MoleculePtr mol_k = *kt1;
								ChargeType npi = _count_pi_electrons(mol_k);
								if (npi <= 18 && npi % 4 == 2) {
									rings.push_back(mol_k);
								}
							}
							
							if (rings.size() == 1) {
								_aromatic_blocks.push_back(rings.front());
							}
							else {
								MoleculePtr v = compine_fragments(block_j, rings.begin(), rings.end());
								rings.clear();
								split(v, rings);
								_aromatic_blocks.insert(_aromatic_blocks.end(), rings.begin(), rings.end());
							}
						}
					}
				}
			}
			_aromatic_atoms.resize(_mol->count_atoms(), false);
			_aromatic_bonds.resize(_mol->count_bonds(), false);
			for (Iterator it = _aromatic_blocks.begin(); it != _aromatic_blocks.end(); ++it) {
				MoleculePtr mol_i = *it;
				for (Index j = 0; j < mol_i->count_atoms(); ++j) {
					AtomPtr atom_j = mol_i->get_atom(j);
					_aromatic_atoms.set(_mol->index(atom_j));
				}
				for (Index j = 0; j < mol_i->count_bonds(); ++j) {
					BondPtr bond_j = mol_i->get_bond(j);
					_aromatic_bonds.set(_mol->index(bond_j));
				}
			}
			_aromatic_tested = true;
		}
	}

	
	bool RingFinder::_is_aromatic_block(MoleculePtr& ring, RingFinder* hints) {
		const Index* const PEND = AROMATIC_ATOMS + ARRAY_LENGTH(AROMATIC_ATOMS);
		ChargeType npi = 0;

		Index natoms = ring->count_atoms();
		for (Index i = 0; i < natoms; ++i) {
			AtomPtr atom_i = ring->get_atom(i);

			if (hints && !hints->is_aromatic(atom_i)) {
				return false;
			}

			Index number = atom_i->number();
			const Index* const p = std::find(AROMATIC_ATOMS, PEND, number);
			if (p == PEND) {
				return false;
			}
			else if (get_hybird_state(_mol, atom_i) != HYBIRD_SP2) {
				return false;
			}

			npi -= atom_i->charge;

			if (number == 6) {
				switch (ring->get_unsaturate_degree(atom_i)) {
				case 1:
					npi += 1;
					break;

				default:
					return false;
				}				
			}
			else {
				switch (_mol->get_unsaturate_degree(atom_i)) {
				case 0:
					npi += 2;
					break;

				case 1:
					npi += 1;
					break;

				default:
					return false;
				}
			}
		}
		return (npi <= 18) && (npi % 4 == 2);
	}

	ChargeType RingFinder::_count_pi_electrons(MoleculePtr& ring) {
		ChargeType npi = 0;

		Index natoms = ring->count_atoms();
		for (Index i = 0; i < natoms; ++i) {
			AtomPtr atom_i = ring->get_atom(i);

			Index number = atom_i->number();
			npi -= atom_i->charge;

			if (number == 6) {
				npi += 1;		
			}
			else {
				switch (_mol->get_unsaturate_degree(atom_i)) {
				case 0:
					npi += 2;
					break;

				case 1:
					npi += 1;
					break;

				default:
					return false;
				}
			}
		}
		return npi;
	}

	MoleculePtr	RingFinder::_get_aromatic_candidates(MoleculePtr& ring, RingFinder* hints) {
		const Index* const PEND = AROMATIC_ATOMS + ARRAY_LENGTH(AROMATIC_ATOMS);
		Index natoms = ring->count_atoms();
		
		std::vector<AtomPtr> atoms; atoms.reserve(natoms);
		
		for (Index i = 0; i < natoms; ++i) {
			AtomPtr atom_i = ring->get_atom(i);

			if (hints && !hints->is_aromatic(atom_i)) {
				continue;
			}

			Index number = atom_i->number();
			const Index* const p = std::find(AROMATIC_ATOMS, PEND, number);
			if (p != PEND && get_hybird_state(_mol, atom_i) == HYBIRD_SP2) {
				atoms.push_back(atom_i);
			}
		}

		if (atoms.size() == natoms) {
			return ring;
		}
		else {
			return MoleculePtr(new Molecule(ring, atoms));
		}
	}


	/*
	Array<BondPtr> RingFinder::get_fused_bonds() {
		Array<BondPtr> result;
		Iterator it, itend;
		for (boost::tie(it, itend) = iter_blocks(); it != itend; ++it) {
			MoleculePtr molI = *it;
			Index nAtoms = molI->count_atoms(), nBonds = molI->count_bonds();
			if (nAtoms < nBonds) {
				Array<Index> rings(nBonds, 0);
				Iterator jt, jtend;
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
	*/

}
