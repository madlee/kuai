#include <cctype>
#include <bitset>

#include <boost/tuple/tuple.hpp>
#include <boost/regex.hpp>

#include <kuai/tools/common.h>
#include <kuai/tools/error.h>
#include <kuai/mols/smiles.h>
#include <kuai/mols/RingFinder.h>
#include <kuai/mols/algo.h>

namespace kuai {

	namespace {
		std::pair<Index, const char*> parse_ring_label(const char* p0) {
			const char* p = p0;
			const char* pend = p0;
			for (;; ++p) {
				switch (*p) {
					case '0':
					case '1':
					case '2':
					case '3':
					case '4':
					case '5':
					case '6':
					case '7':
					case '8':
					case '9': 
						break;

					default:
						pend = p;
						goto end_label;
				}
			}

		end_label:
			return std::make_pair(lexical_cast<Index>(String(p0, pend)), --p);
		}

		boost::tuples::tuple<Atom, bool, String, const char*> parse_atom(const char* p0) {
			Atom result;
			bool is_aromatic = islower(*p0);
			String stereo;
			if (*p0 == '[') {
				static const boost::regex ATOM_PATTERN("^(\\d+)?([\\u|\\l]\\l?)(@{1,2}H?)?(H\\d*)?((\\+{2,})|(\\-{2,})|((\\+|\\-)\\d*))?$");

				const char* p1 = strchr(p0, ']');
				if (p1 != NULL) {
					boost::cmatch tokens;
					if (boost::regex_match(p0+1, p1, tokens, ATOM_PATTERN)) {
						if (tokens[1].matched) {
							result.isotopic_mass = lexical_cast<Index>(tokens[1].str());
						}
						assert (tokens[2].matched);
						result.set_symbol(tokens[2].str());

						if (tokens[3].matched) {
							stereo = tokens[3].str();
						}
						if (tokens[4].matched) {
							String v = tokens[4].str();
							if (v.size() == 1) {
								result.num_iso_H[0] = 1;
							}
							else {
								result.num_iso_H[0] = lexical_cast<Index>(v.substr(1));
							}
						}
						if (tokens[6].matched) {
							result.charge = tokens[6].length();
						}
						if (tokens[7].matched) {
							result.charge = -tokens[7].length();
						}
						if (tokens[8].matched) {
							if (tokens[8].length() == 1) {
								if (*tokens[8].first == '+') {
									result.charge = 1;
								}
								else {
									assert (*tokens[8].first == '-');
									result.charge = -1;
								}
							}
							else {
								result.charge = lexical_cast<ChargeType>(tokens[8].str());
							}
						}

						p0 = p1;
					}
					else {
						throw error("Unknonw Atom token %1%", String(p0+1, p1-1));
					}
				}
				else {
					throw error("Can not find \']\' for symbol %1%", p0);
				}
			}
			else {
				String symbol;
				if (islower(*p0)) {
					symbol = toupper(*p0);
				}
				else {
					symbol = *p0;
				}
				if (symbol == "B" && *(p0+1) == 'r' || symbol == "C" && *(p0+1) == 'l' ) {
					symbol += *(++p0);
				}
				result.set_symbol(symbol);
			}
			return boost::tuples::make_tuple(result, is_aromatic, stereo, p0);
		}
	}

	MoleculePtr parse_smiles(const String& v) {
		Array<Atom> atoms; atoms.reserve(v.size());
		Array<Bond> bonds; bonds.reserve(v.size());
		std::stack<Index> branches;
		HashMap<Index, Index> ring_labels;
		Bond last_bond(INVALID_INDEX, INVALID_INDEX, UNKNOWN_BOND);

		BitSet aromatic_flags; 
		Array<Index> left_bonds, right_bonds;
		Array<ChiralCenter> stereos; 

		char vBondStereo = ' ';

		const char* p0;
		for (p0 = v.c_str(); *p0; ++p0) {
			switch (*p0) {
			case '(':
				branches.push(last_bond.atom1);
				break;

			case ')':
				if (branches.empty()) {
					throw error("Incorrect pair of brace near %1%", p0);
				}
				last_bond.atom1 = branches.top();
				branches.pop();
				break;

			case '-':
				last_bond.order = SINGLE_BOND;
				break;

			case '=':
				last_bond.order = DOUBLE_BOND;
				break;

			case '#':
				last_bond.order = TRIPLE_BOND;
				break;

			case ':':
				last_bond.order = PARTIAL_BOND;
				break;

			case '.':
				last_bond.atom1 = last_bond.atom2 = INVALID_INDEX;
				last_bond.order = UNKNOWN_BOND;
				break;

			case '\\':
			case '/':
				vBondStereo = *p0;
				break;

			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9': 
			case '%': {
					Index ring_label;
					if (*p0 == '%') {
						boost::tie(ring_label, p0) = parse_ring_label(++p0);
					}
					else {
						ring_label = *p0-'0';
					}

					HashMap<Index, Index>::iterator it = ring_labels.find(ring_label);
					if (it == ring_labels.end()) {
						ring_labels.insert(std::make_pair(ring_label, atoms.size()-1));
					}
					else {
						last_bond.atom2 = it->second;
						bonds.push_back(last_bond);
						last_bond.order = UNKNOWN_BOND;
						ring_labels.erase(it);
					}
				}
				break;

			default:
				if (isalpha(*p0) || *p0 == '[' || *p0 == '*') {
					Atom a;
					bool is_aromatic;
					String stereo;
					boost::tie(a, is_aromatic, stereo, p0) = parse_atom(p0);
					aromatic_flags.push_back(is_aromatic);

					atoms.push_back(a);
					if (!stereo.empty()) {
						ChiralCenter center;
						center.type = ChiralCenter::STEREO_TYPE_TETRAHEDRAL;
						center.center_atom = atoms.size()-1;
						center.center_atom = atoms.size()-1;

						if (last_bond.atom1 != INVALID_INDEX) {
							center.neighbor[0] = last_bond.atom1;
						}
						else {
							center.neighbor[0] = INVALID_INDEX;
						}
						for (Index j = 1; j < 4; ++j) {
							center.neighbor[j] = INVALID_INDEX;
						}
						
						stereos.push_back(center);
					}

					if (last_bond.atom1 != INVALID_INDEX) {
						last_bond.atom2 = atoms.size()-1;
						if (vBondStereo == '/') {
							left_bonds.push_back(bonds.size());
						}
						else if (vBondStereo == '\\') {
							right_bonds.push_back(bonds.size());
						}
						bonds.push_back(last_bond);
						vBondStereo = ' ';
					}
					last_bond.atom1 = atoms.size()-1;
					last_bond.order = UNKNOWN_BOND;
				}
				else {
					goto end_label;
				}
				break;
			} // switch
		}	// for


	end_label:
		if (atoms.empty()) {
			if (boost::trim_copy(v).empty()) {
				throw error("Can not parse molecule in empty line.");
			}
			else {
				throw error("Can not parse molecule in line %1%.", v);
			}
		}

		for (BondArray::iterator 
			it = bonds.begin(); it != bonds.end(); ++it)
		{
			Bond& bondI = *it;
			if (bondI.order == UNKNOWN_BOND) {
				if (aromatic_flags.test(bondI.atom1) && aromatic_flags.test(bondI.atom2)) {
					bondI.order = PARTIAL_BOND;
				}
				else {
					bondI.order = SINGLE_BOND;
				}
			}
		}
		MoleculePtr result(new Molecule(atoms, bonds, true));
		result->text_info["SMILES"] = String(v.c_str(), p0);
		result->name = boost::trim_copy(String(p0));
		return result;
	}

	namespace {

		static const String ORGANIC_SUBSET_SYMBOLS[] = {"C", "O", "N", "P", "S", "F", "Cl", "Br", "I", "B", ""};

		static String get_atom_tag(AtomPtr atom, bool aromatic) {
			String symbol = atom->symbol();
			bool need_brace = true;
			if (is_in_array(symbol, ORGANIC_SUBSET_SYMBOLS)) {
				if (atom->charge == 0 && atom->isotopic_mass == 0 && atom->num_iso_H[0] == INVALID_INDEX) {
					need_brace = false;
				}
			}
			if (aromatic) {
				boost::algorithm::to_lower(symbol);
			}

			
			if (need_brace) {
				std::ostringstream stream;
				stream << "[";
			
				if (atom->isotopic_mass != 0) {
					stream << atom->isotopic_mass;
				}
				
				stream << symbol;

				{
					// TODO: Add stereo center here.
				}

				if (atom->num_iso_H[0] > 0 && atom->num_iso_H[0] != INVALID_INDEX) {
					stream << "H";
					if (atom->num_iso_H[0] > 1) {
						stream << atom->num_iso_H[0];
					}
				}

				if (atom->charge > 0) {
					stream << "+";
					if (atom->charge > 1) {
						stream << atom->charge;
					}
				}
				else if (atom->charge < 0) {
					stream << "-";
					if (atom->charge > 1) {
						stream << -atom->charge;
					}
				}
				
				stream << "]";
				return stream.str();
			}
			else {
				return symbol;
			}
		}

		class SmilesCreator
			: public BasicMolecularVisitor
		{
		public:
			typedef boost::tuples::tuple<AtomPtr, AtomPtr, Index> BackTraceItem;
			typedef std::bitset<256> RingLabelBitSet;

			class PredItemOrder {
			public:
				explicit PredItemOrder(const SmilesCreator& marker)
					: _marker(marker)
				{ }

				bool operator()(const SmilesCreator::BackTraceItem& v1, const SmilesCreator::BackTraceItem& v2) const {
					Index i1 = _marker.order(boost::tuples::get<0>(v1));
					Index i2 = _marker.order(boost::tuples::get<0>(v2));
					
					if (i1 < i2) {
						return true;
					}
					else if (i1 > i2) {
						return false;
					}
					else {
						return _marker.order(boost::tuples::get<1>(v1)) < _marker.order(boost::tuples::get<1>(v2));
					}
				}


			private:
				const SmilesCreator& _marker;
			};

		public:
			SmilesCreator(MoleculePtr mol) 
				: _mol(mol), _orders(mol->count_atoms(), NULL), _flags(mol->count_atoms(), WHITE_FLAG)
			{
				_backtrace.reserve(mol->count_bonds());
			}
			void start(AtomPtr p) 
			{ 
				Index i = _atom2order.size();
				_atom2order[p] = i;
				_orders[i] = p;
			}
		
			void back(AtomPtr target, AtomPtr source) 
			{ 
				Index i = _mol->index(target);
				if (_flags[i] != BLACK_FLAG) {
					Index n = _backtrace.size() / 2;
					_backtrace.push_back(boost::make_tuple(target, source, n));
					_backtrace.push_back(boost::make_tuple(source, target, n));
				}
			}

			void modify_tag(MoleculePtr mol, StringArray& tags, RingFinder* rings) {
				if (!_backtrace.empty()) {
					PredItemOrder pred_order(*this);
					std::sort(_backtrace.begin(), _backtrace.end(), pred_order);

					std::vector<Index> labels(_backtrace.size()/2, INVALID_INDEX);
					RingLabelBitSet ringsets(0);
					for (Array<BackTraceItem>::const_iterator 
						it = _backtrace.begin(); it != _backtrace.end(); ++it)
					{
						if (is_open(*it)) {
							Index label = lowest_0(ringsets);
							AtomPtr atom1 = boost::tuples::get<0>(*it);
							AtomPtr atom2 = boost::tuples::get<1>(*it);
							Index label_i = boost::tuples::get<2>(*it);

							labels[label_i] = label;
							ringsets.set(label);

							String slabel = ring_label(label);
							tags[mol->index(atom1)] += slabel;
							tags[mol->index(atom2)] += bond_label(mol, mol->get_bond(atom1, atom2), rings) + slabel;
						}
						else {
							Index label_i = boost::tuples::get<2>(*it);
							assert (labels[label_i] != INVALID_INDEX);
							ringsets.reset(labels[label_i]);
						}
					}
				}
			}

			Index lowest_0(const RingLabelBitSet& v) {
				for (Index i = 0; i < v.size(); ++i) {
					if (!v.test(i)) {
						return i;
					}
				}
				throw std::runtime_error("The ring system is too complex!");
			}

			Index order(AtomPtr p) const {
				KuaiMap<AtomPtr, Index>::const_iterator it = _atom2order.find(p);
				assert (it != _atom2order.end());
				return it->second;
			}

			bool is_open(const BackTraceItem& v) const {
				Index order1 = order(boost::tuples::get<0>(v));
				Index order2 = order(boost::tuples::get<1>(v));
				assert (order1 != order2);
				return order1 < order2;
			}

			static String ring_label(Index i) {
				i += 1;
				if (i > 9) {
					return "%" + str(i);
				}
				else {
					return str(i);
				}
			}

			static String make_smiles(MoleculePtr mol, AtomPtr atom, const Array<String>& tags, RingFinder* rings, VisitedFlag flags[]) {
				Index id = mol->index(atom);
				flags[id] = GRAY_FLAG;
				String result = tags[id];
				Index deg = mol->degree(atom);
				if (deg > 0) {
					std::vector<String> branches; branches.reserve(deg);
					for (Index i = 0; i < deg; ++i) {
						AtomPtr neighbor = mol->get_neighbor_atom(atom, i);
						id = mol->index(neighbor);
						if (flags[id] == WHITE_FLAG) {
							flags[id] = GRAY_FLAG;
							branches.push_back(bond_label(mol, mol->get_bond(atom, neighbor), rings) + make_smiles(mol, neighbor, tags, rings, flags));
						}
					}

					if (!branches.empty()) {

						for (Index i = 0; i < branches.size()-1; ++i) {
							result += "(";
							result += branches[i];
							result += ")";
						}
						result += branches.back();

					}
				}
				flags[mol->index(atom)] = BLACK_FLAG;

				return result;
			}

			static String bond_label(MoleculePtr mol, BondPtr bond, RingFinder* rings) {
				if (rings) {
					switch (bond->order) {
					case SINGLE_BOND: 
						{
							AtomPtr atom1 = mol->get_atom(bond, 0);
							AtomPtr atom2 = mol->get_atom(bond, 1);
							if (!rings->is_aromatic(bond) && rings->is_aromatic(atom1) && rings->is_aromatic(atom2))  {
								return "-";
							}
							else {
								return "";
							}
						}

					case DOUBLE_BOND:
						if (rings->is_aromatic(bond)) {
							return "";
						}
						else {
							return "=";
						}

					case TRIPLE_BOND:
						return "#";

					case PARTIAL_BOND:
						if (rings->is_aromatic(bond)) {
							return "";
						}
						else {
							return ":";
						}

					default:
						return "?";
					}
				}
				else {
					switch (bond->order) {
					case SINGLE_BOND:
						return "";

					case DOUBLE_BOND:
						return "=";

					case TRIPLE_BOND:
						return "#";

					case PARTIAL_BOND:
						return ":";

					default:
						return "?";
					}
				}
			}

			VisitedFlag* flags() {
				return &_flags[0];
			}


		private:
			MoleculePtr _mol;
			Array<AtomPtr> _orders;
			KuaiMap<AtomPtr, Index> _atom2order;
			Array<BackTraceItem> _backtrace;
			Array<VisitedFlag> _flags;
		};
	}


	String smiles(MoleculePtr& mol, bool aromatic) {
		Index nAtoms = mol->count_atoms();
		StringArray tags; tags.reserve(nAtoms);
		boost::shared_ptr<RingFinder> rings;
		if (aromatic) {
			rings = boost::shared_ptr<RingFinder>(new RingFinder(mol));
		}
		for (Index i = 0; i < nAtoms; ++i) {
			AtomPtr atom_i = mol->get_atom(i);
			tags.push_back(get_atom_tag(atom_i, aromatic && rings->is_aromatic(atom_i)));
		}


		SmilesCreator creator(mol);
		dft(mol, creator, creator.flags());

		creator.modify_tag(mol, tags, rings.get());

		std::vector<VisitedFlag> flags(nAtoms, WHITE_FLAG);
		std::vector<String> result;
		for (Index i = 0; i < nAtoms; ++i) {
			if (flags[i] == WHITE_FLAG) {
				result.push_back(SmilesCreator::make_smiles(mol, mol->get_atom(i), tags, rings.get(), &flags[0]));
			}
		}

		return boost::algorithm::join(result, ".");
	}

	String unique_smiles(MoleculePtr& mol, bool aromatic) {
		return smiles(mol, aromatic);	
	}

}
