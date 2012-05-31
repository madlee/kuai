#include <cctype>
#include <boost/tuple/tuple.hpp>
#include <boost/regex.hpp>
#include <kuai/tools/common.h>
#include <kuai/mols/smiles.h>
#include <kuai/mols/algo.h>
#include <kuai/mols/AromaticFinder.h>


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

		Array<bool> aromatic_flags; aromatic_flags.reserve(v.size());
		Array<Index> left_bonds, right_bonds;
		Array<ChiralCenter> stereos; 

		char vBondStereo = ' ';

		const char* p0;
		for (p0 = v.c_str(); *p0; ++p0) {
			switch (*p0) {
			case '(':
				branches.push(atoms.size()-1);
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
						Bond b(it->second, atoms.size()-1, last_bond.order);
						if (vBondStereo == '/') {
							left_bonds.push_back(bonds.size());
						}
						else if (vBondStereo == '\\') {
							right_bonds.push_back(bonds.size());
						}
						bonds.push_back(b);
						vBondStereo = ' ';
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
				if (aromatic_flags[bondI.atom1] && aromatic_flags[bondI.atom2]) {
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

	bool SmilesReader::read(std::istream& stream, Mutex& mutex, const FileName* filename, size_t count, Reference data) {
		String line;
		mutex.lock();
		bool result = getline(stream, line);
		mutex.unlock();
		if (result) {
			data = parse_smiles(line);
		}
		return result;
	}

	namespace {

		class SmilesCreator
			: public BasicMolecularVisitor
		{
		public:
			SmilesCreator(MoleculePtr mol, bool kekule = true) 
				: _mol(mol)
			{ 
				Index nAtoms = mol->count_atoms();
				_sequence.reserve(nAtoms);
				_tags.reserve(nAtoms);
				for (Index i = 0; i < nAtoms; ++i) {
					AtomPtr atomI = mol->get_atom(i);
					_tags.push_back(get_atom_tag(atomI));
				}
				if (!kekule) {
					AromaticFinder finder(mol);
					
				}
			}

			String get_result() {
				std::ostringstream stream;
				for (Array<Index>::const_iterator
					i = _sequence.begin(); i != _sequence.end(); ++i)
				{
					stream << _tags[*i];
				}

				return stream.str();
			}

			String get_atom_tag(AtomPtr atom) {
				assert (_mol->index(atom) != INVALID_INDEX);

				static const String ORGANIC_SUBSET_SYMBOLS[] = {"C", "O", "N", "P", "S", "F", "Cl", "Br", "I", "B", ""};
				String symbol = atom->symbol();
				bool need_brace = true;
				if (is_in_array(symbol, ORGANIC_SUBSET_SYMBOLS)) {
					if (atom->charge == 0 && atom->isotopic_mass == 0 && atom->num_iso_H[0] == 0) {
						need_brace = true;
					}
				}

				if (need_brace) {
					std::ostringstream stream;
					stream << "[";
				
					if (atom->isotopic_mass != 0) {
						stream << atom->isotopic_mass;
					}

					if (true) {	// TODO: judge aromatic here.
						stream << symbol;
					}

					{
						// TODO: Add stereo center here.
					}

					if (atom->num_iso_H[0] != 0) {
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

			void finish(AtomPtr atom) {
			}
			void find(AtomPtr target, AtomPtr source) { 
				Index target_id=_mol->index(target), source_id=_mol->index(source);
				_sequence.push_back(target_id);
			}
			void back(AtomPtr target, AtomPtr source) 
			{ }

		private:
			Array<Index>		_sequence;
			Array<String>		_tags;
			MoleculePtr			_mol;

		private:
			void _dump_atom_tag(AtomPtr atom, std::ostringstream& stream) {
				
			}
		};
	}

	String smiles(MoleculePtr mol) {
		SmilesCreator creator(mol);
		dft(mol, creator);

		return creator.get_result();
	}

}
