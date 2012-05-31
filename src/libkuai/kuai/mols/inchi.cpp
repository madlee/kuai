#include <kuai/tools/strtool.h>
#include <kuai/mols/inchi.h>

namespace kuai {

	namespace {
		void check_return(int code) {
			switch (code) {
			case inchi_Ret_SKIP:
				throw std::runtime_error("not used in InChI library");

			case inchi_Ret_EOF:
				throw std::runtime_error("no structural data has been provided");

			case inchi_Ret_OKAY:	/* Success; no errors or warnings */
			case inchi_Ret_WARNING: /* Success; warning(s) issued */
				return; 

			case inchi_Ret_ERROR:   
				throw std::runtime_error("no InChI has been created");

			case inchi_Ret_FATAL:
				throw std::runtime_error("Severe error: no InChI has been created (typically, memory allocation failure)");

			case inchi_Ret_BUSY:
				throw std::runtime_error("Previuos call to InChI has not returned yet");

			default:
				throw std::runtime_error("Unknown program error");
			}
		}

		void calc_inchi(std::vector<inchi_Atom>& atoms, String options, InchiOutput& output) {
			inchi_Input input = {0};
			input.atom = &atoms[0];      
			if (!options.empty()) {
				input.szOptions = const_cast<char*>(options.c_str());
			}
			input.num_atoms = atoms.size();

			int result = GetStdINCHI(&input, &output);
			check_return(result);
		}

		void convert_to_inchi_array(MoleculePtr& mol, std::vector<inchi_Atom>& atoms) 
		{
			Index nAtoms = mol->count_atoms();
			atoms.resize(nAtoms);
			memset(&atoms[0], 0, sizeof(inchi_Atom)*nAtoms);

			for (Index i = 0; i < nAtoms; ++i) {
				AtomPtr atomI = mol->get_atom(i);
				Index degree = mol->degree(atomI);
				atoms[i].x = atomI->coords.x;
				atoms[i].y = atomI->coords.y;
				atoms[i].z = atomI->coords.z;
				String symbol = atomI->symbol();
				strcpy(atoms[i].elname, symbol.c_str());
				atoms[i].num_bonds = degree;
				std::copy(atomI->num_iso_H, atomI->num_iso_H+NUM_H_ISOTOPES+1, atoms[i].num_iso_H);
				atoms[i].isotopic_mass = atomI->isotopic_mass;
				atoms[i].radical = atomI->radical;
				#ifdef USE_RATIONAL_CHARGE
					if (atomI->charge.denominator() == 1) {
						atoms[i].charge = atomI->charge.numerator();
					}
					else {
						// TODO: How to handel rational charge?
					}
				#else
					atoms[i].charge = atomI->charge;
				#endif
				for (Index j = 0; j < degree; ++j) {
					AtomPtr atomJ = mol->get_neighbor_atom(atomI, j);
					BondPtr bondJ = mol->get_neighbor_bond(atomI, j);

					atoms[i].neighbor[j] = mol->index(atomJ);
					switch (bondJ->order) {
					case SINGLE_BOND:
						atoms[i].bond_type[j] = INCHI_BOND_TYPE_SINGLE;
						break;

					case DOUBLE_BOND:
						atoms[i].bond_type[j] = INCHI_BOND_TYPE_DOUBLE;
						break;

					case TRIPLE_BOND:
						atoms[i].bond_type[j] = INCHI_BOND_TYPE_TRIPLE;
						break;

					case PARTIAL_BOND:
						atoms[i].bond_type[j] = INCHI_BOND_TYPE_ALTERN;
						break;

					default:
						atoms[i].bond_type[j] = INCHI_BOND_TYPE_NONE;
						break;
					}

					if (bondJ->stereo != BOND_STEREO_NONE) {
						if (mol->get_atom(bondJ, 0) == atomI) {
							atoms[i].bond_stereo[j] = bondJ->stereo;
						}					
						else {
							atoms[i].bond_stereo[j] = -bondJ->stereo;
						}
					}
				}
			}
		}
	}


	MoleculePtr parse_inchi(String v, String options, InchiOutputStruct& output) {
		inchi_InputINCHI input = {const_cast<char*>(v.c_str()), const_cast<char*>(options.c_str())};
		int code = GetStructFromINCHI(&input, &output );
		check_return(code);
		return output.get_result();
	}

	void inchi(MoleculePtr mol, const String& options, InchiOutput& output) {
		std::vector<inchi_Atom> atoms;
		convert_to_inchi_array(mol, atoms);
		calc_inchi(atoms, options, output);
	}

	String inchi_key(const String& inchi_code) {
		String result(28, ' ');
		int code = GetStdINCHIKeyFromStdINCHI(inchi_code.c_str(), &result[0]);
		switch (code) {
		case INCHIKEY_VALID_STANDARD:
			return result;

		case INCHIKEY_VALID_NON_STANDARD:
			throw std::runtime_error("VALID NON STANDARD");

		case INCHIKEY_INVALID_LENGTH:
			throw std::runtime_error("INVALID LENGTH");

		case INCHIKEY_INVALID_LAYOUT:
			throw std::runtime_error("INVALID LAYOUT");

		case INCHIKEY_INVALID_VERSION:
			throw std::runtime_error("INVALID VERSION");

		default:
			throw std::runtime_error("Unknown program error");
		}
	}

	MoleculePtr InchiOutputStruct::get_result() const {
		if (num_atoms) {
			std::vector<Atom> atoms; atoms.reserve(num_atoms);
			int nBonds = 0;
			for (int i = 0; i < num_atoms; ++i) {
				nBonds += atom[i].num_bonds;
			}
			std::vector<Bond> bonds; bonds.reserve(nBonds / 2);
			for (int i = 0; i < num_atoms; ++i) {
				atoms.push_back(Atom(atom[i].elname));
				Atom& atomI = atoms.back();
				atomI.coords.x = atom[i].x;
				atomI.coords.y = atom[i].y;
				atomI.coords.z = atom[i].z;

				std::copy(atom[i].num_iso_H, atom[i].num_iso_H+4, atomI.num_iso_H);
				atomI.isotopic_mass = atom[i].isotopic_mass;
				atomI.radical = Radical(atom[i].radical);
				atomI.charge = atom[i].charge;

				for (AT_NUM j = 0; j < atom[i].num_bonds; ++j) {
					if (atom[i].neighbor[j] > i) {
						bonds.push_back(Bond(i, atom[i].neighbor[j]));
						Bond& bondI = bonds.back();
						switch (atom[i].bond_type[j]) {
						case INCHI_BOND_TYPE_SINGLE:
							bondI.order = SINGLE_BOND;
							break;

						case INCHI_BOND_TYPE_DOUBLE:
							bondI.order = DOUBLE_BOND;
							break;

						case INCHI_BOND_TYPE_TRIPLE:
							bondI.order = TRIPLE_BOND;
							break;

						case INCHI_BOND_TYPE_ALTERN:
							bondI.order = PARTIAL_BOND;
							break;

						default:
							bondI.order = UNKNOWN_BOND;
							break;
						}

						if (atom[i].bond_stereo[j]>=0) {
							bondI.atom1 = i;
							bondI.atom2 = atom[i].neighbor[j];
							bondI.stereo = BondStereo(atom[i].bond_stereo[j]);
						}
						else {
							bondI.atom2 = i;
							bondI.atom1 = atom[i].neighbor[j];
							bondI.stereo = BondStereo(-atom[i].bond_stereo[j]);
						}
					}
				}
			}
			MoleculePtr result = MoleculePtr(new Molecule(atoms, bonds, true));
			return result;
		}
		else {
			return MoleculePtr();
		}
	}


	bool InchiReader::read(std::istream& stream, Mutex& mutex, const FileName* filename, size_t count, Reference data) {
		String line;

		mutex.lock();
		if (_next_inchi.empty()) {
			std::ostringstream os;
			_to_next_record(stream, os);
		}
		
		std::ostringstream os;
		os << _next_inchi << std::endl;
		bool result = _to_next_record(stream, os);
		mutex.unlock();

		data = parse_inchi(os.str(), _option);
		return result;
	}

	bool InchiReader::_to_next_record(std::istream& stream, std::ostream& os) {
		String line;
		_next_name = _next_inchi = "";
		while (getline(stream, line)) {
			if (starts_with(line, "Structure:")) {
				std::swap(_next_name, line);
			}
			else {
				os << line << std::endl;
				if (starts_with(line, "InChI=")) {
					std::swap(_next_inchi, line);
					break;
				}
			}
		}
		return !_next_inchi.empty();
	}

}
