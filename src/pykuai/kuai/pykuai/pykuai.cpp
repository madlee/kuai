#include <kuai/pykuai/pykuai.h>

namespace kuai {
	extern PyTypeObject* PYKUAI_TYPE_XYZ = NULL;
	extern PyTypeObject* PYKUAI_TYPE_ATOM = NULL;
	extern PyTypeObject* PYKUAI_TYPE_BOND = NULL;
	extern PyTypeObject* PYKUAI_TYPE_MOLECULE = NULL;

	PyTypeObject* import_type(PyObject* target, PyObject* source, const char* name) {
		PyObject* type = PyObject_GetAttrString(source, name);
		if (!type) {
			return NULL;
		}
		if (!PyType_Check(type)) {
			PyErr_Format(PyExc_TypeError, "%s is not a type ", name);
			return NULL;
		}

		if (PyObject_SetAttrString(target, name, type) < 0) {
			return NULL;
		}

		Py_DECREF(type);
		return (PyTypeObject*)type;
	}

	bool import_pykuai (PyObject* model) {
		PyObject* kuai_mol = PyImport_ImportModule("kuai.mol");
		if (!kuai_mol) {
			return false;
		}
		PYKUAI_TYPE_XYZ = import_type(model, kuai_mol, "XYZ");
		if (!PYKUAI_TYPE_XYZ) {
			return false;
		}

		PYKUAI_TYPE_ATOM  = import_type(model, kuai_mol, "Atom");
		if (!PYKUAI_TYPE_ATOM) {
			return false;
		}

		PYKUAI_TYPE_BOND  = import_type(model, kuai_mol, "Bond");
		if (!PYKUAI_TYPE_BOND) {
			return false;
		}

		PYKUAI_TYPE_MOLECULE  = import_type(model, kuai_mol, "Molecule");
		if (!PYKUAI_TYPE_MOLECULE) {
			return false;
		}

		Py_DECREF(kuai_mol);
		return true;
	}

	PyObject* convert_atom(const AtomPtr atom, bool convert_coords) {
		PyObject* result = PyObject_CallFunction((PyObject*)PYKUAI_TYPE_ATOM, "si", atom->symbol().c_str(), atom->charge);
		if (result == NULL) {
			return NULL;
		}
		if (convert_coords) {
			PyObject* xyz = PyObject_CallFunction((PyObject*)PYKUAI_TYPE_XYZ, "ddd", atom->coords.x, atom->coords.y, atom->coords.z);
			if (xyz == NULL) {
				Py_DECREF(result);
				return NULL;
			}
			if (PyObject_SetAttrString(result, "coords", xyz) == -1) {
				Py_DECREF(result);
				Py_DECREF(xyz);
				return NULL;
			}
			else {
				Py_DECREF(xyz);
			}
		}

		if (atom->isotopic_mass) {
			// TODO:
		}

		if (atom->radical) {
			// TODO:
		}

		return result;
	}

	Atom convert_atom(PyObject* atom, bool convert_coords) {
		PyObject* o = PyObject_GetAttrString(atom, "symbol");
		if (o == NULL) {
			throw py::PythonException();
		}
		String symbol(PyString_AsString(o));
		Py_DECREF(o);
		if (symbol.empty()) {
			throw py::PythonException();
		}

		Atom result(symbol);
		
		if (PyObject_HasAttrString(atom, "charge")) {
			o = PyObject_GetAttrString(atom, "charge");
			if (o == NULL) {
				throw py::PythonException();
			}
			int charge = PyInt_AsLong(o);
			Py_DECREF(o);
			if (PyErr_Occurred()) {
				throw py::PythonException();
			}
			result.charge = charge;
		}

		if (convert_coords && PyObject_HasAttrString(atom, "coords")) {
			o = PyObject_GetAttrString(atom, "coords");
			if (o == NULL) {
				throw py::PythonException();
			}
			// TODO:
		}

		return result;
	}

	PyObject* convert_bond(const BondPtr bond, PyObject* atom1, PyObject* atom2) {
		PyObject* result = PyObject_CallFunction((PyObject*)PYKUAI_TYPE_BOND, "OOi", atom1, atom2, int(bond->order));
		if (result == NULL) {
			return NULL;
		}

		return result;
	}

	Bond convert_bond(PyObject* bond, const KuaiMap<PyObject*, Index>& atoms) {
		char buffer[256];
		
		PyObject* o = PyObject_GetAttrString(bond, "atom1");
		if (o == NULL) {
			throw py::PythonException();
		}
		KuaiMap<PyObject*, Index>::const_iterator it_atom1 = atoms.find(o);
		Py_DECREF(o);
		if (it_atom1 == atoms.end()) {
			sprintf(buffer, "Atom 1 of Bond %d is not existed in the atom list", bond);
			throw py::PythonException(buffer);
		}

		o = PyObject_GetAttrString(bond, "atom2");
		if (o == NULL) {
			throw py::PythonException();
		}
		KuaiMap<PyObject*, Index>::const_iterator it_atom2 = atoms.find(o);
		Py_DECREF(o);
		if (it_atom2 == atoms.end()) {
			sprintf(buffer, "Atom 2 of Bond %d is not existed in the atom list", bond);
			throw py::PythonException(buffer);
		}

		o = PyObject_GetAttrString(bond, "order");
		if (o == NULL) {
			throw py::PythonException();
		}
		int order = PyInt_AsLong(o);
		Py_DECREF(o);
		if (PyErr_Occurred()) {
			throw py::PythonException();
		}

		Bond result(it_atom1->second, it_atom2->second, (BondOrder)(order));
		return result;
	}

	PyObject* convert_mol(const MoleculePtr& pmol, bool convert_coords) {
		Index nAtoms = pmol->count_atoms(), nBonds = pmol->count_bonds();
		PyObject* atoms = PyList_New(nAtoms);
		if (atoms == NULL) {
			return NULL;
		}

		for (Index i = 0; i < nAtoms; ++i) {
			PyObject* atomI = convert_atom(pmol->get_atom(i), convert_coords);
			if (atomI == NULL) {
				Py_DECREF(atoms);
			}

			PyList_SET_ITEM(atoms, i, atomI);
		}

		PyObject* bonds = PyList_New(nBonds);
		if (bonds == NULL) {
			Py_DECREF(atoms);
			return NULL;
		}
		for (Index i = 0; i < nBonds; ++i) {
			BondPtr bondI = pmol->get_bond(i);

			PyObject* atom1 = PyList_GET_ITEM(atoms, pmol->index(pmol->get_atom(bondI, 0)));
			PyObject* atom2 = PyList_GET_ITEM(atoms, pmol->index(pmol->get_atom(bondI, 1)));

			PyObject* bondII = convert_bond(bondI, atom1, atom2);
			if (bondII == NULL) {
				Py_DECREF(atoms);
				Py_DECREF(bonds);
				return NULL;
			}
			PyList_SET_ITEM(bonds, i, bondII);
		}

		PyObject* result = PyObject_CallFunction((PyObject*)PYKUAI_TYPE_MOLECULE, "OO", atoms, bonds);
		Py_DECREF(atoms);
		Py_DECREF(bonds);

		if (result == NULL) {
			return NULL;
		}
		if (!pmol->name.empty()) {
			PyObject* name = Py_BuildValue("s", pmol->name.c_str());
			if (name == NULL) {
				PyErr_Clear();
			}
			else {
				if (PyObject_SetAttrString(result, "name", name) == -1) {
					PyErr_Clear();
				}
			}
		}

		return result;
	}

	MoleculePtr convert_mol(PyObject* mol, bool convert_coords) {
		PyObject* atoms = PyObject_GetAttrString(mol, "atoms");
		if (atoms == NULL) {
			throw py::PythonException();
		}
		Py_ssize_t natoms = PySequence_Length(atoms);
		if (natoms == -1) {
			Py_DECREF(atoms);
			throw py::PythonException();
		}
		PyObject* bonds = PyObject_GetAttrString(mol, "bonds");
		if (bonds == NULL) {
			Py_DECREF(atoms);
			throw py::PythonException();
		}
		Py_ssize_t nbonds = PySequence_Length(bonds);
		if (nbonds == -1) {
			Py_DECREF(atoms);
			Py_DECREF(bonds);
			throw py::PythonException();
		}

		try {
			Array<Atom> atoms_v; atoms_v.reserve(natoms);

			KuaiMap<PyObject*, Index> atom_mappings;
			for (Index i = 0; i < natoms; ++i) {
				PyObject* atom_i = PySequence_GetItem(atoms, i);
				atom_mappings[atom_i] = i;
				if (atom_i == NULL) {
					throw py::PythonException();
				}
				try {
					atoms_v.push_back(convert_atom(atom_i, convert_coords));
					Py_DECREF(atom_i);
				}
				catch (py::PythonException&) {
					Py_DECREF(atom_i);
					throw;
				}
			}

			Array<Bond> bonds_v; bonds_v.reserve(nbonds);
			for (Index i = 0; i < nbonds; ++i) {
				PyObject* bond_i = PySequence_GetItem(bonds, i);

				if (bond_i == NULL) {
					throw py::PythonException();
				}
				try {
					bonds_v.push_back(convert_bond(bond_i, atom_mappings));
					Py_DECREF(bond_i);
				}
				catch (py::PythonException&) {
					Py_DECREF(bond_i);
					throw;
				}
			}

			// TODO:
			MoleculePtr result(new Molecule(atoms_v, bonds_v));
			return result;
		}
		catch (py::PythonException&) {
			Py_DECREF(atoms);
			Py_DECREF(bonds);
			throw;
		}
	}

}
