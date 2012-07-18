#include <kuai/pykuai/pykuai.h>
#include <kuai/mols/smiles.h>

static PyObject * parse_smiles(PyObject *self, PyObject *args);
static PyObject * get_smiles(PyObject *self, PyObject *args);
static PyObject * get_unique_smiles(PyObject *self, PyObject *args);


static PyMethodDef KuaiExtMethods[] = {
    {"parse_smiles",  parse_smiles, METH_VARARGS, "Parse SMILES to a molecule."},
	{"get_smiles",  get_smiles, METH_VARARGS, "Return SMILES of a molecule."},
	{"get_unique_smiles",  get_unique_smiles, METH_VARARGS, "Return unique SMILES of a molecule."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initkuaicore(void)
{
    PyObject *m;

    m = Py_InitModule("kuaicore", KuaiExtMethods);
    if (m == NULL)
        return;

	kuai::import_pykuai(m);
}


PyObject* parse_smiles(PyObject *self, PyObject *args) {
	const char *input;

    if (!PyArg_ParseTuple(args, "s", &input))
        return NULL;

	try {
		kuai::MoleculePtr mol = kuai::parse_smiles(input);
		return kuai::convert_mol(mol, false);
	}
	catch (std::exception& error) {
		PyErr_SetString(PyExc_ValueError, error.what());
		return NULL;
	}
}


PyObject * get_smiles(PyObject *self, PyObject *args) {
	PyObject* input;
	unsigned char aromatic = 0;
	if (!PyArg_ParseTuple(args, "O|b", &input, &aromatic))
        return NULL;

	try {
		kuai::MoleculePtr mol = kuai::convert_mol(input, false);
		kuai::String smiles = kuai::smiles(mol, aromatic);
		return Py_BuildValue("s#", smiles.c_str(), smiles.size());
	}
	catch (kuai::py::PythonException&) {
		return NULL;
	}
}

PyObject * get_unique_smiles(PyObject *self, PyObject *args) {
	PyObject* input;
	unsigned char aromatic = 0;
	if (!PyArg_ParseTuple(args, "O|b", &input, &aromatic))
        return NULL;

	try {
		kuai::MoleculePtr mol = kuai::convert_mol(input, false);
		kuai::String smiles = kuai::unique_smiles(mol, aromatic);
		return Py_BuildValue("s#", smiles.c_str(), smiles.size());
	}
	catch (kuai::py::PythonException&) {
		return NULL;
	}
}
