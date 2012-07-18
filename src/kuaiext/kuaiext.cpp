#include <kuai/pykuai/pykuai.h>
#include <kuai/mols/smiles.h>

static PyObject * parse_smiles(PyObject *self, PyObject *args);


static PyMethodDef KuaiExtMethods[] = {
    {"parse_smiles",  parse_smiles, METH_VARARGS, "Parse SMILES to a molecule."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC init_pykuai(void)
{
    PyObject *m;

    m = Py_InitModule("_kuaiext", KuaiExtMethods);
    if (m == NULL)
        return;

	kuai::import_pykuai(m);
}


PyObject* parse_smiles(PyObject *self, PyObject *args) {
	const char *smiles;

    if (!PyArg_ParseTuple(args, "s", &smiles))
        return NULL;

	kuai::MoleculePtr mol = kuai::parse_smiles(smiles);
	return kuai::convert_mol(mol);
}
