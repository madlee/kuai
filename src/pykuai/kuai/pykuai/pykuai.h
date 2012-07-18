#include <Python.h>
#include <kuai/mols/Molecule.h>

#ifndef _PYKUAI_H_
#define _PYKUAI_H_

namespace kuai {

	extern PyTypeObject* PYKUAI_TYPE_XYZ;
	extern PyTypeObject* PYKUAI_TYPE_ATOM;
	extern PyTypeObject* PYKUAI_TYPE_BOND;
	extern PyTypeObject* PYKUAI_TYPE_MOLECULE;

	PyTypeObject* import_type(PyObject* target, PyObject* source, const char* name);
	bool import_pykuai (PyObject* model);

	PyObject* convert_atom(const AtomPtr atom, bool convert_coords);
	Atom convert_atom(PyObject* atom, bool convert_coords);

	PyObject* convert_mol(const MoleculePtr& pmol, bool convert_coords);

	

	MoleculePtr convert_mol(PyObject* mol, bool convert_coords);

	namespace py {

		class PythonException
			: public std::runtime_error
		{ 
		public:
			PythonException()
				: std::runtime_error("")
			{ }

			explicit PythonException(const Char msg[])
				: std::runtime_error(msg)
			{ 
				PyErr_SetString(PyExc_RuntimeError, msg);
			}
		};

		class Object {
		public:
			explicit Object(PyObject* p) 
				: _p(p)
			{ }
			virtual ~Object() {
				dec();
			};

			explicit Object(const Object& p) {
				_p = p._p;
				inc();
			}

			Object& operator=(const Object& p) {
				dec();
				_p = p._p;
				inc();
			}

			void inc() {
				Py_XINCREF(_p);
			}
			void dec() {
				Py_XDECREF(_p);
			}

			PyObject* c_ptr() const {
				return _p;
			}

		private:
			PyObject* _p;
		};
	}
}


#endif
