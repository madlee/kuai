#include <algorithm>
#include <kuai/mols/ConnectTable.h>

namespace kuai  {

	// A connect table is an array of index.
	// It has following parts sequencely.
	// nA, nB, atoms, bonds, degrees, offsets, adjs, ratoms, rbonds
	// nA        is the number of atoms
	// nB        is the number of bonds
	// atoms     are the index of atom in original atom array. Each atom has one element.
	// bonds     are the index of bonds and the atoms connected by it. So each bond has three element.
	//           The 1st is the index of bond in original bond array. The 2nd and the 3rd are the 
	//           indexes of the atom connected by it in the CURRENT connect table.
	// degrees   are the number of atoms connect to current atom. Each atom has one element.
	// offsets   are the start position of neighbor list of atoms in the vArray.
	// adjs      are the packed neighbor list. It should be accessed with help of offsets. 
	// ratoms    are the reverse atom index. It record the index of atom in CURRENT ConnectTable for each 
	//           atom in the original atom array. And the INVALID_INDEX if the original atom is not 
	//           included in this ConnectTable
	// rbonds    are the reverse bond index. It record the index of atom in CURRENT ConnectTable for each 
	//           bond in the original atom array. And the INVALID_INDEX if the original bond is not 
	//           included in this ConnectTable
	// All indexes are started from 0.
	void ConnectTable::setup(const Index v[], Index max_atoms /*= INVALID_INDEX */, Index max_bonds /*= INVALID_INDEX */) {
		Index nAtoms = v[0], nBonds = v[1];
		Index nBuffer = 3 * nAtoms + 7 * nBonds + 2;
		if (max_atoms == INVALID_INDEX) {
			max_atoms = *std::max_element(v+2, v+2+nAtoms) + 1;
		}
		if (max_bonds == INVALID_INDEX) {
			max_bonds = v[2+nAtoms];
			for (const Index* p0 = v+2+nAtoms+3; p0 < v+2+nAtoms+3*nBonds; p0 += 3) {
				if (*p0 > max_bonds) {
					max_bonds = *p0;
				}
			}
			max_bonds += 1;
		}
		nBuffer += max_atoms + max_bonds;
		
		_vArray.resize(nBuffer);
		std::copy(v, v+2+nAtoms + 3*nBonds, _vArray.begin());

		_atoms = &_vArray[2];
		_bonds = _atoms + v[0];
		_degrees = _bonds + 3 * nBonds;
		_offsets = _degrees + nAtoms;
		_adjs = _offsets + nAtoms;
		_ratoms = _adjs + 4*nBonds;
		_rbonds = _ratoms + max_atoms;
		std::fill(_ratoms, _ratoms+max_atoms+max_bonds, INVALID_INDEX);

		for (Index i = 0; i < nAtoms; ++i) {
			_degrees[i] = _offsets[i] = 0;
			_ratoms[_atoms[i]] = i;
		}
		for (Index i = 0; i < nBonds; ++i) {
			_degrees[_bonds[i*3+1]] += 1;
			_degrees[_bonds[i*3+2]] += 1;
			_rbonds[_bonds[i*3]] = i;
		}
		_offsets[0] = 0;
		for (size_t i = 1; i < nAtoms; ++i) {
			_offsets[i] = _offsets[i-1] + _degrees[i-1] * 2;
			_degrees[i-1] = 0;
		}
		_degrees[nAtoms-1] = 0;
		for (size_t i = 0; i < nBonds; ++i) {
			Index atom1 = _bonds[i*3+1];
			Index atom2 = _bonds[i*3+2];

			_adjs[_offsets[atom1] + _degrees[atom1]*2] = atom2;
			_adjs[_offsets[atom1] + _degrees[atom1]*2+1] = i;
			_adjs[_offsets[atom2] + _degrees[atom2]*2] = atom1;
			_adjs[_offsets[atom2] + _degrees[atom2]*2+1] = i;

			_degrees[atom1] += 1;
			_degrees[atom2] += 1;
		}
	} ;

}
