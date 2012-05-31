#include <kuai/typedef.h>
#include <kuai/mols/PBC.h>

#include <boost/tuple/tuple.hpp>

#ifndef _KUAI_MOLS_BOND_H_
#define _KUAI_MOLS_BOND_H_

namespace kuai {

	enum BondOrder {
		UNKNOWN_BOND = 0,
		SINGLE_BOND = 10,
		DOUBLE_BOND = 20,
		TRIPLE_BOND = 30, 
		PARTIAL_BOND = 15
	};

	enum BondStereo {
	   /* stereocenter-related; positive: the sharp end points to this atom  */
	   BOND_STEREO_NONE           =  0,
	   BOND_STEREO_SINGLE_1UP     =  1,
	   BOND_STEREO_SINGLE_1EITHER =  4,
	   BOND_STEREO_SINGLE_1DOWN   =  6,

	   /* stereobond-related */
	   BOND_STEREO_DOUBLE_EITHER  =  3 /* unknown stereobond geometry */
	};

	struct Bond {
		Bond()
		{ }
		Bond(Index v_atom1, Index v_atom2, BondOrder v_order = SINGLE_BOND, BondStereo v_stereo=BOND_STEREO_NONE)
			: atom1(v_atom1), atom2(v_atom2), order(v_order), stereo(v_stereo)
		{ }
		Index atom1, atom2;
		BondOrder order;
		BondStereo stereo;
	};

	typedef Bond* BondPtr;
	typedef std::vector<Bond> BondArray;
	typedef SharedPtr<BondArray> BondArrayPtr;

	template<typename OutputIterator>
	void near_neighbors(const PBC& pbc, Index natoms, const XYZ atoms[], RealNumber cutoff, OutputIterator it) {
		RealNumber v = cutoff / 3;
		int nx = int(pbc.a()/v+.5), ny = int(pbc.b()/v+.5), nz = int(pbc.c()/v+.5);
		RealNumber boxX = 1.0/nx, boxY = 1.0/ny, boxZ = 1.0/nz;

		KuaiMap<boost::tuple<int, int, int>, std::dequeue<Index> > points_in_box;
		for (Index i = 0; i < natoms; ++i) {
			XYZ frac = pbc.fraction(pbc.image(atoms[i]));
			int x = int(frac.x * boxX), y = int(frac.y * boxY), z = int(frac.z*boxZ);
			points_in_box[boost::make_tuple(x, y, z)].push_back(i);
		}
		for (KuaiMap<boost::tuple<int, int, int>, std::dequeue<Index> >::const_iterator
			it = points_in_box.begin(); it != points_in_box.end(); ++it)
		{

		}
	}
}

#endif
