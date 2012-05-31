
#include <kuai/typedef.h>

#ifndef _KUAI_MOLS_ELEMENT_H_
#define _KUAI_MOLS_ELEMENT_H_

namespace kuai {

	enum {
		METAL = 1,          /* definition of an element: lowest valence */
		METAL2 = 3,         /* definition of an element: lowest and next to it valence */
		IS_METAL = 3        /* metal bitmap */
	};

	class Element {
	public:
		enum {
			MIN_ATOM_CHARGE = (-2),
			MAX_ATOM_CHARGE = 2,
			NEUTRAL_STATE = -MIN_ATOM_CHARGE,
			NUM_ATOM_CHARGES = (MAX_ATOM_CHARGE - MIN_ATOM_CHARGE + 1),
			MAX_NUM_VALENCES = 5                /* max. number + 1 to provide zero termination */
		};

		Index		number;
		const char*	symbol;
		int			nAtMass;			/* Avg. atomic mass from the Periodic Chart of the Elements (Fisher cat. no. 05-702-10) */
		int			nNormAtMass;		/* Atomic mass of the most abundant isotope */
		RealNumber	weight;				/* exact mw of the most abundant isotope */
		int			type;				/* METAL or METAL2 */
		RealNumber  electro_negativity; /* Pauling electronegativity; 0 => unknown */
		bool		do_not_add_H;		/* Does not add implicit H to atoms that have bDoNotAddH != 0 */
		char		valence[NUM_ATOM_CHARGES][MAX_NUM_VALENCES];

		bool is_metal() const {
			return (type & IS_METAL) != 0;
		}

		static const Element* get(Index i);
		static const Element* get(const Char symbol[]);
		static const Element* get(const String& symbol);
	};

}


#endif
