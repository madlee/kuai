#include <kuai/typedef.h>
#include <kuai/tools/XYZ.h>
#include <kuai/mols/Element.h>

#ifndef _KUAI_MOLS_ATOM_H_
#define _KUAI_MOLS_ATOM_H_

namespace kuai {

#ifdef USE_RATIONAL_CHARGE
	typedef Rational ChargeType;
#else 
	typedef Integer ChargeType;
#endif

	/* radical definitions */
	enum Radical {
		RADICAL_NONE    = 0,
		RADICAL_SINGLET = 1,
		RADICAL_DOUBLET = 2,
		RADICAL_TRIPLET = 3
	};

	class Atom {
	public:
		explicit Atom(const Index& number = 6) 
			: coords(0, 0, 0), isotopic_mass(0), radical(RADICAL_NONE), charge(0)
		{ 
			set_element(number);
			num_iso_H[0] = INVALID_INDEX;
			num_iso_H[1] = num_iso_H[2] = num_iso_H[3] = 0;
		}
		explicit Atom(const String& symbol) 
			: _symbol(symbol), coords(0, 0, 0), isotopic_mass(0), radical(RADICAL_NONE), charge(0)
		{
			num_iso_H[0] = INVALID_INDEX;
			num_iso_H[1] = num_iso_H[2] = num_iso_H[3] = 0;
		}
		XYZ		coords;
		Index	isotopic_mass;				/* 0 => non-isotopic; isotopic mass or  */
											/* ISOTOPIC_SHIFT_FLAG + mass - (average atomic mass) */
	    Radical radical;					/* inchi_Radical */
		ChargeType charge;					/* positive or negative; 0 => no charge */

		Index num_iso_H[4];					/* implicit hydrogen atoms */
											/* [0]: number of implicit non-isotopic H
											   (exception: num_iso_H[0]=-1 means INCHI
											   adds implicit H automatically),
											   [1]: number of implicit isotopic 1H (protium),
											   [2]: number of implicit 2H (deuterium),
											   [3]: number of implicit 3H (tritium) */

		const Element* element() const {
			return Element::get(_symbol);
		}
		String symbol() const {
			return _symbol;
		}
		Index number() const {
			if (const Element* e = element()) {
				return e->number;
			}
			else {
				return INVALID_INDEX;
			}
		}
		RealNumber weight() const {
			if (const Element* e = element()) {
				return e->weight;
			}
			else {
				return 0;
			}
		}

		void set_element(Index n) {
			if (const Element* e = Element::get(n)) {
				_symbol = e->symbol;
			}
			else {
				_symbol = "?";
			}
		}

		void set_symbol(const String& v0) {
			_symbol = v0;
		}

	private:
		String	_symbol;
	};

	typedef Atom* AtomPtr;
	typedef std::vector<Atom> AtomArray;
	typedef boost::shared_ptr<AtomArray> AtomArrayPtr;
}

#endif
