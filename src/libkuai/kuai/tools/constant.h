#include <kuai/typedef.h>


#ifndef _KUAI_TOOLS_CONSTANT_HPP_
#define _KUAI_TOOLS_CONSTANT_HPP_

namespace kuai
{
	extern const Byte		BIT_FLAGS[8];
	extern const RealNumber PI;

	extern const RealNumber EPSILON;	// small number
	extern const RealNumber LARGE_NUMBER;

	/// degree-radian transform.
	inline RealNumber deg2rad(RealNumber v0) 
	{
		return v0 * PI / 180;
	}

	/// radian-degree transform.
	inline RealNumber rad2deg(RealNumber v0) 
	{
		return v0 / PI * 180;
	}

	namespace SI {

		extern const RealNumber LIGHT_OF_SPEED;		//	m / s	(in vacuum);
		extern const RealNumber E_CHARGE;			//	COULOMB (charge of electron);
		 

		extern const RealNumber GRAVITATIONAL_CONSTANT;		// (m*m*m) / (kg*s)

		extern const RealNumber E_MASS;			//	kg (mass of electron);
		extern const RealNumber MASS_OF_ATOM_UNIT;			// kg (mass of C12 / 12)
		extern const RealNumber MASS_OF_PROTON;			// kg
		extern const RealNumber MASS_OF_NEUTRON;			// kg

        extern const RealNumber PLANCK_CONSTANT;				// m*m*kg/s
		extern const RealNumber REDUCED_PLANCK_CONSTANT;	// m*m*kg/s (h/2/PI)
		extern const RealNumber BOHR_RADIUS;				// m
		extern const RealNumber BOLTZMANN_CONSTANT;			// J/K
		extern const RealNumber PERMITTIVITY_OF_FREE_SPACE;	// C*C / J m 

		extern const RealNumber ONE_CALORIE;				// J

		extern const RealNumber AVOGADRO_CONSTANT;

		extern const RealNumber MOLAR_GAS_CONSTANT;				// J / (mol*K)

		extern const RealNumber CELSIUS_0;							// K

		extern const RealNumber COULOMB_CONST;
	}

	extern BigIndex PRIMES_1000[];
}

#endif
