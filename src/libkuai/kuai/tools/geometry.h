#include <ostream>
#include <kuai/tools/phalanx.h>
#include <kuai/tools/xyz.h>

#ifndef _KUAI_TOOLS_GEOMETRY_H_2007_5_22_
#define _KUAI_TOOLS_GEOMETRY_H_2007_5_22_

namespace kuai { 
	typedef Phalanx<RealNumber, 3> HessianItem;
	typedef Phalanx<RealNumber, 3> TensionType;


	inline std::ostream& operator<<(std::ostream& stream, const HessianItem& item) {
		Char buf[64];
		sprintf(buf, ("[%16.6f, %16.6f, %16.6f]\n"), item[0][0], item[0][1], item[0][2]);
		stream << buf;
		sprintf(buf, ("[%16.6f, %16.6f, %16.6f]\n"), item[1][0], item[1][1], item[1][2]);
		stream << buf;
		sprintf(buf, ("[%16.6f, %16.6f, %16.6f]\n"), item[2][0], item[2][1], item[2][2]);
		stream << buf;
		return stream;
	
	}

	RealNumber calcAtomPair(const XYZ& v1);
	RealNumber calcAtomPair(const XYZ& v1, XYZ dx[1]);
	RealNumber calcAtomPair(const XYZ& v1, XYZ dx[1], HessianItem dx2[1]);

	RealNumber calcAtomPair(const XYZ& v1, const XYZ& v2);
	RealNumber calcAtomPair(const XYZ& v1, const XYZ& v2, XYZ dx[2]);
	RealNumber calcAtomPair(const XYZ& v1, const XYZ& v2, XYZ dx[2], HessianItem dx2[4]);

	RealNumber calcCosAngle(const XYZ& v1, const XYZ& v2);
	RealNumber calcCosAngle(const XYZ& v1, const XYZ& v2, XYZ dx[2]);
	RealNumber calcCosAngle(const XYZ& v1, const XYZ& v2, XYZ dx[2], HessianItem dx2[4]);

	RealNumber calcCosAngle(const XYZ& v1, const XYZ& v2, const XYZ& v3);
	RealNumber calcCosAngle(const XYZ& v1, const XYZ& v2, const XYZ& v3, XYZ dx[3]);
	RealNumber calcCosAngle(const XYZ& v1, const XYZ& v2, const XYZ& v3, XYZ dx[3], HessianItem dx2[9]);
	
	RealNumber calcAngle(const XYZ& v1, const XYZ& v2);
	RealNumber calcAngle(const XYZ& v1, const XYZ& v2, XYZ dx[2]);
	RealNumber calcAngle(const XYZ& v1, const XYZ& v2, XYZ dx[2], HessianItem dx2[4]);
	
	RealNumber calcAngle(const XYZ& v1, const XYZ& v2, const XYZ& v3);
	RealNumber calcAngle(const XYZ& v1, const XYZ& v2, const XYZ& v3, XYZ dx[3]);
	RealNumber calcAngle(const XYZ& v1, const XYZ& v2, const XYZ& v3, XYZ dx[3], HessianItem dx2[9]);

	RealNumber calcCosTorsion(const XYZ& v1, const XYZ& v2, const XYZ& axis);
	RealNumber calcCosTorsion(const XYZ& v1, const XYZ& v2, const XYZ& axis, XYZ result[3]);
	RealNumber calcCosTorsion(const XYZ& v1, const XYZ& v2, const XYZ& axis, XYZ result[3], HessianItem dx2[9]);

	RealNumber calcTorsion(const XYZ& v1, const XYZ& v2, const XYZ& axis);
	RealNumber calcTorsion(const XYZ& v1, const XYZ& v2, const XYZ& axis, XYZ result[3]);
	RealNumber calcTorsion(const XYZ& v1, const XYZ& v2, const XYZ& axis, XYZ result[3], HessianItem dx2[9]);

	RealNumber calcCosTorsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4);
	RealNumber calcCosTorsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4, XYZ dx[4]);
	RealNumber calcCosTorsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4, XYZ dx[4], HessianItem dx2[16]);

	RealNumber calcTorsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4);
	RealNumber calcTorsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4, XYZ dx[4]);
	RealNumber calcTorsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4, XYZ dx[4], HessianItem dx2[16]);

	RealNumber calcOffPlane(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4);
	RealNumber calcOffPlane(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4, XYZ dx[4]);

	
	const HessianItem diag(RealNumber v11, RealNumber v22, RealNumber v33);
	inline const HessianItem diag(RealNumber vdiag) {
		return diag(vdiag, vdiag, vdiag);
	}
	inline const HessianItem diag(const XYZ& v0) {
		return diag(v0.x, v0.y, v0.z);
	}

	const HessianItem prod(const XYZ& v1);
	const HessianItem prod(const XYZ& v1, const XYZ& v2);
	

}

#endif
