#include <kuai/typedef.h>
#include <kuai/tools/strtool.h>
#include <kuai/tools/Phalanx.h>

#ifndef _KUAI_XYZ_H_
#define _KUAI_XYZ_H_

namespace kuai {

	struct XYZ {
		inline XYZ() 
		{ }
		inline XYZ(RealNumber x0, RealNumber y0, RealNumber z0)
			: x(x0), y(y0), z(z0)
		{ }
		inline XYZ(const XYZ& v1, const XYZ& v2)
			: x(v1.x-v2.x), y(v1.y-v2.y), z(v1.z-v2.z)
		{ }

		RealNumber x, y, z;
		
		XYZ& operator+=(const XYZ& v1);
		XYZ& operator-=(const XYZ& v1);
		XYZ& operator*=(const XYZ& v1);
		XYZ& operator*=(const RealNumber& v1);
		XYZ& operator/=(const RealNumber& v1);
		
		RealNumber abs() const;
		XYZ operator-() const;
		
		XYZ operator+(const XYZ& v2) const;
		XYZ operator-(const XYZ& v2) const;
		XYZ operator*(const XYZ& v2) const;
		XYZ operator*(const RealNumber& v2) const;
		XYZ operator/(const RealNumber& v2) const;
		
		bool operator==(const XYZ& v2) const;
		bool operator!=(const XYZ& v2) const;
	};	

	inline XYZ operator*(const RealNumber& v1, const XYZ& v2) {
		return v2 * v1;
	}
	
	inline RealNumber abs(const XYZ& v1) {
		return v1.abs();
	};
	


	RealNumber dot(const XYZ& v1, const XYZ& v2);
	
	RealNumber angle(const XYZ& v1, const XYZ& v2);
	RealNumber angle(const XYZ& v1, const XYZ& v2, const XYZ& v3);
	RealNumber torsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4);

	inline bool shorter2(const XYZ& v, const RealNumber r2, RealNumber& result) {
		result = dot(v, v);
		return result < r2;
	}

	inline bool shorter2(const XYZ& v, const RealNumber r2) {
		RealNumber result;
		return shorter2(v, r2, result);
	}

	inline bool shorter(const XYZ& v, const RealNumber r, RealNumber& result) {
		return shorter2(v, r*r, result);
	}

	inline bool shorter(const XYZ& v, const RealNumber r) {
		return shorter2(v, r*r);
	}

	inline std::ostream& operator<<(std::ostream& stream, const XYZ& xyz) {
		char buf[128];
		sprintf(buf, "%16.6f %16.6f %16.6f", xyz.x, xyz.y, xyz.z);
		return stream << buf;
	}

	inline std::istream& operator>>(std::istream& stream, XYZ& xyz) {
		return stream >> xyz.x >> xyz.y >> xyz.z;
	}

	template<>
	inline XYZ lexical_cast<XYZ, const String&>(const String& v0) {
		StringArray tokens = split(v0);
		if (tokens.size() == 3) {
			try {
				XYZ result;
				result.x = lexical_cast<RealNumber>(tokens[0]);
				result.y = lexical_cast<RealNumber>(tokens[1]);
				result.z = lexical_cast<RealNumber>(tokens[2]);
				return result;
			}
			catch (BadLexicalCast&) 
			{ }
		}
		throw BadLexicalCast(v0, typeid(XYZ));
	}

	template<>
	inline XYZ lexical_cast<XYZ, String>(String v0) {
		return lexical_cast<XYZ, const String&>(v0);
	}


	typedef Phalanx<RealNumber, 4> TranslateMatrix;

	inline XYZ operator*(const XYZ& v1, const TranslateMatrix& v2)
	{
		XYZ result;
		result.x = v1.x * v2(0, 0) + v1.y * v2(1, 0) + v1.z * v2(2, 0) + v2(3, 0);
		result.y = v1.x * v2(0, 1) + v1.y * v2(1, 1) + v1.z * v2(2, 1) + v2(3, 1);
		result.z = v1.x * v2(0, 2) + v1.y * v2(1, 2) + v1.z * v2(2, 2) + v2(3, 2);
		RealNumber d = v1.x * v2(0, 3) + v1.y * v2(1, 3) + v1.z * v2(2, 3) + v2(3, 3);
		result /= d;
		return result;
	}

	inline XYZ operator*(const TranslateMatrix& v2, const XYZ& v1)
	{
		XYZ result;
		result.x = v1.x * v2(0, 0) + v1.y * v2(0, 1) + v1.z * v2(0, 2) + v2(0, 3);
		result.y = v1.x * v2(1, 0) + v1.y * v2(1, 1) + v1.z * v2(1, 2) + v2(1, 3);
		result.z = v1.x * v2(2, 0) + v1.y * v2(2, 1) + v1.z * v2(2, 2) + v2(2, 3);
		RealNumber d = v1.x * v2(3, 0) + v1.y * v2(3, 1) + v1.z * v2(3, 2) + v2(3, 3);
		result /= d;
		return result;
	}

	inline TranslateMatrix TransMatrix(const XYZ& v0)
	{
		TranslateMatrix result(0.0);
		result(0, 0) = 1;
		result(1, 1) = 1;
		result(2, 2) = 1;
		result(3, 3) = 1;
		result(3, 0) = v0.x;
		result(3, 1) = v0.y;
		result(3, 2) = v0.z;
		return result;
	}

	inline TranslateMatrix RotateMatrix(const XYZ& vAxis, RealNumber angle)
	{
		RealNumber ll = 1/abs(vAxis);
		RealNumber a = vAxis.x * ll;
		RealNumber b = vAxis.y * ll;
		RealNumber c = vAxis.z * ll;

		RealNumber sina = sin(angle);
		RealNumber cosa = cos(angle);

		RealNumber result[] = {
			(1-a*a)*cosa + a*a, 
			c*sina + a*b*(1 - cosa), 
			-b*sina + a*c*(1 - cosa), 
			0,

			a*b*(1-cosa) - c*sina, 
			b*b+cosa*(1-b*b),
			b*c*(1-cosa) + a*sina,
			0,

			a*c*(1-cosa) + b*sina,
			b*c*(1-cosa) - a*sina,
			c*c + cosa*(1-c*c),
			0,

			0, 0, 0, 1
		};
		

		return TranslateMatrix(result);
	}

	inline TranslateMatrix RotateMatrix(const XYZ& vAxis1, const XYZ& vAxis2, RealNumber angle)
	{
        TranslateMatrix t1 = TransMatrix(-vAxis1);
		TranslateMatrix m = RotateMatrix(vAxis2 - vAxis1, angle);
		TranslateMatrix t2 = TransMatrix(vAxis1);
        TranslateMatrix result = t1;
        result *= m;
        result *= t2;
		return result;
	}
	

	inline XYZ rotate(const XYZ& v0,  const XYZ& vAxis, RealNumber angle)
	{
		return v0 * RotateMatrix(vAxis, angle);
	}

	inline XYZ rotate(const XYZ& v0,  const XYZ& vAxis1, const XYZ& vAxis2, RealNumber angle)
	{
		return rotate(v0-vAxis1, vAxis2 - vAxis1, angle)+vAxis1;
	}
	
}

#endif
