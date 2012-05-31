#include <cmath>

#include <kuai/tools/XYZ.h>
#include <kuai/tools/constant.h>

namespace kuai {

	XYZ& XYZ::operator+=(const XYZ& v1) {
		x += v1.x;
		y += v1.y;
		z += v1.z;
		return *this;
	}
	
	XYZ& XYZ::operator-=(const XYZ& v1) {
		x -= v1.x;
		y -= v1.y;
		z -= v1.z;
		return *this;
	}
	XYZ& XYZ::operator*=(const XYZ& v1) {
		*this = *this * v1;
		return *this;
	}
	XYZ& XYZ::operator*=(const RealNumber& v1) {
		x *= v1;
		y *= v1;
		z *= v1;
		return *this;
	}
	XYZ& XYZ::operator/=(const RealNumber& v1) {
		x /= v1;
		y /= v1;
		z /= v1;
		return *this;
	}
	
	XYZ XYZ::operator+(const XYZ& v2) const {
		XYZ result(*this);
		result += v2;
		return result;
	}
	XYZ XYZ::operator-(const XYZ& v2) const {
		XYZ result(*this);
		result -= v2;
		return result;
	}
	XYZ XYZ::operator*(const XYZ& v2) const {
		XYZ result(*this);
		result.x = y * v2.z - z * v2.y;
		result.y = z * v2.x - x * v2.z;
		result.z = x * v2.y - y * v2.x;
		return result;
	}
	XYZ XYZ::operator*(const RealNumber& v2) const {
		XYZ result(*this);
		result *= v2;
		return result;
	}
	XYZ XYZ::operator/(const RealNumber& v2) const {
		XYZ result(*this);
		result /= v2;
		return result;
	}
		
	RealNumber XYZ::abs() const {
		return sqrt(dot(*this, *this));
	}
	
	XYZ XYZ::operator-() const {
		XYZ result(*this);
		result.x = -x;
		result.y = -y;
		result.z = -z;
		return result;
	}
	
	bool XYZ::operator==(const XYZ& v2) const {
		return x == v2.x && y == v2.y && z == v2.z; 
	}
	bool XYZ::operator!=(const XYZ& v2) const {
		return !(*this== v2);
	}
	
	RealNumber dot(const XYZ& v1, const XYZ& v2) {
		return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	}
	
	RealNumber angle(const XYZ& v1, const XYZ& v2) {
		RealNumber r12 = sqrt(dot(v1, v1) * dot(v2, v2));
		RealNumber cosA = dot(v1, v2) / r12;
		if (cosA <= -1) {
			return -PI;
		}
		else if (cosA >= 1) {
			return 0;
		}
		else {
			return acos(cosA);
		}
	}
	RealNumber angle(const XYZ& v1, const XYZ& v2, const XYZ& v3) {
		return angle(v2-v1, v2-v3);
	}
	RealNumber torsion(const XYZ& v1, const XYZ& v2, const XYZ& v3, const XYZ& v4) {
		XYZ v12 = v1 - v2;
		XYZ axis = v2 - v3;
		XYZ v43 = v4 - v3;
		
		XYZ vt = v12 * axis;
		XYZ vu = v43 * axis;
		
		RealNumber result = angle(vt, vu);

		if (dot(v12, vu) < 0) { 
			result = -result;
		}
		return result;
	}
}
