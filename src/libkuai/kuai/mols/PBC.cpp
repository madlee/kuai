#include <kuai/mols/PBC.h>


namespace kuai {

	PBC::PBC() {
		ax = by = cz = LARGE_NUMBER;
		bx = cx = cy = 0;
	}

	PBC::PBC(const RealNumber& a, const RealNumber& b, const RealNumber& c) 
	{
		ax = a;
		by = b;
		cz = c;
		bx = cx = cy = 0;
	}

	PBC::PBC(const RealNumber& a, const RealNumber& b, const RealNumber& c, const RealNumber& alpha, const RealNumber& beta, const RealNumber& gamma)
	{
		RealNumber cosA = cos(alpha), cosB = cos(beta), cosG = cos(gamma);
		RealNumber sinG = sin(gamma);
		ax = a;

		bx = b * cosG;
		by = b * sinG;

		cx = c * cosB;
		cy = c * (cosA - cosB*cosG) / sinG;
		cz = sqrt(c*c - cx*cx - cy*cy);
	}

	const RealNumber PBC::a() const {
		return ax;
	}
	const RealNumber PBC::b() const {
		return sqrt(bx*bx + by*by);
	}
	const RealNumber PBC::c() const {
		return sqrt(cx*cx + cy*cy + cz*cz);
	}

	const XYZ PBC::va() const {
		return XYZ(ax, 0, 0);
	}

	const XYZ PBC::vb() const {
		return XYZ(bx, by, 0);
	}

	const XYZ PBC::vc() const {
		return XYZ(cx, cy, cz);
	}

	const RealNumber PBC::alpha() const {
		return angle(vb(), vc());
	};
	const RealNumber PBC::beta() const {
		return acos(cx/c());
	}
	const RealNumber PBC::gamma() const {
		return atan2(by, bx);
	}

	void PBC::update(XYZ& v) const {
		while (v.z < -cz/2) {
			v.x += cx;
			v.y += cy;
			v.z += cz;
		}
		while (v.z > cz/2) {
			v.x -= cx;
			v.y -= cy;
			v.z -= cz;
		}
		while (v.y < -by/2) {
			v.x += bx;
			v.y += by;
		}
		while (v.y > by/2) {
			v.x -= bx;
			v.y -= by;
		}
		while (v.x < -ax/2) {
			v.x += ax;
		}
		while (v.x > ax/2) {
			v.x -= ax;
		}
	}

	const XYZ PBC::image(const XYZ& v) const {
		XYZ result(v);
		update(result);
		return result;
	}

	RealNumber PBC::volumn() const {
		return ax * by * cz;
	};

	void PBC::scale_x(RealNumber factor) {
		ax *= factor; bx*= factor; cx *= factor;
	};
	void PBC::scale_y(RealNumber factor) {
		by *= factor; cy *= factor;
	};
	void PBC::scale_z(RealNumber factor) {
		cz *= factor;
	};

	void PBC::scale(RealNumber x, RealNumber y, RealNumber z) {
		scale_x(x);
		scale_y(y);
		scale_z(z);
	}
	
	void PBC::scale(RealNumber v) {
		scale(v, v, v);
	}

	PBC PBC::parse(const String& v0) {
		StringArray tokens = split(v0);
		if (tokens.size() >= 6) {
			try {
				Index offset = tokens.size() - 6;
				RealNumber x = lexical_cast<RealNumber>(tokens[offset+0]);
				RealNumber y = lexical_cast<RealNumber>(tokens[offset+1]);
				RealNumber z = lexical_cast<RealNumber>(tokens[offset+2]);

				RealNumber alpha = lexical_cast<RealNumber>(tokens[offset+3]);
				RealNumber beta  = lexical_cast<RealNumber>(tokens[offset+4]);
				RealNumber gamma = lexical_cast<RealNumber>(tokens[offset+5]);

				alpha = deg2rad(alpha);
				beta  = deg2rad(beta);
				gamma = deg2rad(gamma);

				return PBC(x, y, z, alpha, beta, gamma);
			}
			catch (BadLexicalCast&) 
			{ }
		}
		throw BadLexicalCast(v0, typeid(PBC));
	}

}
