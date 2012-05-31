#include <cmath>

#include <kuai/typedef.h>
#include <kuai/tools/xyz.h>
#include <kuai/tools/constant.h>

#ifndef _KUAI_PBC_H_
#define _KUAI_PBC_H_

namespace kuai {

	class PBC {
	public:
		explicit PBC();
		explicit PBC(const RealNumber& a, const RealNumber& b, const RealNumber& c);
		explicit PBC(const RealNumber& a, const RealNumber& b, const RealNumber& c, const RealNumber& alpha, const RealNumber& beta, const RealNumber& gamma);

	public:
		const RealNumber a() const;
		const RealNumber b() const;
		const RealNumber c() const;

		const XYZ va() const;
		const XYZ vb() const;
		const XYZ vc() const;

		const RealNumber alpha() const;
		const RealNumber beta() const;
		const RealNumber gamma() const;

		void update(XYZ& v) const;
		const XYZ image(const XYZ& v) const;

		const XYZ fraction(const XYZ& v) const;

		RealNumber volumn() const;

		void scale_x(RealNumber factor);
		void scale_y(RealNumber factor);
		void scale_z(RealNumber factor);

		void scale(RealNumber x, RealNumber y, RealNumber z);
		void scale(RealNumber v);

		static PBC parse(const String& line);

	private:
		RealNumber ax, bx, by, cx, cy, cz;
	};

	template<>
	inline String str(const PBC& pbc) 
	{
		char buf[128];
		sprintf(buf, "PBC: %12.4f %12.4f %12.4f %6.2f %6.2f %6.2f", pbc.a(), pbc.b(), pbc.c(), pbc.alpha(), pbc.beta(), pbc.gamma());
		return buf;
	}

	inline std::istream& operator>>(std::istream& stream, PBC& pbc) {
		String line;
		getline(stream, line);	
		pbc = PBC::parse(line);
		return stream;
	}

	inline std::ostream& operator<<(std::ostream& stream, const PBC& pbc) {
		return stream << str(pbc);
	}
}

#endif
