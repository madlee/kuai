#include <vector>
#include <algorithm>
#include <kuai/typedef.h>

#ifndef _KUAI_TOOLS_PHALANX_H_2005_7_26_
#define _KUAI_TOOLS_PHALANX_H_2005_7_26_


namespace kuai
{
	/// an n*n Phalanx 
	template<typename TElement, size_t n>
		class Phalanx
	{
	public:
		typedef TElement value_type;
		typedef size_t size_type;
		typedef const TElement* const_iterator;
		typedef TElement* iterator;
		typedef TElement& value_reference;
		typedef const TElement& const_value_reference;
		
	public:
		explicit Phalanx(value_type v = value_type()) { 
			std::fill(_data, _data + n * n,v);
		}
		explicit Phalanx(const_iterator pV) { 
			std::copy(pV, pV + n * n, _data);
		}
		
		value_type operator()(size_type i1, size_type i2) const {
			return _data[i1*n+i2];
		}
		value_reference operator()(size_type i1, size_type i2) {
			return _data[i1*n+i2];
		}
		
		const_iterator begin() const {
			return _data;
		}
		iterator begin() {
			return _data;
		}
		const_iterator end() const {
			return _data + n * n;
		}
		iterator end() {
			return _data + n * n;
		}

		size_type size() const {
			return n * n;
		}
		size_type rank() const {
			return n;
		}

		const Phalanx operator*(const Phalanx& v0) const {
			Phalanx result;
			for (size_t i = 0; i < n; ++i) {
				for (size_t j = 0; j < n; ++j) {
					for (size_t k = 0; k < n; ++k) {
						result(i, j) += (*this)(i, k) * v0(k, j);
					}
				}
			}
			return result;
		}

		Phalanx& operator*=(const Phalanx& v0) {
            Phalanx p((*this) * v0);
            (*this) = p;
			return (*this);
		}

		Phalanx& operator*=(const_value_reference v0) {
            for (size_t i = 0; i < n*n; ++i) { 
				_data[i] *= v0;
			}
			return (*this);
		}
		const Phalanx operator*(const_value_reference v0) const {
			Phalanx result((*this));
			result *= v0;
			return result;
		}

		Phalanx& operator+=(const Phalanx& v0) {
			for (int i = 0; i < n * n; ++i) {
				_data[i] += v0._data[i];
			}
			return (*this);
		}

		Phalanx operator+(const Phalanx& v0) const {
			Phalanx result((*this));
			result += v0;
			return result;
		}

		Phalanx& operator-=(const Phalanx& v0) {
			for (int i = 0; i < n * n; ++i) {
				_data[i] -= v0._data[i];
			}
			return (*this);
		}

		Phalanx operator-(const Phalanx& v0) const {
			Phalanx result((*this));
			result -= v0;
			return result;
		}

		Phalanx& operator/=(value_type v0) {
			for (int i = 0; i < n*n; ++i) {
				_data[i] /= v0;
			}
			return (*this);
		}

		Phalanx operator/(value_type v0) const {
			Phalanx result((*this));
			result /= v0;
			return result;
		}

		bool operator==(const Phalanx& v0) const {
			for (const_iterator i = begin(), j = v0.begin(); i != end(); ++i, ++j) {
				if (*i != *j) {
					return false;
				}
			}
			return true;
		}
		bool operator!=(const Phalanx& v0) const {
			return ! (*this) == v0;
		}

		const Phalanx operator-() const { 
			Phalanx result((*this));
			for (size_t i = 0; i < n*n; ++i) {
				result._data[i] = -result._data[i];
			}
			return result;
		}
		
		const_iterator operator[](size_t i) const { 
			return &_data[i*n];
		}
		iterator operator[](size_t i) { 
			return &_data[i*n];
		}

		void trans() {
			for (size_type i = 1; i < rank(); ++i) {
				for (size_type j = 0; j < i; ++j) {
					std::swap((*this)(i, j), (*this)(j, i));
				}
			}
		}

	public:
		static Phalanx unit_phalanx(value_type v0 = value_type(1)) {
			Phalanx result;
			for (size_type i = 0; i < result.rank(); ++i) {
				result(i, i) = v0;
			}
			return result;
		}

	private: 
		TElement _data[n*n];
	};

	template<typename TElement, size_t n>
	const Phalanx<TElement, n> operator*(const TElement& v1, const Phalanx<TElement, n>& v2) {
		Phalanx<TElement, n> result(v2);
		result *= v1;
		return result;
	}

	inline RealNumber abs(const Phalanx<RealNumber, 3>& v) { 
		return v[0][0]*v[1][1]*v[2][2] + v[0][1]*v[1][2]*v[2][0]
			+ v[1][0]*v[2][1]*v[0][2] - v[0][2]*v[1][1]*v[2][0]
			- v[0][1]*v[1][0]*v[2][2] - v[0][0]*v[1][2]*v[2][1];
	}

	template<typename TElement, size_t n>
	const Phalanx<TElement, n> trans(const Phalanx<TElement, n>& v1) {
		Phalanx<TElement, n> result(v1);
		result.trans();
		return result;
	}
}


#endif // _KUAI_PHALANX_H_2005_7_26_

