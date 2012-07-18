#include <algorithm>

#include <kuai/typedef.h>


#ifndef _KUAI_TOOLS_COMMON_H_
#define _KUAI_TOOLS_COMMON_H_

#define ARRAY_LENGTH(a) (sizeof(a) / sizeof(a[0]))

namespace kuai {

	template<typename T>
	inline bool is_in_array(const T& value, const Array<T>& arr) {
		return std::find(arr.begin(), arr.end(), value) != arr.end();
	}

	template<typename T>
	inline bool is_in_array(const T& value, const T arr[]) {
		static const T VALUE0;
		for (Index i = 0; arr[i] != VALUE0; ++i) {
			if (arr[i] == value) {
				return true;
			}
		}
		return false;
	}

	template<typename Container>
		void unique(Container& v) 
	{
		std::sort(v.begin(), v.end());
		typename Container::iterator it = std::unique(v.begin(), v.end());
		v.erase(it, v.end());
	}

}

#endif
