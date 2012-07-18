#include <string>
#include <vector>
#include <map>
#include <set>

#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/utility.hpp> 
#include <boost/rational.hpp>
#include <boost/dynamic_bitset.hpp>


#ifndef _KUAI_TYPEDEF_H_
#define _KUAI_TYPEDEF_H_

#define Array std::vector

#define OrderedMap std::map
#define OrderedSet std::set

#define HashMap boost::unordered_map
#define HashSet boost::unordered_set

#define KuaiMap boost::unordered_map
#define KuaiSet boost::unordered_set

#define SharedPtr boost::shared_ptr


namespace kuai {

	typedef double			RealNumber;
	typedef boost::uint32_t	Index;
	typedef boost::int32_t	Integer;
	typedef boost::uint64_t	BigIndex;
	typedef boost::int64_t	BigInteger;

	typedef void*			Handle;
	typedef char			Char;
	typedef boost::int16_t	Short;

	typedef unsigned char	Byte;

	typedef std::vector<Byte> ByteArray;

	typedef std::string String;
	typedef std::vector<String> StringArray;
	typedef std::pair<String, String> StringPair;

	typedef std::map<String, String> StringMap;
	typedef std::map<String, String> StringMap;

	extern const Index INVALID_INDEX;

	template<typename T>
	inline bool isZero(const std::vector<T>& v0) {
		for (size_t i = 0; i < v0.size(); ++i) {
			if (v0[i] != 0) {
				return false;
			}
		}
		return true;
	}

	typedef std::vector<Index> IndexArray;
	typedef std::pair<size_t, size_t> SizePair;

	typedef boost::rational<Integer> Rational;

	typedef boost::noncopyable Noncopyable;

	typedef boost::dynamic_bitset<> BitSet;

}


#endif
