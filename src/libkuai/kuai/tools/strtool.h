#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <kuai/typedef.h>

#ifndef _KUAI_TOOLS_STRTOOL_H_
#define _KUAI_TOOLS_STRTOOL_H_

namespace kuai {

	template <typename T>
		String str(const T& v0)
	{
		return boost::lexical_cast<String>(v0);
	};

	template<>
	inline	String str<String>(const String& v0)
	{
		return v0;
	};

	extern const Char WHITE_SPACES[];

	const std::pair<String, String> split_pair(const String& v0, Char delim='=');
	const std::pair<String, String> split_pair(const String& v0, const Char delims[]);

	const StringArray split(const String& v0, const Char delims[] = WHITE_SPACES);
	const StringArray split(const String& v0, const size_t sizes[], const size_t n);
	const StringArray split(const String& v0, const size_t sizes[]);
	const StringArray split(const String& v0, const size_t size);

	const StringArray split_n(const String& v0, size_t n, const Char delim);
	const StringArray split_n(const String& v0, size_t n, const Char delims[]);

	class StringTemplate
	{
	public:
		class Handle
		{
		private:
			Handle(StringTemplate& vParent, const String& vName)
				: _pParent(&vParent), _s(vName)
			{ }
			Handle(StringTemplate& vParent, const Char szName[])
				: _pParent(&vParent), _s(szName)
			{ }
			
		public:
			Handle& operator=(const String& sString)
			{
				String::size_type i = _pParent->str().find(_s);
				while (i != String::npos)
				{
					_pParent->str() = _pParent->str().substr(0, i)
						+ sString + _pParent->str().substr(i+_s.size());
					i = _pParent->str().find(_s, i + sString.size());
				}
				return (*this);
			}
			Handle& operator=(const Char v[])
			{
				String vv(v);
				return (*this) = vv;
			}

		private: 
			StringTemplate* _pParent;
			String _s;
			friend class StringTemplate;
		};

	public:
		explicit StringTemplate(const Char szv0[])
			: data(szv0)
		{ }
		explicit StringTemplate(const String& v0)
			: data(v0)
		{ }

	public:
		const Char* c_str() const { 
			return data.c_str();
		}
		const String& str() const { 
			return data;
		}
		String& str() {
			return data;
		}

	public:
		StringTemplate::Handle operator[](const Char szName[])
		{
			return StringTemplate::Handle((*this), szName);
		}
		StringTemplate::Handle operator[](const String& vName)
		{
			return StringTemplate::Handle((*this), vName);
		}

	private:
		String data;
	};

	String encode_utf8(const String& v);

	class BadLexicalCast
		: public boost::bad_lexical_cast
	{
	public:
		BadLexicalCast(const String& s0, const std::type_info &target_type_arg) 
			: source(s0), bad_lexical_cast(typeid(String), target_type_arg)
		{ 
			descript = "Can not convert " + s0 + " into " + target_type_arg.name();
		}

		virtual ~BadLexicalCast() throw()
        { }

		const char* what() const throw() {
			return descript.c_str();
		}

	private:
		String source;
		String descript;
	};

	template<typename Target, typename Source>
		Target lexical_cast(Source arg)
	{
		try {
			return boost::lexical_cast<Target>(arg);
		}
		catch (boost::bad_lexical_cast&)
		{
			throw BadLexicalCast(arg, typeid(Target));
		}
	}

	template<>
	inline	bool lexical_cast<bool, const String&>(const String& item)
	{
		String v = boost::algorithm::to_upper_copy(item);
		if (v == "YES" || v == "Y" || v == "TRUE" || v == "T" || v == ".Y." || v == ".T." || v == "1" || v == "ON")  {
			return true;
		}
		else if (v == "NO" || v == "N" || v == "FALSE" || v == "F" || v == ".N." || v == ".F." || v == "0" || v == "OFF") {
			return false;
		}
		else {
			throw BadLexicalCast(item, typeid(bool));
		}
	};
	template<>
	inline	bool lexical_cast<bool, String>(String item)
	{
		return lexical_cast<bool, const String&>(item);
	};

	class MD5Coder
	{
	public:
		MD5Coder();
		explicit MD5Coder(const Byte start[], const Byte end[]);
		explicit MD5Coder(const Byte start[], const boost::uint64_t len);
		explicit MD5Coder(const String& data);

		void add(const Byte start[], const Byte end[]);
		void add(const Byte data[], boost::uint64_t len)
		{
			add(data, data+len);
		};
		void add(const Char data[], const boost::uint64_t n)
		{
			add(reinterpret_cast<const Byte*>(data), sizeof(Char)*n);
		}
		void add(const String& data)
		{
			add(data.data(), data.size());
		}

		const ByteArray digest() const;
		const String hexDigest() const;

	private:
		boost::uint64_t		_len;
		boost::uint32_t		_state[4];
		ByteArray	_stub;	// the remainds of last add.
	};

	inline const String md5(const Byte data[], const boost::uint64_t length)
	{
		MD5Coder coder(data, length);
		ByteArray result = coder.digest();
		return String(reinterpret_cast<const char*>(&result[0]), result.size());
	}
	inline const String md5(const String& data)
	{
		MD5Coder coder(data);
		ByteArray result = coder.digest();
		return String(reinterpret_cast<const char*>(&result[0]), result.size());
	}

	inline const String md5hex(const Byte data[], const boost::uint64_t length)
	{
		MD5Coder coder(data, length);
		return coder.hexDigest();
	}
	inline const String md5hex(const String& data)
	{
		MD5Coder coder(data);
		return coder.hexDigest();
	}


	const String
		encodeHEX(const Byte data[], size_t len);
	const String
		encodehex(const Byte data[], size_t len);

	template<typename Type>
		inline const String HEX(const Type v0)
	{ 
		return encodeHEX(reinterpret_cast<const Byte*>(&v0), sizeof(v0));
	}
	template<>
		inline const String HEX<String>(const String v0)
	{ 
		return encodeHEX(reinterpret_cast<const Byte*>(v0.data()), v0.size() * sizeof(Char));
	}
	template<>
		inline const String HEX<Byte>(const Byte v0)
	{
		Char buf[2] = {v0 / 16, v0 % 16};
		buf[0] =  (buf[0]<10)?(buf[0]+'0'):(buf[0]+'A'-10);
		buf[1] =  (buf[1]<10)?(buf[1]+'0'):(buf[1]+'A'-10);
		return String(buf, 2);
	}
	
	template<typename Type>
		inline const String hex(const Type v0)
	{ 
		return encodehex(reinterpret_cast<const Byte*>(&v0), sizeof(v0)/sizeof(Byte));
	}

	template<>
		inline const String hex<String>(const String v0)
	{ 
		return encodehex(reinterpret_cast<const Byte*>(v0.data()), v0.size() * sizeof(Char));
	}
	template<>
		inline const String hex<Byte>(const Byte v0)
	{ 
		Char buf[2] = {v0 / 16, v0 % 16};
		buf[0] =  (buf[0]<10)?(buf[0]+'0'):(buf[0]+'a'-10);
		buf[1] =  (buf[1]<10)?(buf[1]+'0'):(buf[1]+'a'-10);
		return String(buf, 2);
	}

	const String
		decodeHex(const Byte data[], size_t len);

	String
		hex(boost::uint16_t v0);
	String
		hex(boost::uint32_t v0);
	String
		hex(boost::uint64_t v0);
	boost::uint32_t
		hex(const Char szString[]);
	boost::uint32_t
		hex(const String& s);

	String
		HEX(boost::uint8_t v0);
	String
		HEX(boost::uint16_t v0);
	String
		HEX(boost::uint32_t v0);
	String
		HEX(boost::uint64_t v0);
	String
		encodeHEX(const Char data[], size_t len);
	String
		hex(boost::uint8_t v0);
	String
		hex(boost::uint16_t v0);
	String
		hex(boost::uint32_t v0);
	String
		hex(boost::uint64_t v0);
	boost::uint32_t
		hex(const Char szString[]);
	

	String
		oct(boost::uint32_t v0);
	String
		oct(boost::uint64_t v0);
	boost::uint32_t
		oct(const Char szString[]);
	boost::uint32_t
		oct(const String& s);

	String
		bin(boost::uint32_t v0);
	String
		bin(boost::uint64_t v0);
	boost::uint32_t
		bin(const Char szString[]);
	boost::uint32_t
		bin(const String& s);

	String
		basenumber(
			boost::uint32_t v0,
			boost::uint32_t base,
			const Char szDigital[]
		);
	String
		basenumber(
			boost::uint64_t v0,
			boost::uint64_t base,
			const Char szDigital[]
		);
	boost::uint32_t
		basenumber(
			const Char szString[],
			boost::uint64_t base,
			const Char szDigital[]
		);
	boost::uint32_t
		basenumber(
			const String& s,
			boost::uint64_t base,
			const Char szDigital[]
		);

	String 
		encodeBase64(
			const Char data[],			// The data to encode.
			boost::uint32_t length,				// The length of the data.
			const Char szCodeBase[]		// The 
		);

	String 
		encodeBase64(
			const Char data[], 
			boost::uint32_t length
		);

	String 
		decodeBase64(
			const Char data[], 
			const Char szCodeBase[]
		);

	String 
		decodeBase64(
			const Char data[]
		);

	const Char* 
		decodeBase64(
			const Char data[], 
			const Char szCodeBase[], 
			String& result
		);

	const Char* 
		decodeBase64(
			const Char data[], 
			String& result
		);

	inline String strip(const String& v0) {
		return boost::algorithm::trim_copy(v0);
	}

	template<typename PODT>
		PODT swapbytes(const PODT& v) 
	{
		int n = sizeof(v);
		PODT result = 0;
		const char* p0 = reinterpret_cast<const char*>(&v);
		char* p1 = reinterpret_cast<char*>(&result);
		for (int i = 0; i < n; ++i) {
			p1[i] = p0[n-i-1];
		}
		return result;
	}

	using namespace boost::algorithm;

}


#endif
