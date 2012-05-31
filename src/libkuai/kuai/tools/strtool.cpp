#include <kuai/tools/strtool.h>

namespace kuai {

	extern const Char WHITE_SPACES[] = " \t\n\r";

	const std::pair<String, String> split_pair(const String& v0, Char delim) { 
		String::size_type i = v0.find(delim);
		if (i == v0.npos) { 
			return std::make_pair(boost::trim_copy(v0), String());
		}
		else { 
			return std::make_pair(boost::trim_copy(v0.substr(0, i)), boost::trim_copy(v0.substr(i+1)));
		}
	}
	const std::pair<String, String> split_pair(const String& v0, const Char delims[]) { 
		String::size_type i = v0.find_first_of(delims);
		if (i == v0.npos) { 
			return std::make_pair(boost::trim_copy(v0), String());
		}
		else { 
			return std::make_pair(boost::trim_copy(v0.substr(0, i)), boost::trim_copy(v0.substr(i+1)));
		}
	}

	const StringArray split(const String& v0, const Char delims[] /*= WHITE_SPACES*/) 
	{
		StringArray result;
		String::size_type i = v0.find_first_not_of(delims);
		while (i != v0.npos) {
			String::size_type j = v0.find_first_of(delims, i+1);
			result.push_back(v0.substr(i, j-i));
			i = v0.find_first_not_of(delims, j);
		}
		return result;
	}

	const StringArray split(const String& v0, const size_t sizes[], const size_t n) { 
		StringArray result;
		for (size_t i = 0, ii = 0; i < n && ii < v0.size(); ++i) { 
			result.push_back(v0.substr(ii, sizes[i]));
			ii += sizes[i];
		}
		return result;
	}
	const StringArray split(const String& v0, const size_t sizes[]) { 
		size_t n = std::find(sizes, sizes+v0.size(), 0) - sizes;
		assert (n < v0.size());
		return split(v0, sizes, n);
	}
	const StringArray split(const String& v0, const size_t size) { 
		StringArray result;
		result.reserve((v0.size() - 1) / size + 1);
		for (size_t i = 0; i < v0.size(); i += size) { 
			result.push_back(v0.substr(i, size));
		}
		return result;
	}

	const StringArray split_n(const String& v0, size_t n, const Char delim) {
		// TODO: Test
		StringArray result; 
		assert (n > 1);
		if (!v0.empty()) {
			result.reserve(n);
			String::size_type i0 = 0, i1 = v0.find(delim);
			while (result.size() < n-1 && i1 != v0.npos) {
				result.push_back(v0.substr(i0, i1-i0));
				i0 = i1+1;
				i1 = v0.find(delim, i0);
			}
			
			result.push_back(v0.substr(i0));
		}
		return result;
	} 

	const StringArray split_n(const String& v0, size_t n, const Char delims[]) {
		// TODO: Test
		StringArray result; 
		assert (n > 1);
		if (!v0.empty()) {
			result.reserve(n);
			String::size_type i = v0.find_first_not_of(delims);
			while (i != v0.npos && result.size() < n-1) {
				String::size_type j = v0.find_first_of(delims, i+1);
				result.push_back(v0.substr(i, j-i));
				i = v0.find_first_not_of(delims, j);
			}
			result.push_back(v0.substr(i));
		}
		return result;
	} 

	String encode_utf8(const String& v) {
		// TODO:
		return v;
	}

	namespace 
	{
		class md5_misc_functions
		{
		public:
			typedef boost::uint32_t WORD_TYPE;	// use int32;
			enum {
				WORD_TYPE_WIDTH = 32			// size in bits.
			};

		private:
			static inline WORD_TYPE
				F(WORD_TYPE x, WORD_TYPE y, WORD_TYPE z)
			{
				return (x & y) | (~x & z);
			}
			static inline WORD_TYPE
				G(WORD_TYPE x, WORD_TYPE y, WORD_TYPE z)
			{
				return (x & z) | (y & ~z);
			}
			static inline WORD_TYPE
				H(WORD_TYPE x, WORD_TYPE y, WORD_TYPE z)
			{
				return x ^ y ^ z;
			}
			static inline WORD_TYPE
				I(WORD_TYPE x, WORD_TYPE y, WORD_TYPE z)
			{
				return y ^ (x | ~z);
			}

			static inline WORD_TYPE cycle_left_shift(WORD_TYPE v0, WORD_TYPE bits)
			{
				return (v0 << bits) | (v0 >> (WORD_TYPE_WIDTH - bits));
			}

			static inline WORD_TYPE
				FF(
					WORD_TYPE a,
					WORD_TYPE b,
					WORD_TYPE c,
					WORD_TYPE d,
					WORD_TYPE x,
					WORD_TYPE s,
					WORD_TYPE ac
				)
			{
				return cycle_left_shift(a + F(b, c, d) + x + ac, s) + b;
			}

			static inline WORD_TYPE
				GG(
					WORD_TYPE a,
					WORD_TYPE b,
					WORD_TYPE c,
					WORD_TYPE d,
					WORD_TYPE x,
					WORD_TYPE s,
					WORD_TYPE ac
				)
			{
				return cycle_left_shift(a + G(b, c, d) + x + WORD_TYPE(ac), s) + b;
			}

			static inline WORD_TYPE
				HH(
					WORD_TYPE a,
					WORD_TYPE b,
					WORD_TYPE c,
					WORD_TYPE d,
					WORD_TYPE x,
					WORD_TYPE s,
					WORD_TYPE ac
				)
			{
				return cycle_left_shift(a + H(b, c, d) + x + WORD_TYPE(ac), s) + b;
			}

			static inline WORD_TYPE
				II(
					WORD_TYPE a,
					WORD_TYPE b,
					WORD_TYPE c,
					WORD_TYPE d,
					WORD_TYPE x,
					WORD_TYPE s,
					WORD_TYPE ac
				)
			{
				return cycle_left_shift(a + I(b, c, d) + x + WORD_TYPE(ac), s) + b;
			}

		public:
			static void encode(const Byte pMessage[], WORD_TYPE state[])
			{
				static const WORD_TYPE S[4][4] =
				{
					7, 12, 17, 22,
					5, 9, 14, 20,
					4, 11, 16, 23,
					6, 10, 15, 21
				};

				const WORD_TYPE* const x = reinterpret_cast<const WORD_TYPE* const>(pMessage);

				WORD_TYPE a = state[0], b = state[1], c = state[2], d = state[3];
				/* Round 1 */
				a = FF (a, b, c, d, x[ 0], S[0][0], 0xd76aa478); /* 1 */
				d = FF (d, a, b, c, x[ 1], S[0][1], 0xe8c7b756); /* 2 */
				c = FF (c, d, a, b, x[ 2], S[0][2], 0x242070db); /* 3 */
				b = FF (b, c, d, a, x[ 3], S[0][3], 0xc1bdceee); /* 4 */
				a = FF (a, b, c, d, x[ 4], S[0][0], 0xf57c0faf); /* 5 */
				d = FF (d, a, b, c, x[ 5], S[0][1], 0x4787c62a); /* 6 */
				c = FF (c, d, a, b, x[ 6], S[0][2], 0xa8304613); /* 7 */
				b = FF (b, c, d, a, x[ 7], S[0][3], 0xfd469501); /* 8 */
				a = FF (a, b, c, d, x[ 8], S[0][0], 0x698098d8); /* 9 */
				d = FF (d, a, b, c, x[ 9], S[0][1], 0x8b44f7af); /* 10 */
				c = FF (c, d, a, b, x[10], S[0][2], 0xffff5bb1); /* 11 */
				b = FF (b, c, d, a, x[11], S[0][3], 0x895cd7be); /* 12 */
				a = FF (a, b, c, d, x[12], S[0][0], 0x6b901122); /* 13 */
				d = FF (d, a, b, c, x[13], S[0][1], 0xfd987193); /* 14 */
				c = FF (c, d, a, b, x[14], S[0][2], 0xa679438e); /* 15 */
				b = FF (b, c, d, a, x[15], S[0][3], 0x49b40821); /* 16 */

				/* Round 2 */
				a = GG (a, b, c, d, x[ 1], S[1][0], 0xf61e2562); /* 17 */
				d = GG (d, a, b, c, x[ 6], S[1][1], 0xc040b340); /* 18 */
				c = GG (c, d, a, b, x[11], S[1][2], 0x265e5a51); /* 19 */
				b = GG (b, c, d, a, x[ 0], S[1][3], 0xe9b6c7aa); /* 20 */
				a = GG (a, b, c, d, x[ 5], S[1][0], 0xd62f105d); /* 21 */
				d = GG (d, a, b, c, x[10], S[1][1], 0x02441453); /* 22 */
				c = GG (c, d, a, b, x[15], S[1][2], 0xd8a1e681); /* 23 */
				b = GG (b, c, d, a, x[ 4], S[1][3], 0xe7d3fbc8); /* 24 */
				a = GG (a, b, c, d, x[ 9], S[1][0], 0x21e1cde6); /* 25 */
				d = GG (d, a, b, c, x[14], S[1][1], 0xc33707d6); /* 26 */
				c = GG (c, d, a, b, x[ 3], S[1][2], 0xf4d50d87); /* 27 */
				b = GG (b, c, d, a, x[ 8], S[1][3], 0x455a14ed); /* 28 */
				a = GG (a, b, c, d, x[13], S[1][0], 0xa9e3e905); /* 29 */
				d = GG (d, a, b, c, x[ 2], S[1][1], 0xfcefa3f8); /* 30 */
				c = GG (c, d, a, b, x[ 7], S[1][2], 0x676f02d9); /* 31 */
				b = GG (b, c, d, a, x[12], S[1][3], 0x8d2a4c8a); /* 32 */

				/* Round 3 */
				a = HH (a, b, c, d, x[ 5], S[2][0], 0xfffa3942); /* 33 */
				d = HH (d, a, b, c, x[ 8], S[2][1], 0x8771f681); /* 34 */
				c = HH (c, d, a, b, x[11], S[2][2], 0x6d9d6122); /* 35 */
				b = HH (b, c, d, a, x[14], S[2][3], 0xfde5380c); /* 36 */
				a = HH (a, b, c, d, x[ 1], S[2][0], 0xa4beea44); /* 37 */
				d = HH (d, a, b, c, x[ 4], S[2][1], 0x4bdecfa9); /* 38 */
				c = HH (c, d, a, b, x[ 7], S[2][2], 0xf6bb4b60); /* 39 */
				b = HH (b, c, d, a, x[10], S[2][3], 0xbebfbc70); /* 40 */
				a = HH (a, b, c, d, x[13], S[2][0], 0x289b7ec6); /* 41 */
				d = HH (d, a, b, c, x[ 0], S[2][1], 0xeaa127fa); /* 42 */
				c = HH (c, d, a, b, x[ 3], S[2][2], 0xd4ef3085); /* 43 */
				b = HH (b, c, d, a, x[ 6], S[2][3], 0x04881d05); /* 44 */
				a = HH (a, b, c, d, x[ 9], S[2][0], 0xd9d4d039); /* 45 */
				d = HH (d, a, b, c, x[12], S[2][1], 0xe6db99e5); /* 46 */
				c = HH (c, d, a, b, x[15], S[2][2], 0x1fa27cf8); /* 47 */
				b = HH (b, c, d, a, x[ 2], S[2][3], 0xc4ac5665); /* 48 */

				/* Round 4 */
				a = II (a, b, c, d, x[ 0], S[3][0], 0xf4292244); /* 49 */
				d = II (d, a, b, c, x[ 7], S[3][1], 0x432aff97); /* 50 */
				c = II (c, d, a, b, x[14], S[3][2], 0xab9423a7); /* 51 */
				b = II (b, c, d, a, x[ 5], S[3][3], 0xfc93a039); /* 52 */
				a = II (a, b, c, d, x[12], S[3][0], 0x655b59c3); /* 53 */
				d = II (d, a, b, c, x[ 3], S[3][1], 0x8f0ccc92); /* 54 */
				c = II (c, d, a, b, x[10], S[3][2], 0xffeff47d); /* 55 */
				b = II (b, c, d, a, x[ 1], S[3][3], 0x85845dd1); /* 56 */
				a = II (a, b, c, d, x[ 8], S[3][0], 0x6fa87e4f); /* 57 */
				d = II (d, a, b, c, x[15], S[3][1], 0xfe2ce6e0); /* 58 */
				c = II (c, d, a, b, x[ 6], S[3][2], 0xa3014314); /* 59 */
				b = II (b, c, d, a, x[13], S[3][3], 0x4e0811a1); /* 60 */
				a = II (a, b, c, d, x[ 4], S[3][0], 0xf7537e82); /* 61 */
				d = II (d, a, b, c, x[11], S[3][1], 0xbd3af235); /* 62 */
				c = II (c, d, a, b, x[ 2], S[3][2], 0x2ad7d2bb); /* 63 */
				b = II (b, c, d, a, x[ 9], S[3][3], 0xeb86d391); /* 64 */

				state[0] += a;
				state[1] += b;
				state[2] += c;
				state[3] += d;
			}

			static inline void Init(WORD_TYPE v[4])
			{
				static const WORD_TYPE magic_number[4] = {0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476};
				memcpy(v, magic_number, sizeof(magic_number));
			}

			static inline void CopySize(Byte* pBuf, boost::uint64_t n)
			{
				boost::uint64_t v = n * 8;
				memcpy(pBuf + 56, &v, sizeof(boost::uint64_t));
			}

			static inline void Finish(Byte* buf, boost::uint64_t nLen, WORD_TYPE p[])
			{
				size_t n = nLen % 64;

				if (n < 56)
				{
					buf[n] = Byte(128);
					if (64 - n - 9 > 0)
					{
						memset(buf + n + 1, 0, 64 - n - 9);
					}
					else
					{
						assert(64 - n - 9 == 0);
					}
					md5_misc_functions::CopySize(buf, nLen);
					md5_misc_functions::encode(buf, p);
				}
				else
				{
					assert(56 <= n && n < 64);

					buf[n] = Byte(128);
					if (64 - n - 1 > 0)
					{
						memset(buf + n + 1, 0, 64 - n - 1);
					}
					else
					{
						assert(64 - n - 1 == 0);
					}
					md5_misc_functions::encode(buf, p);

					memset(buf, 0, 56);
					md5_misc_functions::CopySize(buf, nLen);
					md5_misc_functions::encode(buf, p);
				}
			}
		};
	}	// namespace

	MD5Coder::MD5Coder()
	{
		_len = 0;
		md5_misc_functions::Init(_state);
	};
	MD5Coder::MD5Coder(const Byte start[], const Byte end[])
	{
		_len = 0;
		md5_misc_functions::Init(_state);
		add(start, end);
	}
	MD5Coder::MD5Coder(const Byte start[], const boost::uint64_t len)
	{
		_len = 0;
		md5_misc_functions::Init(_state);
		add(start, len);
	}
	MD5Coder::MD5Coder(const String& data)
	{
		_len = 0;
		md5_misc_functions::Init(_state);
		add(data);
	};
		
	void MD5Coder::add(const Byte start[], const Byte end[])
	{
		const Byte* p = start;
		const Byte* p2 = end;
		if (!_stub.empty())
		{
			size_t n = size_t(64 - _stub.size());
			if (p+n < end)
			{
				_stub.insert(_stub.end(), p, p + n);
                assert(_stub.size() == 64);
				md5_misc_functions::encode(&(_stub[0]), _state);
				p += n;
				_len += n;
			}
			else
			{
				
				_stub.insert(_stub.end(), p, p2);
				assert(_stub.size() < 64);
				p = end;
				_len += end-p;
			}
		}

		for (; p + 64 < end; p +=64)
		{
            md5_misc_functions::encode(p, _state);
			_len += 64;
		}
		_stub = ByteArray(p, p2);
		_len += _stub.size();
	};

	const ByteArray MD5Coder::digest() const
	{
		md5_misc_functions::WORD_TYPE v[4];
		std::copy(_state, _state+4, v);
		Byte buffer[64];
		std::copy(_stub.begin(), _stub.end(), buffer);
		md5_misc_functions::Finish(buffer, _len, v);

		const Byte* p1 = reinterpret_cast<const Byte*>(v);
		const Byte* p2 = p1 + sizeof(v);
		return ByteArray(p1, p2);
	};
	const String MD5Coder::hexDigest() const
	{
		ByteArray result = digest();
		return encodeHEX(&(result[0]), result.size());
	}

	const String
		encodeHEX(const Byte data[], size_t len)
	{ 
		String result;
		result.reserve(len * 2);
		for (size_t i = 0; i < len; ++i)
		{
			result += HEX(data[i]);
		}
		return result;
	}
	const String
		encodehex(const Byte data[], size_t len)
	{ 
		String result;
		result.reserve(len * 2);
		for (size_t i = 0; i < len; ++i)
		{
			result += HEX(data[i]);
		}
		return result;
	}


	extern String
		basenumber(
			unsigned int v0,
			unsigned int base,
			const char szDigital[]
		)
	{
		assert(base > 0);
		assert(strlen(szDigital) >= base);

		String result;
		while(v0 > 0)
		{
			result += szDigital[v0 % base];
			v0 /= base;
		}
		std::reverse(result.begin(), result.end());

		return result;
	}

	extern String
		hex(unsigned int v0)
	{
		const static char szDigital[] = "0123456789ABCDEF";
		return basenumber(v0, 16, szDigital);
	}

	extern String
		oct(unsigned int v0)
	{
		const static char szDigital[] = "01234567";
		return basenumber(v0, 8, szDigital);
	}

	extern String
		bin(unsigned int v0)
	{
		const static char szDigital[] = "01";
		return basenumber(v0, 2, szDigital);
	}

	extern String
		addslash(const char* pStr, size_t n)
	{
		String result;
		for (const char* const pEnd = pStr + n; pStr != pEnd; ++pStr)
		{
            switch(*pStr)
			{
			case '\n':		// New line
				result += "\\n";
				break;
			case '\r':		// Return
				result += "\\r";
				break;
			case '\0':		// Nil
				result += "\\0";
				break;
			case '\032':	// Ctrl+Z; Sub
				result += "\\Z";
				break;

			case '\'':
			case '\"':
			case '\\':
				result += '\\';
				result += *pStr;
				break;

			default:
				result += *pStr;
				break;
			}
		}

		return result;
	}

	
	namespace 
	{
		static const Char defaultBase64Code[] 
			= "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/=";

		typedef std::map<Char, boost::uint8_t> Base64CodeMap;

		Base64CodeMap createCodeMap(const Char szCode[])
		{
			Base64CodeMap result;
			for (boost::uint8_t i = 0; i < 64; ++i)
			{
				result[szCode[i]] = i;
			}
			assert (result.size() == 64);
			return result;
		}


		void encodeBase64(const boost::uint8_t code[3], const Char szCodeBase[], String& result)
		{
			result += szCodeBase[code[0] & 0x3f];
			result += szCodeBase[(code[0] >> 6) | ((code[1] & 0xf) << 2)];
			result += szCodeBase[(code[1] >> 4) | (code[2] & 0x3) << 4];
			result += szCodeBase[code[2] >> 2];
		}

		void decodeBase64(const boost::uint8_t code[4], Char result[3])
		{
			assert (code[0] < 64);
			assert (code[1] < 64);
			assert (code[2] < 64);
			assert (code[3] < 64);

			result[0] = code[0] | code[1] << 6;
			result[1] = code[1] >> 2 | code[2] << 4;
			result[2] = code[2] >> 4 | code[3] << 2;
		}

		char decodeBase64(const Char code[4], const Base64CodeMap& codeMap, Char result[3])
		{
			boost::uint8_t number[4] = {0, 0, 0, 0};
			char i = 0;
			for (; i < 4; ++i)
			{
				Base64CodeMap::const_iterator it = codeMap.find(code[i]);
				if (it != codeMap.end())
				{
					number[i] = it->second;
				}
				else
				{
					break;
				}
			}
			decodeBase64(number, result);
			return i;
		}
	}	// namespace 

	String encodeBase64(const Char sData[], boost::uint32_t length, const Char szCodeBase[])
	{
		String result;
		boost::uint8_t buffer[3];

		const boost::uint32_t len3 = (length / 3) * 3;
		for (const Char* p = sData; p != sData+len3; p += 3)
		{
			buffer[0] = boost::uint8_t(*p); 
			buffer[1] = boost::uint8_t(*(p+1)); 
			buffer[2] = boost::uint8_t(*(p+2)); 
			encodeBase64(buffer, szCodeBase, result);
		}
		buffer[0] = buffer[1] = buffer[2] = 0;
		for (boost::uint32_t i = len3; i < length; ++i)
		{
			buffer[i - len3] = sData[i];
		}
		String tail;
		encodeBase64(buffer, szCodeBase, tail);
		if (length % 3 == 1)
		{
			tail[2] = tail[3] = szCodeBase[64];
			result += tail;
		}
		else if (length % 3 == 2)
		{
			tail[3] = szCodeBase[64];
			result += tail;
		}
		
		return result;
	}

	String encodeBase64(const Char data[], boost::uint32_t length)
	{
		return encodeBase64(data, length, defaultBase64Code);
	}

	const Char* decodeBase64(const Char data[], const Char szCodeBase[], String& result)
	{
		result.clear();

		Base64CodeMap code2number = createCodeMap(szCodeBase);
		const Char* next = data;
		
		char n;
		do
		{
			Char v[3];
			n = decodeBase64(next, code2number, v);
			next += n;

			switch(n)
			{
			case 1:
			case 2:
				result += v[0];
				break;

			case 3:
				result += v[0];
				result += v[1];
				break;

			case 4:
				result += v[0];
				result += v[1];
				result += v[2];
				break;

			default:
				break;
			}

		} while (n == 4);


		return next;
	}

	const Char* decodeBase64(const Char data[], String& result)
	{
		return decodeBase64(data, defaultBase64Code, result);
	}

	String decodeBase64(const Char data[], const Char code[])
	{
		String result;
		decodeBase64(data, code, result);
		return result;
	}

	String decodeBase64(const Char data[])
	{
		return decodeBase64(data, defaultBase64Code);
	}

	String
		HEX(boost::uint8_t v0)
	{
		Char buf[2] = {v0 / 16, v0 % 16};
		for (Char* p = buf; p != buf+2; ++p)
		{
			if (*p < 10)
			{
				*p  += '0';
			}
			else
			{
				*p += 'A' - 10;
			}
		}
		
		return String(buf, 2);
	}
	String
		encodeHEX(const Char data[], size_t size)
	{
		String result;
		for (size_t i = 0; i < size; ++i)
		{
			result += HEX(Byte(data[i]));
		}
		return result;
	}
	String
		HEX(boost::uint16_t v0);
	String
		HEX(boost::uint32_t v0);
	String
		HEX(boost::uint64_t v0);

	String
		hex(boost::uint8_t v0)
	{
		Char buf[2] = {v0 / 16, v0 % 16};
		for (Char* p = buf; p != buf+2; ++p)
		{
			if (*p < 10)
			{
				*p  += '0';
			}
			else
			{
				*p += 'a' - 10;
			}
		}
		
		return String(buf, 2);
	}

}
