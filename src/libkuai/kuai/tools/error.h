#include <stdexcept>
#include <kuai/tools/strtool.h>
#include <kuai/tools/FileName.h>

#ifndef _KUAI_ERROR_H_
#define _KUAI_ERROR_H_

namespace kuai {

	inline void dumpError(std::exception& error, std::ostream& stream=std::cerr) {
		stream << "\n#Error: " << error.what() << "\n" << std::endl;
		stream.flush();
	}

	inline void dumpWarning(const char* error, std::ostream& stream=std::cout) {
		stream << "\n#Warning: " << error << "\n" << std::endl;
		stream.flush();
	}

	inline void dumpWarning(std::exception& error, std::ostream& stream=std::cout) {
		dumpWarning(error.what(), stream);
	}

	extern const Char SZ_UNKNOWN_ERROR[]; // = "Unknown program error";

	inline std::runtime_error error() {
		return std::runtime_error(SZ_UNKNOWN_ERROR);
	}

	inline std::runtime_error error(const String& error) {
		return std::runtime_error(error);
	}

	template<typename T1>
	inline std::runtime_error error(const String& error, const T1& v1) {
		StringTemplate tt(error);
		tt["%1%"] = str(v1);
		return std::runtime_error(tt.str());
	}


	template<typename T1, typename T2>
	inline std::runtime_error error(const String& error, const T1& v1, const T2& v2) {
		StringTemplate tt(error);
		tt["%1%"] = str(v1);
		tt["%2%"] = str(v2);
		return std::runtime_error(tt.str());
	}

	template<typename T1, typename T2, typename T3>
	inline std::runtime_error error(const String& error, const T1& v1, const T2& v2, const T3& v3) {
		StringTemplate tt(error);
		tt["%1%"] = str(v1);
		tt["%2%"] = str(v2);
		tt["%3%"] = str(v3);
		return std::runtime_error(tt.str());
	}

	void check_last_error();

	void tryOpen(std::ifstream& stream, const FileName& file);
	void tryOpen(std::ofstream& stream, const FileName& file);


	class NotImplementedError
		: public std::runtime_error
	{ 
	public:
		NotImplementedError(const String& cls, const String& methods)
			: std::runtime_error(cls + ":" + methods + " is not implemented yet.")
		{ }
	};

	#define NOT_IMPLEMENTED(methods) NotImplementedError(typeid(*this).name(), #methods)

	class UndefinedMember
		: public std::runtime_error
	{ 
	public:
		UndefinedMember(const String& cls, const String& member)
			: std::runtime_error(cls + "::" + member + " is not defined.")
		{ }
	};
	#define UNDEFINED_MEMBER(member) UndefinedMember(typeid(*this).name(), member)

	inline std::ostream& operator<<(std::ostream& stream, const std::exception& err) {
		return stream << err.what();
	}

}

#endif
