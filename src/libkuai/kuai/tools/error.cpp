#include <boost/system/system_error.hpp>


#ifdef WIN32
	#include <windows.h>
	#include <kuai/tools/error.h>
	namespace kuai {
		void check_last_error() {
			int err = ::GetLastError();
			if (err != 0) {
				boost::system::error_code code(err, boost::system::system_category());
				throw boost::system::system_error(code);
			}
		};
	}

#else 
	#include <error.h>
	#include <kuai/tools/error.h>
	namespace kuai {
		void check_last_error() {
			if (errno != 0) {
				boost::system::error_code code(errno, boost::system::system_category());
				throw boost::system::system_error(code);
			}
		};
	}

#endif 

#include <fstream>

namespace kuai {
	
	extern const Char SZ_UNKNOWN_ERROR[] = "Unknown program error";

	void tryOpen(std::ifstream& stream, const FileName& file) {
		if (file.empty()) {
			throw error("File name is empty.");
		}
		else if (!file.is_file()) {
			throw error("File %1% does not existed.", file);
		}
		else {
			stream.clear();
			stream.open(file.c_str());
			if (!stream) {
				check_last_error();
			}
		}
	}
	void tryOpen(std::ofstream& stream, const FileName& file) {
		if (file.empty()) {
			throw error("File name is empty.");
		}
		else {
			stream.clear();
			stream.open(file.c_str());
			if (!stream) {
				check_last_error();
			}
		}
	}

}

