#include <boost/filesystem.hpp>
#include <kuai/typedef.h>
#include <kuai/tools/strtool.h>


#ifndef _KUAI_TOOLS_FILE_NAME_H_
#define _KUAI_TOOLS_FILE_NAME_H_

namespace kuai {

	class FileName
		: public boost::filesystem::path
	{
	public:
		FileName()
		{ }
		explicit FileName(const Char szName[]) 
			: boost::filesystem::path(szName)
		{ }

		explicit FileName(const String& s) 
			: boost::filesystem::path(s)
		{ }

		FileName(const boost::filesystem::path& v)
			: boost::filesystem::path(v)
		{ }

	public:
		bool is_file() const {
			return boost::filesystem::is_regular_file(*this);
		}

		bool is_dir() const {
			return boost::filesystem::is_directory(*this);
		}

		bool empty() const {
			return string().empty();
		}

		const String extname() const {
			return boost::filesystem::extension(*this);		
		}
	};

	template<>
	inline	String str<FileName>(const FileName& v0)
	{
		return v0.string();
	};

}

#endif
