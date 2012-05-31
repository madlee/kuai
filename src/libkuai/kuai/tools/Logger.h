#include <cstring>
#include <cassert>

#include <ostream>

#include <kuai/typedef.h>
#include <kuai/tools/FileName.h>
#include <kuai/tools/thread.h>


#ifndef _KUAI_LOGGER_H_
#define _KUAI_LOGGER_H_

namespace kuai {

	enum LogLevel {
		LOG_NONE, LOG_CRITICAL, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DETAIL, LOG_DEBUG, END_OF_LOG
	};

	class Logger {
	public:
		Logger();

		virtual ~Logger()
		{ }
		void set_stream(const FileName& filename, bool append=true);
		void set_stream(std::ostream& stream);

		bool lock(LogLevel level);
		bool free(LogLevel level);

		std::ostream* get_stream(LogLevel level);

		LogLevel get_current_level() {
			return _level;
		}
		void set_current_level(LogLevel level) {
			_level = level;
		}

	public:
		static Logger& get_instance();

	private:
		LogLevel _level;
		Mutex _mutex;

		std::ostream* _stream;
		SharedPtr<std::ofstream> _logfile;

		static const String _names[END_OF_LOG];
	};

	class LoggerProxy {
	public:
		LoggerProxy(LogLevel level)
			: _level(level)
		{ }

		std::ostream* get_stream();

		LoggerProxy& lock(const String& head = "");
		LoggerProxy& free();
		void dump(const String& line);


	private:
		LogLevel _level;
	};

	template<typename T>
		inline LoggerProxy& operator<<(LoggerProxy& log, const T& v) 
	{
		if (std::ostream* stream = log.get_stream()) {
			*stream << v;
		}
		return log;
	}


	// typedef std::ostream& (*StreamFunction)(std::ostream&) ;

	inline LoggerProxy& operator<< (LoggerProxy& log, std::ostream& (*v)(std::ostream&)) 
	{
		if (std::ostream* stream = log.get_stream()) {
			if (v == std::endl) {
				log.free();
			}
			else {
				*stream << v;
			}
		}
		return log;
	}

	class DumbLogger{
	public:
		std::ostream* get_stream() {
			return NULL;
		}
		DumbLogger& lock(const String& head = "") 
		{
			return *this;
		}
		DumbLogger& free()
		{
			return *this;
		}
		void dump(const String& line) 
		{ }
	};

	template<typename T>
		inline DumbLogger& operator<<(DumbLogger& log, const T& v) 
	{
		return log;
	}

	extern LoggerProxy log_critical;
	extern LoggerProxy log_error;
	extern LoggerProxy log_warning;
	extern LoggerProxy log_info;
	extern LoggerProxy log_detail;


#ifdef NDEBUG
	extern DumbLogger log_debug;
#else
	extern LoggerProxy log_debug;
#endif

}

#endif
