#include <kuai/tools/time.h>
#include <kuai/tools/Logger.h>

namespace kuai {

	const String Logger::_names[END_OF_LOG] = {
		"", 
		"[CRITICAL] ", 
		"[ ERROR! ] ", 
		"[WARNING!] ", 
		"[  INFO  ] ",
		"[ DETAIL ] ", 
		"[ DEBUG! ] "
	};

	
	Logger::Logger()
		: _level(LOG_INFO), _stream(&std::clog) 
	{ }
	
	void Logger::set_stream(const FileName& filename, bool append) {
		if (append) {
			_logfile = SharedPtr<std::ofstream>(new std::ofstream(filename.c_str(), std::ios::app));
		}
		else {
			_logfile = SharedPtr<std::ofstream>(new std::ofstream(filename.c_str()));
		}
		_stream = _logfile.get();
	};
	void Logger::set_stream(std::ostream& stream) {
		_stream = &stream;
		_logfile = SharedPtr<std::ofstream>();
	};

	Logger& Logger::get_instance() {
		static Logger instance;
		return instance;
	};

	
	bool Logger::lock(LogLevel level) {
		if (std::ostream* stream = get_stream(_level)) {
			_mutex.lock();
			*stream << now() << _names[_level];
			return true;
		}
		else {
			return false;
		}
	}

	bool Logger::free(LogLevel level) {
		if (std::ostream* stream = get_stream(_level)) {
			*stream << std::endl;
			_mutex.unlock();
			return true;
		}
		else {
			return false;
		}
	}

	std::ostream* Logger::get_stream(LogLevel level) {
		if (level <= _level) {
			return _stream;
		}	
		else {
			return NULL;
		}
	}

	
	std::ostream* LoggerProxy::get_stream() {
		return Logger::get_instance().get_stream(_level);
	}

	LoggerProxy& LoggerProxy::lock(const String& head) {
		if (Logger::get_instance().lock(_level)) {
			*this << head;
		};
		return *this;
	}
		
	LoggerProxy& LoggerProxy::free() {
		Logger::get_instance().free(_level);
		return *this;
	};
	
	void LoggerProxy::dump(const String& line) {
		this->lock() << line << std::endl;
	};

	
	LoggerProxy log_critical(LOG_CRITICAL);
	LoggerProxy log_error(LOG_ERROR);
	LoggerProxy log_warning(LOG_WARNING);
	LoggerProxy log_info(LOG_INFO);
	LoggerProxy log_detail(LOG_DETAIL);


#ifdef NDEBUG
	DumbLogger log_debug;
#else
	LoggerProxy log_debug(LOG_DEBUG);
#endif


}
