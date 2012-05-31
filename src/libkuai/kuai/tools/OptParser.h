#include <ostream>
#include <boost/tuple/tuple.hpp>
#include <kuai/typedef.h>

#ifndef _KUAI_TOOLS_OPT_PARSTER_H_
#define _KUAI_TOOLS_OPT_PARSTER_H_

namespace kuai {


	class OptParserHandle {
	public:
		virtual ~OptParserHandle();
		virtual Index operator()(const String& tag, char* p[]) const = 0;

		virtual void help(std::ostream& stream) const = 0;
	};

	template<typename T> 
		class OptParserHandleSetVar
			: public OptParserHandle
	{
	public:
		explicit OptParserHandleSetVar(T* target, const Char szHelp[])
			: _target(target), _help(szHelp)
		{ }

	public:
		virtual Index operator()(const String& tag, char* p[]) const {
			*_target = lexical_cast<T>(String(*p));
			return 1;
		}
		virtual void help(std::ostream& stream) const {
			stream << _help << std::endl;
		}

	private:
		T* _target;
		String _help;
	};

	template<typename TData, typename TContainer> 
		class OptParserHandlePushVar
			: public OptParserHandle
	{
	public:
		explicit OptParserHandlePushVar(TContainer* target, const Char szHelp[])
			: _target(target), _help(szHelp)
		{ }

	public:
		virtual Index operator()(const String& tag, char* p[]) const {
			_target->push_back(lexical_cast<TData>(*p));
			return 1;
		}
		virtual void help(std::ostream& stream) const {
			stream << _help << std::endl;
		}

	private:
		TContainer* _target;
		String _help;
	};

	template<typename T> 
		class OptParserHandleSetEnum
			: public OptParserHandle
	{
	public:
		explicit OptParserHandleSetEnum(T* target, const Char szHelp[])
			: _target(target), _help(szHelp)
		{ }

	public:
		virtual Index operator()(const String& tag, char* p[]) const {
			String v(*p);
			typename KuaiMap<String, T>::const_iterator it = _mapping.find(v);
			if (it == _mapping.end()) {
				throw error("Unknonw Option %1%", v);
			}
			else {
				*_target = it->second;
			}
			return 1;
		}
		virtual void help(std::ostream& stream) const {
			stream << _help << std::endl;
		}

	protected:
		T* _target;
		KuaiMap<String, T> _mapping;
		String _help;
	};


	class OptParser {
	public:
		OptParser();
		virtual ~OptParser();
		void add_handle(const String& short_tag, const String& long_tag, boost::shared_ptr<OptParserHandle> handle);

		template<typename T>
			void add_handle_set_var(const String& short_tag, const String& long_tag, T* target, const Char szHelp[]="")
		{
			 boost::shared_ptr<OptParserHandle> handle(new OptParserHandleSetVar<T>(target, szHelp));
			 add_handle(short_tag, long_tag, handle);
		}

		void add_handle_reset_var(const String& short_tag, const String& long_tag, bool* target, const Char help[]=NULL);	
		
		template<typename T>
		void add_handle_push_var(const String& short_tag, const String& long_tag, Array<T>* target, const Char szHelp[]=NULL) {
			boost::shared_ptr<OptParserHandle> handle(new OptParserHandlePushVar<T, Array<T> >(target, szHelp));
			add_handle(short_tag, long_tag, handle);
		};	
		
		bool parse(int argc, char* argv[], std::ostream& output = std::cout);

		const StringArray& arguments() const {
			return _args;
		}
		const String& arguments(size_t i) const {
			return _args[i];
		}

		virtual void on_parse_success() = 0;
		virtual void check_parse() const = 0;
		virtual void basic_usage(std::ostream& stream) const = 0;

		virtual void on_help(std::ostream& stream) const;
		virtual void on_help(std::ostream& stream, const String& par) const;

	protected:
		String        _command;
		Array<String> _args;
		Array<boost::tuple<char, String, SharedPtr<OptParserHandle> > > _handles;

	private:
		SharedPtr<OptParserHandle> _find(const String& tag) const;
		SharedPtr<OptParserHandle> _find(char tag) const;
	};


	template<>
		void OptParser::add_handle_set_var<bool>(const String& short_tag, const String& long_tag, bool* target, const Char szHelp[]);

}

#endif
