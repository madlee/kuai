#include <kuai/tools/error.h>
#include <kuai/tools/strtool.h>
#include <kuai/tools/OptParser.h>


namespace kuai {

	OptParserHandle::~OptParserHandle()
	{ }

	namespace {
		
		class OptParserHandleSetBoolean
			: public OptParserHandle
		{
		public:
			explicit OptParserHandleSetBoolean(bool* target, bool value, const Char szHelp[])
				: _target(target), _value(value), _help(szHelp)
			{ }

		public:
			virtual Index operator()(const String& tag, char* p[]) const {
				*_target = _value;
				return 0;
			}

			virtual void help(std::ostream& stream) const {
				stream << _help << std::endl;
			}

		private:
			bool* _target;
			bool _value;
			String _help;
		};

		class OptParserHandleHelp
			: public OptParserHandle
		{
		public:
			explicit OptParserHandleHelp(OptParser* parser)
				: _parser(parser)
			{ }

			virtual Index operator()(const String& tag, char* p[]) const {
				if (*p == NULL) {
					_parser->on_help(std::cout);
					return 0;
				}
				else {
					Index n = 0;
					for (n = 0; p[n] != NULL; ++n) {
						_parser->on_help(std::cout, p[n]);
					}
					return n;
				}
			}

			virtual void help(std::ostream& stream) const {
				stream << "Show help message" << std::endl;
			}

		private:
			OptParser* _parser;
		};
	}

	OptParser::OptParser() {
		add_handle("-h", "--help", SharedPtr<OptParserHandle>(new OptParserHandleHelp(this)));
	}
	OptParser::~OptParser() 
	{ }

	template<>
		void OptParser::add_handle_set_var<bool>(const String& short_tag, const String& long_tag, bool* target, const Char szHelp[]) 
	{
		add_handle(short_tag, long_tag, boost::shared_ptr<OptParserHandle>(new OptParserHandleSetBoolean(target, true, szHelp)));
	}

	
	void OptParser::add_handle_reset_var(const String& short_tag, const String& long_tag, bool* target, const Char szHelp[]) {
		add_handle(short_tag, long_tag, boost::shared_ptr<OptParserHandle>(new OptParserHandleSetBoolean(target, false, szHelp)));
	}

	void OptParser::add_handle(const String& short_tag, const String& long_tag, boost::shared_ptr<OptParserHandle> handle) {
		char tag1 = '\0';
		if (!short_tag.empty()) {
			assert (short_tag.size() == 2 && short_tag[0] == '-');
			tag1 = short_tag[1];
		}
		String tag2;
		if (!long_tag.empty()) {
			assert (long_tag.size() > 2 && long_tag[0] == '-' && long_tag[1] == '-');
			tag2 = long_tag.substr(2);
		}
		_handles.push_back(boost::make_tuple(tag1, tag2, handle));
	}

	SharedPtr<OptParserHandle> OptParser::_find(const String& tag) const
	{
		for (Array<boost::tuple<char, String, SharedPtr<OptParserHandle> > >::const_iterator 
			i = _handles.begin(); i != _handles.end(); ++i)
		{
			if (i->get<1>() == tag) {
				return i->get<2>();
			}
		}
		throw error("Option %1% is unknown.", "--"+tag);
	}
	SharedPtr<OptParserHandle> OptParser::_find(char tag) const
	{
		for (Array<boost::tuple<char, String, SharedPtr<OptParserHandle> > >::const_iterator 
			i = _handles.begin(); i != _handles.end(); ++i)
		{
			if (i->get<0>() == tag) {
				return i->get<2>();
			}
		}
		throw error("Option %1% is unknown.", "-"+tag);
	}

	bool OptParser::parse(int argc, char* argv[], std::ostream& output) {
		_args.reserve(argc);
		_command = argv[0];
		try {
			for (Index i = 1; i < argc; ++i) {
				String v(argv[i]);
				if (v.size() > 2 && v[0] == '-' && v[1] == '-') {
					StringPair vv = split_pair(v, "=");
					String tag = vv.first.substr(2);
					SharedPtr<OptParserHandle> pfunc = _find(tag);
					char* p[2] = {NULL, NULL};
					if (!vv.second.empty()) {
						p[0] = &vv.second[0];
					}
					(*pfunc)(vv.first, p);
				}
				else if (v.size() > 1 && v[0] == '-') {
					i += 1;
					for (size_t j = 1; j < v.size(); ++j) {
						char tag = v[j];
						SharedPtr<OptParserHandle> pfunc = _find(tag);
						i += (*pfunc)(v, argv+i);
					}
					i -= 1;
				}
				else {
					_args.push_back(argv[i]);
				}
			}
			check_parse();
			on_parse_success();
			return true;
		}
		catch (std::exception& error) {
			dumpError(error);
			on_help(std::cerr);
			return false;
		}
	}

	void OptParser::on_help(std::ostream& stream) const {
		basic_usage(stream);
		stream << "[Options]\n";
		for (Array<boost::tuple<char, String, SharedPtr<OptParserHandle> > >::const_iterator
			it = _handles.begin(); it != _handles.end(); ++it)
		{
			stream << "    ";
			char c = it->get<0>();
			if (c != char(0)) {
				stream << "-" << c << ", ";
			}
			String s = it->get<1>();
			if (!s.empty()) {
				stream << "--" << s << ",";
			}
			stream << std::endl;
			it->get<2>()->help(stream);
		}
	}

	void OptParser::on_help(std::ostream& stream, const String& tag) const {
		try {
			SharedPtr<OptParserHandle> handle;
			if (tag.size() == 1) {
				handle = _find(tag[0]);
			}
			else {
				handle = _find(tag);
			}
			handle->help(stream);
		}
		catch (std::exception& error) {
			dumpError(error, stream);
		}
	}

}

