#ifndef SWIG
	#include <boost/optional.hpp>
	#include <kuai/typedef.h>
	#include <kuai/tools/sqlite3db.h>
	#include <kuai/tools/ProviderConsumerBuffer.h>
	#include <kuai/mols/Molecule.h>
#endif

#ifndef _KUAICHOP_H_
#define _KUAICHOP_H_

class KuaiChopJob {
public:
	explicit KuaiChopJob(const char database[], size_t nThreads=1, size_t buffer_size = 1024);
	virtual ~KuaiChopJob();

public:
	void feed(const char smiles[]);
	void close();

	bool cut_single_chain() const  {
		return _cut_single_chain;
	}
	void cut_single_chain(bool v) {
		_cut_single_chain = v;
	}

#ifndef SWIG

public:
	kuai::String pop() {
		return _buffer.pop();
	};
	
	void log_error(const std::exception& error);
	void log_error(const kuai::String& error);
	void log_warning(const kuai::String& warning);
	void log_info(const kuai::String& info);

	template<typename T1>
		void log_info(const kuai::String& info, const T1& v1)
	{
		StringTemplate t(info);
		t["%1%"] = str(v1);
		log_info(t.str());
	}

	template<typename T1>
		void log_warning(const kuai::String& info, const T1& v1)
	{
		StringTemplate t(info);
		t["%1%"] = str(v1);
		log_warning(t.str());
	}


	kuai::Sqlite3RowID add_molecule(kuai::MoleculePtr& mol);

private:
	bool											_cut_single_chain;
	bool											_cut_fused_rings;
	bool											_cut_spiro_rings;
	boost::optional<kuai::Index>					_cut_big_ring;
	kuai::Index										_min_fragment_size, _max_fragment_size;
	kuai::Index										_max_cuts;

	boost::optional<int>							_isotope_base;
	int												_isotope_shift;
	bool											_fix_cutpoint_hydrogen;

	kuai::Sqlite3DatabasePointer					_database;
	kuai::ProviderConsumerBuffer<kuai::String>		_buffer;
	kuai::ThreadGroup								_groups;

	kuai::Mutex										_access_db;
	kuai::Mutex										_access_log;

#endif
};

#endif
