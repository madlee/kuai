#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <kuai/typedef.h>
#include <kuai/tools/sqlite3.h>
#include <kuai/tools/strtool.h>
#include <kuai/tools/error.h>


#ifndef _KUAI_SQLITE3_HPP_
#define _KUAI_SQLITE3_HPP_

namespace kuai { 

	typedef sqlite_int64 Sqlite3RowID;

	class Sqlite3Statement 
		: boost::noncopyable
	{
	private:
		explicit Sqlite3Statement(sqlite3_stmt* v0) 
			: stmt(v0)
		{ }
	public:
		~Sqlite3Statement() {
			sqlite3_finalize(stmt);
		}

		void reset() {
			sqlite3_reset(stmt);			
		}

		void bind(int i, int value) {
			check_bind(i, sqlite3_bind_int(stmt, i, value));
		}
		void bind(int i, unsigned int value) {
			check_bind(i, sqlite3_bind_int64(stmt, i, value));
		}
		void bind(int i, sqlite3_int64 value) {
			check_bind(i, sqlite3_bind_int64(stmt, i, value));
		}
		void bind(int i, RealNumber value) {
			check_bind(i, sqlite3_bind_double(stmt, i, value));
		}
		void bind(int i, const String& value) {
			check_bind(i, sqlite3_bind_text(stmt, i, value.c_str(), value.size(), SQLITE_TRANSIENT));
		}
		void bind(int i, const char* blob, size_t n) {
			check_bind(i, sqlite3_bind_blob(stmt, i, blob, n, SQLITE_TRANSIENT));
		}

		void bind_null(int i) {
			check_bind(i, sqlite3_bind_null(stmt, i));
		}

		int step() {
			int result = sqlite3_step(stmt);
			switch (result) {
			case SQLITE_OK:
			case SQLITE_DONE:
			case SQLITE_ROW:
				return result;

			default:
				{
					sqlite3* handle = sqlite3_db_handle(stmt);
					const char* error = sqlite3_errmsg(handle);
					throw std::runtime_error(error);
				}
			}
		}

		boost::int64_t getInt64(int icol) {
			return sqlite3_column_int64(stmt, icol);
		}

		const char* c_str() const {
			return sqlite3_sql(stmt);
		}

		size_t countColumn() const {
			return sqlite3_column_count(stmt);
		}

		const Char* columnName(size_t i) const {
			return sqlite3_column_name(stmt, i);
		}

		const Char* column(size_t i) const {
			return reinterpret_cast<const Char*>(sqlite3_column_text(stmt, i));
		}

		int columnSize(size_t i) const {
			return sqlite3_column_bytes(stmt, i);
		}
		std::pair<int, const char*> getBLOB(int icol) {
			int n = columnSize(icol);
			const char* p = (const char*)(sqlite3_column_blob(stmt, icol));
			return std::make_pair(n, p);
		}

	private:
		void check_bind(int i, int result);
		sqlite3_stmt* stmt;
		friend class Sqlite3Database;
	};

	template<>
	inline String str<Sqlite3Statement>(const Sqlite3Statement& v0)
	{
		return v0.c_str();
	}

	class ErrorCanNotFindRow
		: public std::runtime_error
	{
	public:
		explicit ErrorCanNotFindRow(const String& row) 
			: std::runtime_error(error("Can not find row for sql %1%", row))
		{ }
	};

	class Sqlite3Database 
		: boost::noncopyable
	{
	public:
		typedef int (*ExecCallBack)(void*,int,char**,char**);  
	private:
		Sqlite3Database(sqlite3* p0);

	public:
		~Sqlite3Database();

	public:
		static boost::shared_ptr<Sqlite3Database> create(const FileName& szName, int mode = SQLITE_OPEN_READWRITE  | SQLITE_OPEN_CREATE);

		int execute(const Char* szsql);
		int execute(const Char* szsql, ExecCallBack func);
		int execute(const Char* szsql, ExecCallBack func, void* par0);

		void begin();
		void commit();

		boost::shared_ptr<Sqlite3Statement> prepare(const String& sql);

		static void check(int result);
		std::runtime_error last_error();

		boost::int64_t getSingleInt(Sqlite3Statement& sql);
		boost::int64_t getSingleInt(const Char szSql[]);
		boost::int64_t getSingleInt(const String& sql);
		void getIntArray(const String& sql, int col, std::vector<int>& result);
		void getIntArray(Sqlite3Statement& stat, int col, std::vector<int>& result);
		void getStringArray(const String& sql, int col, StringArray& result);
		void getStringArray(Sqlite3Statement& sql, int col, StringArray& result);

		boost::int64_t insert(Sqlite3Statement& sql);
		boost::int64_t insert(const String& sql);
		boost::int64_t insert(const Char szsql[]);

		void getTable(const Char szSql[], std::vector<StringArray>& result);
	private:
		sqlite3* pdb;
	};

	typedef boost::shared_ptr<Sqlite3Database> Sqlite3DatabasePointer;
	typedef boost::shared_ptr<Sqlite3Statement> Sqlite3StatementPointer;
}


#endif
