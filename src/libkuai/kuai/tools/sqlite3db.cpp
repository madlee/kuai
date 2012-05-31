#include <kuai/tools/error.h>
#include <kuai/tools/sqlite3db.h>

namespace kuai { 

	Sqlite3Database::Sqlite3Database(sqlite3* p0) 
		: pdb(p0)
	{ }

	Sqlite3Database::~Sqlite3Database() {
		sqlite3_close(pdb);
	}

	boost::shared_ptr<Sqlite3Database> Sqlite3Database::create(const FileName& szName, int mode) 
	{
		String name = encode_utf8(szName.string());
		sqlite3* pdb = NULL;
		if (sqlite3_open_v2(name.c_str(), &pdb, mode, NULL) == SQLITE_OK) {
			return boost::shared_ptr<Sqlite3Database>(new Sqlite3Database(pdb));
		}
		else if (pdb == NULL) {
			throw error("Failed To Open sqlite3 database %1%", szName);
		}
		else {
			const char* sze = sqlite3_errmsg(pdb);
			std::runtime_error error(sze);
			sqlite3_close(pdb);
			throw error;
		};
	};

	int Sqlite3Database::execute(const Char* szsql) {
		char* error;
		int result = sqlite3_exec(pdb, szsql, NULL, NULL, &error);
		if (result == SQLITE_OK) {
			return SQLITE_OK;
		}
		else if (result == SQLITE_ABORT) {
			sqlite3_free(error);
			return result;
		}
		else {
			std::runtime_error err(error);
			sqlite3_free(error);
			throw err;
		}
	};

	int Sqlite3Database::execute(const Char* szsql, ExecCallBack func) {
		char* error;
		int result = sqlite3_exec(pdb, szsql, func, this, &error);
		if (result == SQLITE_OK) {
			return SQLITE_OK;
		}
		else if (result == SQLITE_ABORT) {
			sqlite3_free(error);
			return result;
		}
		else {
			std::runtime_error err(error);
			sqlite3_free(error);
			throw err;
		}
	};

	int Sqlite3Database::execute(const Char* szsql, ExecCallBack func, void* par0) {
		char* error;
		int result = sqlite3_exec(pdb, szsql, func, par0, &error);
		if (result == SQLITE_OK) {
			return SQLITE_OK;
		}
		else if (result == SQLITE_ABORT) {
			sqlite3_free(error);
			return result;
		}
		else {
			std::runtime_error err(error);
			sqlite3_free(error);
			throw err;
		}
	};

	void Sqlite3Database::begin() {
		execute("BEGIN TRANSACTION");
	}
	void Sqlite3Database::commit() {
		execute("COMMIT TRANSACTION");
	}

	void Sqlite3Database::check(int result) {
		// TODO:
		if (result != SQLITE_OK) {
			throw error("SQLITE_ERROR %1%", result);
		}
	}

	boost::shared_ptr<Sqlite3Statement> Sqlite3Database::prepare(const String& sql) {
		sqlite3_stmt* stmt;
		const char* pzTail;
		int result = sqlite3_prepare_v2(pdb, sql.c_str(), sql.size(), &stmt, &pzTail);
		
		if (result == SQLITE_OK) {
			return boost::shared_ptr<Sqlite3Statement>(new Sqlite3Statement(stmt));
		}
		else {
			const char* what = sqlite3_errmsg(pdb);
			throw std::runtime_error(what);	
		}
	}

	void Sqlite3Statement::check_bind(int i, int result) {
		if (result == SQLITE_OK) {
		}
		else if (result == SQLITE_RANGE) {
			throw error("There is no %1% parameter in SQL.");
		}
		else if (result == SQLITE_NOMEM) {
			throw std::bad_alloc();
		}
		else {
			sqlite3* handle = sqlite3_db_handle(stmt);
			const char* error = sqlite3_errmsg(handle);
			throw std::runtime_error(error);
		}
	}

	boost::int64_t Sqlite3Database::getSingleInt(const String& sql) {
		return getSingleInt(sql.c_str());
	}

	boost::int64_t Sqlite3Database::getSingleInt(const Char szSql[]) {
		boost::shared_ptr<Sqlite3Statement> p = prepare(szSql);
		return getSingleInt(*p);
	}

	boost::int64_t Sqlite3Database::getSingleInt(Sqlite3Statement& sql) {
		sql.reset();
		int code = sql.step();
		switch (code) {
		case SQLITE_ROW:
			if (sqlite3_column_count(sql.stmt) == 0) {
				throw error("The sql %1% return nothing.", sql);
			}
			if (sqlite3_column_count(sql.stmt) > 1) {
				throw error("The sql %1% return more than one column.", sql);
			}
			if (sqlite3_column_type(sql.stmt, 0) != SQLITE_INTEGER) {
				throw error("The sql %1% dose not return an integer.", sql);
			}
			else {
				boost::int64_t result = sql.getInt64(0);
				code = sql.step();
				if (code == SQLITE_DONE) {
					return result;
				}
				else {
					throw error("The sql %1% return more than one line.", sql);
				}
				return result;
			}
			break;

		case SQLITE_DONE:
			throw ErrorCanNotFindRow(str(sql));
			break;

		default:
			{
				throw last_error();
			}
		}
	}

	boost::int64_t Sqlite3Database::insert(Sqlite3Statement& sql) {
		sql.step();
		return sqlite3_last_insert_rowid(pdb);
	}

	boost::int64_t Sqlite3Database::insert(const String& sql) {
		return insert(sql.c_str());
	}
	
	boost::int64_t Sqlite3Database::insert(const Char szsql[]) {
		return insert(*prepare(szsql));
	}

	void Sqlite3Database::getTable(const Char szSql[], std::vector<StringArray>& result) {
		char** presult;
		int row, col;
		char* pzErrmsg = NULL;
		int code = sqlite3_get_table(pdb, szSql, &presult, &row, &col, &pzErrmsg);

		if (code == SQLITE_OK) {
			result.resize(row+1);
			for (int i = 0; i < row+1; ++i) {
				result[i].resize(col);
				for (int j = 0; j < col; ++j) {
					result[i][j] = presult[i*col+j];		
				}
			}
			sqlite3_free_table(presult);
		}
		else {
			std::runtime_error result = error(pzErrmsg);
			sqlite3_free(pzErrmsg);
			throw result;
		}
	}

	void Sqlite3Database::getIntArray(const String& sql, int col, std::vector<int>& result) {
		Sqlite3StatementPointer pst = prepare(sql);	
		getIntArray(*pst, col, result);
	}
	
	void Sqlite3Database::getIntArray(Sqlite3Statement& stat, int col, std::vector<int>& result) {
		while (stat.step() == SQLITE_ROW) {
			result.push_back(int(stat.getInt64(col)));
		}
	}

	void Sqlite3Database::getStringArray(const String& sql, int col, StringArray& result) {
		Sqlite3StatementPointer pst = prepare(sql);	
		getStringArray(*pst, col, result);
	}
	
	void Sqlite3Database::getStringArray(Sqlite3Statement& stat, int col, StringArray& result) {
		while (stat.step() == SQLITE_ROW) {
			result.push_back(stat.column(col));
		}
	}

	std::runtime_error Sqlite3Database::last_error() {
		const char* error = sqlite3_errmsg(pdb);
		return std::runtime_error(error);
	}

}
