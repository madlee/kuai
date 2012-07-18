#include <kuai/tools/time.h>
#include <kuai/mols/Molecule.h>
#include <kuai/mols/smiles.h>
#include <kuai/mols/inchi.h>
#include <kuai/mols/RingFinder.h>
#include "kuaichop.h"


using namespace kuai;

namespace {

	static const Char* SZZ_CREATE_DATABASE_SQL[] = {
		"CREATE TABLE IF NOT EXISTS Molecule (\n"
			"id INTEGER PRIMARY KEY AUTOINCREMENT,\n"
			"natoms INTEGER,\n"
			"nbonds INTEGER,\n"
			"formula CHAR(50),\n"
			"weight DOUBLE,\n"
			"inchi_key CHAR(50) UNIQUE,\n"
			"smiles TEXT DEFAULT ''\n"
		");\n",
		"CREATE INDEX IF NOT EXISTS index_molecule_natoms ON Molecule (natoms);\n",
		"CREATE INDEX IF NOT EXISTS index_molecule_nbonds ON Molecule (nbonds);\n",
		"CREATE INDEX IF NOT EXISTS index_molecule_weight ON Molecule (weight);\n",

		"CREATE TABLE IF NOT EXISTS Fragment (\n"
			"id INTEGER PRIMARY KEY AUTOINCREMENT,\n"
			"natoms INTEGER,\n"
			"nbonds INTEGER,\n"
			"inchi_key CHAR(50) UNIQUE,\n"
			"smiles TEXT DEFAULT ''\n"
		");\n",
		"CREATE INDEX IF NOT EXISTS index_fragment_natoms ON Fragment (natoms);\n",
		"CREATE INDEX IF NOT EXISTS index_fragment_nbonds ON Fragment (nbonds);\n",

		"CREATE TABLE IF NOT EXISTS FragInMol (\n"
			"mol_id INTEGER REFERENCES Molecule (id),\n"
			"frag_id INTEGER REFERENCES Fragment (id),\n"
			"remains INTEGER REFERENCES Fragment (id) \n"
			"UNIQUE (mol_id, frag_id, remains) \n"
		");\n",
		"CREATE INDEX IF NOT EXISTS index_fraginmol_mol_id ON FragInMol (mol_id);\n",
		"CREATE INDEX IF NOT EXISTS index_fraginmol_frag_id ON FragInMol (frag_id);\n",
		"CREATE INDEX IF NOT EXISTS index_fraginmol_remains ON FragInMol (remains);\n",

		"CREATE TABLE IF NOT EXISTS FragTree (\n"
			"parent INTEGER REFERENCES Fragment (id),\n"
			"child INTEGER REFERENCES Fragment (id) \n"
			"UNIQUE (parent, child) \n"
		");\n",
		"CREATE INDEX IF NOT EXISTS index_fragtree_parent ON FragTree (parent);\n",
		"CREATE INDEX IF NOT EXISTS index_fragtree_parent ON FragTree (child);\n",

		""
	};


	const Char* SZ_SQL_FIND_MOLECULE = "SELECT id FROM Molecule WHERE inchi_key = ?;"; 
	const Char* SZ_SQL_INSERT_MOLECULE = "INSERT INTO Molecule (natoms, nbonds, formula, weight, inchi_key, smiles) VALUES (?, ?, ?, ?, ?, ?);"; 


	boost::shared_ptr<Sqlite3Database> create_database(const FileName& database) {
		boost::shared_ptr<Sqlite3Database> result = Sqlite3Database::create(database);
		for (int i = 0; SZZ_CREATE_DATABASE_SQL[i][0]; ++i) {
			result->execute(SZZ_CREATE_DATABASE_SQL[i]);
		}
		result->commit();
		return result;
	}

	class KuaiChop {
	private:
		KuaiChop(Sqlite3RowID id, KuaiChopJob* job, MoleculePtr mol);
	
		Sqlite3RowID	_id;
		KuaiChopJob*	_job;
		MoleculePtr		_mol;

		Timer			_timer;

		void do_chop(KuaiChopJob* job, MoleculePtr mol, RingFinder& rings, Index cutted_bonds);

	public:
		static void do_chop(Sqlite3RowID id, KuaiChopJob* job, MoleculePtr mol);
	};

	KuaiChop::KuaiChop(Sqlite3RowID id, KuaiChopJob* job, MoleculePtr mol)
		: _id(id), _mol(mol), _job(job)
	{ }

	void KuaiChop::do_chop(Sqlite3RowID id, KuaiChopJob* job, MoleculePtr mol) {
		RingFinder rings(mol);
		KuaiChop chop(id, job, mol);
		chop.do_chop(job, mol, rings, 0);
	}

	void KuaiChop::do_chop(KuaiChopJob* job, MoleculePtr mol, RingFinder& rings, Index cutted_bonds) {

		if (job->cut_single_chain()) {
			Index nbonds = mol->count_bonds();
			for (Index i = 0; i < nbonds; ++i) {
				Bond* bond = mol->get_bond(i);
				if (rings.rank(bond) == RingFinder::NOT_ON_RING) {

				}
			}
		}
	}


	struct do_chop {
		do_chop(KuaiChopJob* job) 
			: _job(job)
		{ }
		void operator()() {
			for (;;) {
				String smiles = _job->pop();
				if (smiles.empty()) {
					break;
				}

				try {
					MoleculePtr mol = parse_smiles(smiles);
					Sqlite3RowID id = _job->add_molecule(mol);
					if (id >= 0) {
						_job->log_info("Molecule %s has been parsed. Start chop...", mol->name);
						KuaiChop::do_chop(id, _job, mol);
					}
					else {
						_job->log_warning("Molecule %s has been in the database. Skip it.", mol->name);
					}
				}
				catch (std::exception& error) {
					_job->log_error(error);
				}
			}
		}

	private:
		KuaiChopJob* _job;
	};
}

KuaiChopJob::KuaiChopJob(const char database[], size_t nThreads, size_t buffer_size)
	: _buffer(buffer_size)
{
	_database = create_database(FileName(database));
	_cut_single_chain = true;
	_cut_fused_rings = false;
	_cut_spiro_rings = false;
	_min_fragment_size = 3;
	_max_fragment_size = 16;
	_max_cuts = 4;

	for (size_t i = 0; i < nThreads; ++i) {
		_groups.add_thread(new Thread(do_chop(this)));
	}
}

KuaiChopJob::~KuaiChopJob() {
	close();
}


void KuaiChopJob::feed(const char v[]) {
	_buffer.push(v);
}


void KuaiChopJob::close() {
	for (size_t i = 0; i < _groups.size(); ++i) {
		_buffer.push("");
	}
	_groups.join_all();
}

kuai::Sqlite3RowID KuaiChopJob::add_molecule(MoleculePtr& molecule) {
	kuai::Sqlite3RowID result = -1;
	String inchi, key;

	boost::tie(inchi, key) = inchi_and_key(molecule, "/AuxNone");
	try {
		boost::mutex::scoped_lock lock(_access_db);
		Sqlite3StatementPointer sqlFindMolType = _database->prepare(SZ_SQL_FIND_MOLECULE);
		sqlFindMolType->bind(1, key);
		_database->getSingleInt(*sqlFindMolType);
	}
	catch (ErrorCanNotFindRow&) {
		_access_db.unlock();

		String formula = molecule->formula();
		RealNumber weight = molecule->weight();
		MoleculePtr unique_mol = parse_inchi(inchi);
		String usmiles = smiles(unique_mol);

		boost::mutex::scoped_lock lock(_access_db);
		Sqlite3StatementPointer sqlInsertMolecule = _database->prepare(SZ_SQL_INSERT_MOLECULE);
		sqlInsertMolecule->bind(1, int(molecule->count_atoms()));
		sqlInsertMolecule->bind(2, int(molecule->count_bonds()));
		sqlInsertMolecule->bind(3, formula);
		sqlInsertMolecule->bind(4, weight);
		sqlInsertMolecule->bind(5, key);
		sqlInsertMolecule->bind(6, usmiles);
		result = _database->insert(*sqlInsertMolecule);
	}

	return result;
}


void KuaiChopJob::log_error(const std::exception& error) {
	log_error(String(error.what()));
}

void KuaiChopJob::log_error(const kuai::String& error) {
	// TODO:
}
void KuaiChopJob::log_warning(const kuai::String& warning) {
	// TODO:
}
void KuaiChopJob::log_info(const kuai::String& info) {
	// TODO:
}
