#include <kuai/tools/Logger.h>
#include <kuai/mols/io.h>
#include <kuai/mols/inchi.h>
#include <kuai/mols/smiles.h>

namespace kuai {

	MoleculeReaderManager::MoleculeReaderManager() {
		put(".smi", shared_reader<SmilesReader>);
		put(".inchi", create_reader<InchiReader>);
	}
	MoleculeReaderManager& MoleculeReaderManager::get_instance() {
		static MoleculeReaderManager instance;
		return instance;
	}

}
