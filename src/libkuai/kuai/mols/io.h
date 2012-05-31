// Basic Molecule IO function were defined in this file.
// For the real implementaions, such as SDF reader, PDB reader, etc 
// are defined in their files.

#include <fstream>

#include <kuai/typedef.h>
#include <kuai/tools/iotools.h>
#include <kuai/mols/Molecule.h>

#ifndef _KUAI_MOLS_IO_
#define _KUAI_MOLS_IO_

namespace kuai {

	typedef Reader<MoleculePtr> BasicMoleculeReader;

	class MoleculeReaderManager 
		: public ReaderManager<MoleculePtr>
	{
	private:
		MoleculeReaderManager();

	public:
		static MoleculeReaderManager& get_instance();
	};

	typedef Writer<MoleculePtr> BasicMoleculeWriter;

	class MoleculeWriterManager 
		: public WriterManager<MoleculePtr>
	{
	private:
		MoleculeWriterManager();

	public:
		static MoleculeWriterManager& get_instance();
	};

	typedef ReaderIterator<MoleculeReaderManager> MoleculeIterator;
}



#endif
