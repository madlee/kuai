#include <fstream>
#include <sstream>
#include <iterator>
#include <boost/optional.hpp>
#include <kuai/typedef.h>
#include <kuai/tools/strtool.h>
#include <kuai/tools/memtool.h>
#include <kuai/tools/error.h>
#include <kuai/tools/FileName.h>
#include <kuai/tools/FunctionManager.h>
#include <kuai/tools/thread.h>



#ifndef _KUAI_TOOL_IOTOOL_H_
#define _KUAI_TOOL_IOTOOL_H_

namespace kuai { 

	std::fstream* try_open(const FileName& filename, std::ios_base::openmode flag);

	template<typename DataType>
		class Reader
	{
	public:
		typedef DataType Value;
		typedef DataType& Reference;

	public:
		virtual ~Reader()
		{ }
		
		/// Implements this function for read the data. Return true if the read is success,
		/// return false if the end of stream is met and throw an exception if some error
		/// occoured. Try to locate the stream to the start of the next record if possible
		/// so that user may ignore the error and read the next record.
		/// Try to lock the mutex when you read from stream and free the lock when you 
		/// finish the read.
		virtual bool read(std::istream& stream, Mutex& mutex, const FileName* filename, size_t count, Reference data) = 0;

		/// Open a file for read. In default it open the file in text mode. You may change
		/// it in derived class.
		virtual std::istream* try_open(const FileName& filename) const {
			return kuai::try_open(filename, std::ios_base::in);
		}
	};

	template<typename ReaderType> 
	inline SharedPtr<Reader<typename ReaderType::Value> > create_reader() 
	{
		return SharedPtr<Reader<typename ReaderType::Value> >(new ReaderType);
	}

	template<typename ReaderType> 
	inline SharedPtr<Reader<typename ReaderType::Value> > shared_reader() 
	{
		static SharedPtr<Reader<typename ReaderType::Value> > instance(new ReaderType);
		return instance;
	}

	template<typename DataType>
		class ReaderManager
			: public FunctionManager<String, SharedPtr<Reader<DataType> > (*)()>
	{ 
	public:
		typedef Reader<DataType> ReaderType;
		typedef SharedPtr<ReaderType>  ReaderPointer;
		typedef typename ReaderType::Value Value;
		typedef typename ReaderType::Reference Reference;
		typedef SharedPtr<Reader<DataType> > (*ReaderCreator)();

	public:
		virtual ReaderPointer get_reader(const String& key) const
		{
			ReaderPointer result;
			String type = boost::algorithm::to_lower_copy(key);
			ReaderCreator creator = get(type);
			if (creator == NULL) { 
				creator = get("");
			}
			if (creator) {
				result = ReaderPointer((*creator)());
			}

			return result;
		}
	};

	template<typename ReaderManager>
		class ReaderIterator
			: public std::iterator<std::forward_iterator_tag, typename ReaderManager::Value>, boost::noncopyable
	{
	public:
		typedef typename ReaderManager::Value Value;
		typedef typename ReaderManager::Reference Reference;
		typedef typename ReaderManager::ReaderType ReaderType;
		typedef typename ReaderManager::ReaderPointer ReaderPointer;

	public:
		explicit ReaderIterator(std::istream* stream, const String& type) 
			: _pstream(stream, false), _count(0)
		{
			_init_reader(type);
		}
		explicit ReaderIterator(const FileName& filename, const String& type) 
			: _count(0), _filename(filename)
		{
			_init_reader(type);
			_pstream.reset(_reader->try_open(filename));
		}
		explicit ReaderIterator(const FileName& filename) 
			: _count(0), _filename(filename)
		{
			_init_reader(filename.extname());
			_pstream.reset(_reader->try_open(filename));
		}

		/// This function try to read a data from a stream. Return true on success and return 
		/// false if it meet the end of the stream. If some other error happend it will throw 
		/// an exception to tell the reason. You can still read the next data when the error 
		/// occoured. 
		bool next(Reference data) {
			if (_pstream) {
				// clear before read so the previous error state will not effect this read.
				_pstream->clear(); 
				if (!read(*_pstream, _mutex, _filename.get_ptr(), ++_count, data)) {
					_pstream.reset();
				}
			}	

			return _pstream;
		}

		/// Return true if the end of the stream is meet. Try to avoid use of this function
		/// in multi-thread program.
		bool end() const {
			return !_pstream;
		}

		operator bool() const {
			return !end();
		}

		ReaderIterator& operator++() {
			next(_buffer);
			return *this;
		}
		pointer operator->() {
			return &_buffer;
		}


	private:
		ReaderPointer					_reader;
		Mutex							_mutex;
		owner_ptr<std::istream>			_pstream;
		boost::optional<FileName>		_filename;
		size_t							_count;
		Value							_buffer;

	private:
		void _init_reader(const String& type) {
			_reader = ReaderManager::get_instance().get_reader(type);
			assert (_reader);
		}
	};

	template<typename DataType>
		class Writer
	{
	public:
		typedef DataType Value;
		typedef const DataType& Reference;

	public:
		virtual ~Writer()
		{ }
		
		/// Implements this function for read the data. Return true if the read is success,
		/// return false if the end of stream is met and throw an exception if some error
		/// occoured. Try to locate the stream to the start of the next record if possible
		/// so that user may ignore the error and read the next record.
		/// Try to lock the mutex when you read from stream and free the lock when you 
		/// finish the read.
		virtual bool write(std::ostream& stream, Mutex& mutex, const FileName* filename, size_t count, Reference data) const = 0;

		/// Open a file for read. In default it open the file in text mode. You may change
		/// it in derived class.
		virtual std::ostream* try_open(const FileName& filename) const {
			return kuai::try_open(filename, std::ios_base::out);
		}
	};

	template<typename DataType>
		class WriterManager
			: public FunctionManager<String, Writer<DataType>*>
	{ 
	public:
		typedef const Writer<DataType>* WriterPointer;
		typedef DataType Value;
		typedef const DataType& Reference;

		virtual WriterPointer get_writer(const String& key) const
		{
			String type = boost::algorithm::to_lower_copy(key);
			if (WriterPointer result = get(type))
			{ 
				return result;
			}
			else
			{
				return get("");
			}
		}
	};



}

#endif
