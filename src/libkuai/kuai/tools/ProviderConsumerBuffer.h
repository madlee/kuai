#include <kuai/typedef.h>
#include <kuai/tools/threads.h>
#include <boost/circular_buffer.hpp>


#ifndef _KUAI_PROVIDER_CONSUMER_BUFFER_H_
#define _KUAI_PROVIDER_CONSUMER_BUFFER_H_

namespace kuai {

	template<typename ValueType>
	class ProviderConsumerBuffer {
	public:
		typedef boost::circular_buffer<ValueType>	Buffer;
		typedef boost::mutex::scoped_lock			Lock;
		typedef boost::condition_variable			Condition;

	public:
		explicit ProviderConsumerBuffer(size_t buffer_size) 
			: _buffer(buffer_size)
		{ }

		void push(const ValueType& v) {
			Lock lock(_monitor);
			while (_buffer.full()){
				_buffer_not_full.wait(lock);
			}
			_buffer.push_back(v);
			_buffer_not_empty.notify_one();
		}
		void pop(ValueType& result) {
			Lock lock(_monitor);
			while (_buffer.empty()){
				_buffer_not_empty.wait(lock);
			}
			swap(result, _buffer.back());
			_buffer.pop_back();
			_buffer_not_full.notify_one();
		}
		ValueType pop() {
			ValueType result;
			pop(result);
			return result;
		}

	private:
		Buffer			_buffer;
		Mutex			_monitor;
		Condition		_buffer_not_full, _buffer_not_empty;
	};
}

#endif

