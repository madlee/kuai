#include <cstring>
#include <cassert>

#include <kuai/typedef.h>


#ifndef _KUAI_THREAD_H_
#define _KUAI_THREAD_H_

	#ifndef USE_DUMN_THREAD
		#include <boost/thread.hpp>
		namespace kuai { 
			typedef boost::thread Thread;
			typedef boost::mutex  Mutex;
			typedef boost::thread_group ThreadGroup;
			typedef boost::barrier Barrier;

			inline Index hardware_concurrency() {
				return boost::thread::hardware_concurrency();
			}

			struct CalculationResource { 
			private:
				CalculationResource() { 
					_nThreads = 1;
				}

			public:
				Index countProcess() const { 
					return 1;
				}

				Index currentProcess() const { 
					return 1;
				}

				Index countThread() const { 
					return _nThreads;
				}

				void setThread(Index thread);

			private:
				Index _nThreads;

			public:
				static CalculationResource instance;
			};

			typedef boost::thread_interrupted ThreadInterrupted;
		}
	#else
		#include <list>
		#include <boost/utility.hpp>

		namespace kuai { 
			class Thread 
				: private boost::noncopyable
			{
			public:
				template<typename Function>
					explicit Thread(Function func)
				{ 
					func();
				}

				bool joinable() const { 
					return true;
				}
				void join() 
				{ }
		
				void detach() 
				{ }
			};
			class Mutex
				: private boost::noncopyable
			{ 
			public:
				void lock() { }
				bool try_lock() {
					return true;
				}
				void unlock() 
				{ }
			};

			class ThreadGroup
				: private boost::noncopyable
			{
			public:
				ThreadGroup()
				{ }
				~ThreadGroup()
				{ }
				template<typename F>
					void create_thread(F threadfunc)
				{ 
					threadfunc();
				}
				void add_thread(Thread* thrd)
				{  }
				void remove_thread(Thread* thrd)
				{  }
				void join_all()
				{ }
				void interrupt_all() 
				{ }
				size_t size() const { 
					return 1;
				}
			};

			class Barrier {
			public:
				explicit Barrier(Index count)
				{
					assert (count == 1);
				};

				~Barrier()
				{ }

				bool wait() 
				{
					return true;
				}
			};

			inline Index hardware_concurrency() {
				return 1;
			}

		}	// namespace kuai { 

	#endif	/* USE_DUMN_THREAD */



#endif	/* _KUAI_THREAD_HPP_ */

