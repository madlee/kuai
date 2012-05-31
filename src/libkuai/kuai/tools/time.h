#include <cstring>
#include <cassert>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <kuai/typedef.h>


#ifndef _KUAI_TIME_H_
#define _KUAI_TIME_H_

namespace kuai {
	typedef boost::posix_time::ptime Time;


	inline Time now() {
		return boost::posix_time::microsec_clock::local_time();
	}

	inline Time utc_now() {
		return boost::posix_time::microsec_clock::universal_time();
	}
}
#endif
