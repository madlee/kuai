#include <kuai/typedef.h>

#ifndef _KUAI_TOOLS_FUNCTION_MANAGER_H_
#define _KUAI_TOOLS_FUNCTION_MANAGER_H_

namespace kuai {

	template<typename Key, typename Function>
	class FunctionManager
		: public KuaiMap<Key, Function> 
	{ 
	public:
		const Function get(const Key& v0) const {
			typename KuaiMap<Key, Function>::const_iterator it = this->find(v0);
			if (it != this->end()) {
				return it->second;
			}
			else {
				return NULL;
			}
		}

		void put(const Key& v0, Function p0) {
			insert(std::make_pair(v0, p0));
		}
	};

}

#endif 
