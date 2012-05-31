#include <memory>

#ifndef _KUAI_MEMTOOL_H_
#define _KUAI_MEMTOOL_H_


namespace kuai
{


//  owner_ptr mimics a built-in pointer except that it guarantees deletion
//  of the object pointed to, either on destruction of the owner_ptr or via
//  an explicit reset(). owner_ptr is a simple solution for simple needs;
//  use shared_ptr if your needs are more complex.

template<class T> 
	class owner_ptr 
		: boost::noncopyable // noncopyable
{
public:
    typedef T element_type;
	typedef owner_ptr<T> this_type;

    explicit owner_ptr( T * p = 0, bool owned=true )
		: px(p) // never throws
    {
		if (owned) {
			py = px;
		}
		else {
			py = NULL;
		}
    }

    ~owner_ptr() // never throws
    {
		delete py;
    }

    void reset(T * p = 0) // never throws
    {
        assert( p == 0 || p != px ); // catch self-reset errors
        this_type(p).swap(*this);
    }

    T & operator*() const // never throws
    {
        assert( px != 0 );
        return *px;
    }

    T * operator->() const // never throws
    {
        assert( px != 0 );
        return px;
    }

    T * get() const // never throws
    {
        return px;
    }

// implicit conversion to "bool"
#include <boost/smart_ptr/detail/operator_bool.hpp>

    void swap(owner_ptr & b) // never throws
    {
		std::swap(px, b.px);
		std::swap(py, b.py);
    }

private:
    T * px;
	T * py;
};

template<class T> inline void swap(owner_ptr<T> & a, owner_ptr<T> & b) // never throws
{
    a.swap(b);
}

// get_pointer(p) is a generic way to say p.get()

template<class T> inline T * get_pointer(owner_ptr<T> const & p)
{
    return p.get();
}

} // namespace kuai

#endif // #ifndef _KUAI_MEMTOOL_H_
