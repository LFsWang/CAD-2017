#include<cstdlib>
std::size_t map_allocator_max;

template<typename T>
struct map_allocator{
	// typedefs
	typedef T value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;
	// convert an map_allocator<T> to map_allocator<U>
	
	pointer mem,top;
	
	template<typename U>
	struct rebind {
		typedef map_allocator<U> other;
	};
	explicit map_allocator() {
		top = mem = reinterpret_cast<pointer>(::operator new(map_allocator_max * sizeof (T)));
	}
	~map_allocator() {
		delete [] mem;
	}
	
	// address
	pointer address(reference r) { return &r; }
	const_pointer address(const_reference r) { return &r; }
	
	pointer allocate(size_type cnt, typename std::allocator<void>::const_pointer = 0) {
		//return reinterpret_cast<pointer>(::operator new(cnt * sizeof (T)));
		//if(cnt!=1) std::cerr<<"allocate error\n";
		return top++;
	}
	void deallocate(pointer p, size_type) { }
	
	void construct(pointer p, const T& t) { new(p) T(t); }
	void destroy(pointer p) { p->~T(); }
};