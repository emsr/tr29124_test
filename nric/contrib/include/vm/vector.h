#ifndef VECTOR_H
#define VECTOR_H

/******************************************************************************
 *
 *	vector.h
 *	
 *	Header file for Vector class declaration and definition.
 *
 *	This header file is intended to be included by vec_mat.h, and it depends
 *	on the following other header files:
 *	
 *	vec_mat.h  vm_traits.h  matrix.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

#ifdef _MSC_VER // Microsoft does not provide the random_access_iterator base-class.
namespace std {
template <class T, class pd>
struct random_access_iterator : public std::iterator<std::random_access_iterator_tag,
	T, pd> {};
}
#endif // _MSC_VER

template <class T> class Matrix;

// Vector class. Views a data object as a vector.
template <class T>
class Vector
{
public:
	typedef typename MemAlloc<T>::alloc_type     alloc_type;
	typedef typename alloc_type::value_type      value_type;
	typedef typename alloc_type::pointer         pointer;
	typedef typename alloc_type::const_pointer   const_pointer;
	typedef typename alloc_type::size_type       size_type;
	typedef typename alloc_type::difference_type difference_type;
	typedef typename VMTraits<T>::real_type      real_type;

protected:
	size_t		nn_;				// Size of vector.
	ptrdiff_t	stride_;			// Stride of vector.
	T			*ptr_;				// Pointer to data for indexing.
	T			*data_;				// Pointer to the data object.
	size_t		*refs_;				// Pointer to the number of references.
	size_t		data_size_;			// Size of data object.
	alloc_type  alloc_;				// Memory allocator.
				
public:
	// Constructors and destructor.
	Vector();
	explicit Vector(size_t n);
	Vector(const NoInit_& noinit, size_t n);
	Vector(const T& a, size_t n);
	Vector(const T* a, size_t n, ptrdiff_t s = 1);
	Vector(const Vector& rhs);
	Vector(const T& strt, const T& step, size_t n);
	Vector(const Vector<real_type>& re, const Vector<real_type>& im);
	Vector(const Vector<real_type>& re, const real_type& im);
	~Vector();

	// Overloaded operators.
	Vector& operator=(const Vector& rhs);
	Vector& operator=(const T& a);
	Vector& operator+=(const Vector& rhs);
	Vector& operator+=(const T& a);
	Vector& operator-=(const Vector& rhs);
	Vector& operator-=(const T& a);
	Vector& operator*=(const Vector& rhs);
	Vector& operator*=(const T& a);
	Vector& operator/=(const Vector& rhs);
	Vector& operator/=(const T& a);
	Vector& operator%=(const Vector& rhs);
	Vector& operator%=(const T& a);
	Vector& operator&=(const Vector& rhs);
	Vector& operator&=(const T& a);
	Vector& operator^=(const Vector& rhs);
	Vector& operator^=(const T& a);
	Vector& operator|=(const Vector& rhs);
	Vector& operator|=(const T& a);
	Vector& operator>>=(const Vector& rhs);
	Vector& operator>>=(const T& a);
	Vector& operator<<=(const Vector& rhs);
	Vector& operator<<=(const T& a);
	Vector operator+() const;
	Vector operator-() const;
	Vector<bool> operator!() const;
	Vector operator~() const;
	T& operator[](size_t i);
	const T& operator[](size_t i) const;

	// Iterator forward declaration.
	class iterator;
	class const_iterator;

	// Member functions.
	Vector& apply(T (*fn)(T x));
	Vector& apply(T (*fn)(const T& x));
	bool isdeep() const;
	void makedeep();
	void reshape(size_t n);
	void free();
	void resize(size_t n);
	Vector copy() const;
	void reference(const Vector& rhs);
	Vector slice(size_t b, size_t n, ptrdiff_t s = 1) const;
	Vector<real_type> real() const;
	Vector<real_type> imag() const;
	Vector<T> rotate(ptrdiff_t n) const;
	Vector<T> delta() const;
	Vector<T> cumsum() const;
	size_t size() const;
	T* data();
	const T* data() const;
	ptrdiff_t stride() const;
	T sum() const;
	T min() const;
	T max() const;
	Matrix<T> matrix(size_t b, size_t n, size_t m, ptrdiff_t rs, ptrdiff_t cs) const;
	Matrix<T> rowmat() const;
	Matrix<T> colmat() const;
	void sort();
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	void read(std::istream& in);
	void read(std::istream& in, size_t ncount);
	void write(std::ostream& out) const;

	template <class G>
	Vector<T>& fill(const G& gen)
	{
		T *p1 = ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) = gen();
			p1 += stride_;
		}
		return *this;
	}

	template <class G>
	Vector& fill2(const G& fn)
	{
		T *p1 = ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) = fn(*p1);
			p1 += stride_;
		}
		return *this;
	}

	template <class U>
	Vector<U> convert(const U& a) const
	{
		Vector<U> tmp(nn_);
		T *p2 = ptr_;
		U *p1 = tmp.data();
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) = static_cast<U>(*p2);
			++p1;
			p2 += stride_;
		}
		return tmp;
	}

	template <class U>
	Vector<U> shallow_cast(const U& a) const
	{
		if (sizeof(U) != sizeof(T))
			vmerror("Size mismatch error in Vector::shallow_cast().");
		Vector<U> tmp;
		Vector<T> *tp    = reinterpret_cast<Vector<T> *>(&tmp);
		(*tp).nn_        = nn_;
		(*tp).stride_    = stride_;
		(*tp).ptr_       = ptr_;
		(*tp).data_      = data_;
		(*tp).refs_      = refs_;
		(*tp).data_size_ = data_size_;
		if (data_) ++(*refs_);
		return tmp;
	}

	friend class Matrix<T>;

	// Iterator class.
	class iterator : public std::random_access_iterator<T, ptrdiff_t>
	{
	protected:
		T			*ptr_;			// Pointer for indexing.
		ptrdiff_t	stride_;		// Stride of data.

	public:
		explicit iterator(T* p = 0, ptrdiff_t s = 1) : ptr_(p), stride_(s) {;}
		iterator(const iterator& rhs) : ptr_(rhs.ptr_), stride_(rhs.stride_) {;}
		~iterator() {;}
		iterator& operator++()
		{
			ptr_ += stride_;
			return *this;
		}
		iterator& operator--()
		{
			ptr_ -= stride_;
			return *this;
		}
		iterator operator++(int)
		{
			iterator tmp = *this;
			ptr_ += stride_;
			return tmp;
		}
		iterator operator--(int)
		{
			iterator tmp = *this;
			ptr_ -= stride_;
			return tmp;
		}
		iterator& operator+=(ptrdiff_t d)
		{
			ptr_ += d * stride_;
			return *this;
		}
		iterator& operator-=(ptrdiff_t d)
		{
			ptr_ -= d * stride_;
			return *this;
		}
		iterator& operator=(const iterator& iter)
		{
			stride_ = iter.stride_;
			ptr_    = iter.ptr_;
			return *this;
		}
		ptrdiff_t operator-(const iterator& iter) const
		{
			return (ptr_ - iter.ptr_) / stride_;
		}
		iterator operator+(ptrdiff_t n) const
		{
			iterator tmp = *this;
			tmp += n;
			return tmp;
		}
		iterator operator-(ptrdiff_t n) const
		{
			iterator tmp = *this;
			tmp -= n;
			return tmp;
		}
		T& operator[](size_t n) const
		{
			return ptr_[n * stride_];
		}
		T& operator*() const
		{
			return *ptr_;
		}
		T *operator->() const
		{
			return ptr_;
		}
		bool operator==(const iterator& iter) const
		{
			return (iter.ptr_ == ptr_);
		}
		bool operator!=(const iterator& iter) const
		{
			return (iter.ptr_ != ptr_);
		}
		bool operator<(const iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ < iter.ptr_) : (iter.ptr_ < ptr_);
		}
		bool operator>(const iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ > iter.ptr_) : (iter.ptr_ > ptr_);
		}
		bool operator<=(const iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ <= iter.ptr_) : (iter.ptr_ <= ptr_);
		}
		bool operator>=(const iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ >= iter.ptr_) : (iter.ptr_ >= ptr_);
		}
		friend class const_iterator;
	};

	// Const_iterator class.
	class const_iterator : public std::random_access_iterator<T, ptrdiff_t>
	{
	protected:
		T			*ptr_;			// Pointer for indexing.
		ptrdiff_t	stride_;		// Stride of data.

	public:
		explicit const_iterator(T* p = 0, ptrdiff_t s = 1) : ptr_(p), stride_(s) {;}
		const_iterator(const const_iterator& rhs) : ptr_(rhs.ptr_), stride_(rhs.stride_) {;}
		const_iterator(const iterator& rhs) : ptr_(rhs.ptr_), stride_(rhs.stride_) {;}
		~const_iterator() {}
		const_iterator& operator++()
		{
			ptr_ += stride_;
			return *this;
		}
		const_iterator& operator--()
		{
			ptr_ -= stride_;
			return *this;
		}
		const_iterator operator++(int)
		{
			const_iterator tmp = *this;
			ptr_ += stride_;
			return tmp;
		}
		const_iterator operator--(int)
		{
			const_iterator tmp = *this;
			ptr_ -= stride_;
			return tmp;
		}
		const_iterator& operator+=(ptrdiff_t d)
		{
			ptr_ += d * stride_;
			return *this;
		}
		const_iterator& operator-=(ptrdiff_t d)
		{
			ptr_ -= d * stride_;
			return *this;
		}
		const_iterator& operator=(const const_iterator& iter)
		{
			stride_ = iter.stride_;
			ptr_    = iter.ptr_;
			return *this;
		}
		const_iterator& operator=(const iterator& iter)
		{
			stride_ = iter.stride_;
			ptr_    = iter.ptr_;
			return *this;
		}
		ptrdiff_t operator-(const const_iterator& iter) const
		{
			return (ptr_ - iter.ptr_) / stride_;
		}
		const_iterator operator+(ptrdiff_t n) const
		{
			const_iterator tmp = *this;
			tmp += n;
			return tmp;
		}
		const_iterator operator-(ptrdiff_t n) const
		{
			const_iterator tmp = *this;
			tmp -= n;
			return tmp;
		}
		const T& operator[](size_t n) const
		{
			return ptr_[n * stride_];
		}
		const T& operator*() const
		{
			return *ptr_;
		}
		const T *operator->() const
		{
			return ptr_;
		}
		bool operator==(const const_iterator& iter) const
		{
			return (iter.ptr_ == ptr_);
		}
		bool operator!=(const const_iterator& iter) const
		{
			return (iter.ptr_ != ptr_);
		}
		bool operator<(const const_iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ < iter.ptr_) : (iter.ptr_ < ptr_);
		}
		bool operator>(const const_iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ > iter.ptr_) : (iter.ptr_ > ptr_);
		}
		bool operator<=(const const_iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ <= iter.ptr_) : (iter.ptr_ <= ptr_);
		}
		bool operator>=(const const_iterator& iter) const
		{
			return (stride_ > 0) ? (ptr_ >= iter.ptr_) : (iter.ptr_ >= ptr_);
		}
	};
};


// Constructors.
template <class T>
inline Vector<T>::Vector() : nn_(0), stride_(1), ptr_(0), data_(0), refs_(0), data_size_(0) {;}

template <class T>
Vector<T>::Vector(size_t n) : nn_(n), stride_(1), ptr_(0), data_(0), refs_(0), data_size_(n)
{
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_) vmerror("Memory allocation error in Vector constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		T zero = Zero<T>();
		T *p1 = ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			alloc_.construct(p1, zero);
			++p1;
		}
	}
}

template <class T>
Vector<T>::Vector(const NoInit_& noinit, size_t n) : nn_(n), stride_(1), ptr_(0), data_(0),
	refs_(0), data_size_(n)
{
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_) vmerror("Memory allocation error in Vector constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		if (!VMTraits<T>::is_simple)
		{
			T *p1 = ptr_;
			for (size_t i = 0; i < nn_; ++i)
			{
				new ((void *)p1) T();
				++p1;
			}
		}
	}
}

template <class T>
Vector<T>::Vector(const T& a, size_t n) : nn_(n), stride_(1), ptr_(0), data_(0), refs_(0),
	data_size_(n)
{
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_) vmerror("Memory allocation error in Vector constructor.");
		ptr_  = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			alloc_.construct(p1, a);
			++p1;
		}
	}
}

template <class T>
Vector<T>::Vector(const T* a, size_t n, ptrdiff_t s) : nn_(n), stride_(1), ptr_(0), data_(0),
	refs_(0), data_size_(n)
{
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_) vmerror("Memory allocation error in Vector constructor.");
		ptr_  = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		const T *p2 = a;
		for (size_t i = 0; i < nn_; ++i)
		{
			alloc_.construct(p1, (*p2));
			++p1;
			p2 += s;
		}
	}
}

template <class T>
Vector<T>::Vector(const Vector<T>& rhs) :nn_(rhs.nn_), stride_(rhs.stride_), ptr_(rhs.ptr_),
	data_(rhs.data_), refs_(rhs.refs_), data_size_(rhs.data_size_)
{
	if (data_) ++(*refs_);
}

template <class T>
Vector<T>::Vector(const T& strt, const T& step, size_t n) : nn_(n), stride_(1), ptr_(0),
	data_(0), refs_(0), data_size_(n)
{
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_) vmerror("Memory allocation error in Vector constructor.");
		ptr_  = data_;
		refs_ = new size_t(1);
		T a = strt;
		T *p1 = ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			alloc_.construct(p1, a);
			++p1;
			a += step;
		}
	}
}

template <class T>
Vector<T>::Vector(const Vector<real_type>& re, const Vector<real_type>& im) :
	nn_(re.size()), stride_(1), ptr_(0), data_(0), refs_(0), data_size_(re.size())
{
	if (im.size() != nn_)
		vmerror("Size mismatch error in complex Vector constructor.");
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_)	vmerror("Memory allocation error in Vector constructor.");
		ptr_  = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		const real_type *p2 = re.data(), *p3 = im.data();
		ptrdiff_t s2 = re.stride(), s3 = im.stride();
		for (size_t i = 0; i < nn_; ++i)
		{
			alloc_.construct(p1, T(*p2, *p3));
			++p1;
			p2 += s2;
			p3 += s3;
		}
	}
}

template <class T>
Vector<T>::Vector(const Vector<real_type>& re, const real_type& im) :
	nn_(re.size()), stride_(1), ptr_(0), data_(0), refs_(0), data_size_(re.size())
{
	if (nn_)
	{
		data_ = alloc_.allocate(nn_, 0);
		if (!data_) vmerror("Memory allocation error in Vector constructor.");
		ptr_  = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		const real_type *p2 = re.data();
		ptrdiff_t s2 = re.stride();
		for (size_t i = 0; i < nn_; ++i)
		{
			alloc_.construct(p1, T(*p2, im));
			++p1;
			p2 += s2;
		}
	}
}

// Destructor.
template <class T>
inline Vector<T>::~Vector()
{
	if (data_)
	{
		if ((*refs_) == 1)
		{
			if (!VMTraits<T>::is_simple)
			{
				T *p1 = data_;
				for (size_t i = 0; i < data_size_; ++i)
				{
					alloc_.destroy(p1);
					++p1;
				}
			}
			alloc_.deallocate(data_, data_size_);
			delete refs_;
		}
		else
			--(*refs_);
	}
}


// Overloaded operators.
template <class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& rhs)
{
	if (nn_ == rhs.nn_)
	{
		if (this == &rhs) return *this;
		if (data_ == rhs.data_) return operator=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) = (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator+=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) += (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator-=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) -= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator*=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator*=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) *= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator/=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator/=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) /= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator%=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator%=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) %= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator&=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator&=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) &= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator^=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator^=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) ^= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator|=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator|=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) |= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator>>=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator>>=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) >>= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator<<=(const Vector<T>& rhs)
{
	size_t n = rhs.nn_;
	if (n == nn_)
	{
		if (data_ == rhs.data_) return operator<<=(rhs.copy());
		T *p1 = ptr_;
		T *p2 = rhs.ptr_;
		for (size_t i = 0; i < nn_; ++i)
		{
			(*p1) <<= (*p2);
			p1 += stride_;
			p2 += rhs.stride_;
		}
	}
	else
		vmerror("Size mismatch error in Vector assignment operator.");
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator+=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) += a;
		p1 += stride_;
	}
	return *this;
}


template <class T>
Vector<T>& Vector<T>::operator-=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) -= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator*=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) *= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator/=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) /= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator%=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) %= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator&=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) &= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator|=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) |= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator>>=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) >>= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::operator<<=(const T& a)
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) <<= a;
		p1 += stride_;
	}
	return *this;
}

template <class T>
inline Vector<T> Vector<T>::operator+() const
{
	Vector<T> result = copy();
	T *p1 = result.ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = +(*p1);
		p1++;
	}
	return result;
}

template <class T>
Vector<T> Vector<T>::operator-() const
{
	Vector<T> result = copy();
	T *p1 = result.ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = -(*p1);
		p1++;
	}
	return result;
}

template <class T>
Vector<bool> Vector<T>::operator!() const
{
	Vector<bool> result(nn_);
	bool *p1 = result.data();
	T    *p2 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = !(*p2);
		p1++;
		p2 += stride_;
	}
	return result;
}

template <class T>
Vector<T> Vector<T>::operator~() const
{
	Vector<T> result = copy();
	T *p1 = result.ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = ~(*p1);
		p1++;
	}
	return result;
}

template <class T>
inline T& Vector<T>::operator[](size_t i)
{
#ifdef VM_BOUNDS_CHECK
	if (i >= nn_)
		vmerror("Index out of bounds error.");
#endif // VM_BOUNDS_CHECK
	return ptr_[i * stride_];
}

template <class T>
inline const T& Vector<T>::operator[](size_t i) const
{
#ifdef VM_BOUNDS_CHECK
	if (i >= nn_)
		vmerror("Index out of bounds error.");
#endif // VM_BOUNDS_CHECK
	return ptr_[i * stride_];
}


// Member functions.
template <class T>
Vector<T>& Vector<T>::apply(T (*fn)(T x))
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = fn(*p1);
		p1 += stride_;
	}
	return *this;
}

template <class T>
Vector<T>& Vector<T>::apply(T (*fn)(const T& x))
{
	T *p1 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = fn(*p1);
		p1 += stride_;
	}
	return *this;
}

template <class T>
inline bool Vector<T>::isdeep() const
{
	if (!data_) return true;
	if (*refs_ == 1 && stride_ == 1 && data_size_ == nn_)
		return true;
	else return false;
}

template <class T>
inline void Vector<T>::makedeep()
{
	if (!isdeep()) reference(copy());
}

template <class T>
inline void Vector<T>::reshape(size_t n)
{
	free();
	if (n)
	{
		Vector<T> tmp(n);
		reference(tmp);
	}
}

template <class T>
inline void Vector<T>::free()
{
	Vector<T> tmp;
	reference(tmp);
}

template <class T>
void Vector<T>::resize(size_t n)
{
	if (n)
	{
		T *new_data = alloc_.allocate(n, 0);
		if (!new_data) vmerror("Memory allocation error in Vector::resize().");
		T *p1 = new_data, *p2 = ptr_;
		ptrdiff_t s = stride_;
		size_t i;
		size_t mn = nn_ < n ? nn_ : n;
		for (i = 0; i < mn; ++i)
		{
			alloc_.construct(p1, (*p2));
			++p1;
			p2 += s;
		}
		T zero = Zero<T>();
		for (i = mn; i < n; ++i)
		{
			alloc_.construct(p1, zero);
			++p1;
		}
		if (data_)
		{
			if (*refs_ == 1)
			{
				if (!VMTraits<T>::is_simple)
				{
					p1 = data_;
					for (i = 0; i < data_size_; ++i)
					{
						alloc_.destroy(p1);
						++p1;
					}
				}
				alloc_.deallocate(data_, data_size_);
			}
			else
			{
				--(*refs_);
				refs_ = new size_t(1);
			}
		}
		else refs_ = new size_t(1);

		data_      = new_data;
		ptr_       = data_;
		nn_        = n;
		stride_    = 1;
		data_size_ = n;
	}
	else free();
}

template <class T>
Vector<T> Vector<T>::copy() const
{
	Vector<T> tmp;
	if (data_)
	{
		alloc_type a2;
		T *new_data = a2.allocate(nn_, 0);
		if (!new_data) vmerror("Memory allocation error in Vector::copy().");
		T *p1 = new_data;
		const T *p2 = ptr_;
		ptrdiff_t s = stride_;
		size_t i;
		for (i = 0; i < nn_; ++i)
		{
			a2.construct(p1, (*p2));
			++p1;
			p2 += s;
		}
		tmp.data_      = new_data;
		tmp.ptr_       = new_data;
		tmp.nn_        = nn_;
		tmp.stride_    = 1;
		tmp.data_size_ = nn_;
		tmp.refs_ = new size_t(1);
	}
	return tmp;
}

template <class T>
void Vector<T>::reference(const Vector<T>& rhs)
{
	if (this != &rhs)
	{
		if (data_)
		{
			if (*refs_ == 1)
			{
				if (!VMTraits<T>::is_simple)
				{
					T *p1 = data_;
					for (size_t i = 0; i < data_size_; ++i)
					{
						alloc_.destroy(p1);
						++p1;
					}
				}
				alloc_.deallocate(data_, data_size_);
				delete refs_;
			}
			else
				--(*refs_);
		}
		data_      = rhs.data_;
		refs_      = rhs.refs_;
		nn_        = rhs.nn_;
		stride_    = rhs.stride_;
		ptr_       = rhs.ptr_;
		data_size_ = rhs.data_size_;
		if (data_) ++(*refs_);
	}
}

template <class T>
inline Vector<T> Vector<T>::slice(size_t b, size_t n, ptrdiff_t s) const
{
	if (b >= nn_ ||
		b + (n - 1) * s < 0 ||
		b + (n - 1) * s >= nn_)
		vmerror("Index out of bounds error.");
	Vector<T> tmp = *this;
	tmp.nn_     = n;
	tmp.stride_ = s * stride_;
	tmp.ptr_    = ptr_ + b * stride_;
	return tmp;
}

template <class T>
Vector<VMTraits<T>::real_type> Vector<T>::real() const
{
	Vector<real_type> tmp(nn_);
	real_type *p1 = tmp.data();
	T  *p2 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = (*p2).real();
		++p1;
		p2 += stride_;
	}
	return tmp;
}

template <class T>
Vector<VMTraits<T>::real_type> Vector<T>::imag() const
{
	Vector<real_type> tmp(nn_);
	real_type *p1 = tmp.data();
	T  *p2 = ptr_;
	for (size_t i = 0; i < nn_; ++i)
	{
		(*p1) = (*p2).imag();
		++p1;
		p2 += stride_;
	}
	return tmp;
}

template <class T>
Vector<T> Vector<T>::rotate(ptrdiff_t n) const
{
	if (!nn_) vmerror("Rotation of zero size vector not allowed.");
	if (abs(n) >= nn_) vmerror("Index out of bounds in Vector::rotate() function.");
	Vector<T> tmp(nn_);
	int i, j;
	for (i = 0; i < nn_; ++i)
	{
		j = i + n;
		if (j < 0) j += nn_;
		else if (j >= nn_) j -= nn_;
		tmp[i] = (*this)[j];
	}
	return tmp;
}

template <class T>
Vector<T> Vector<T>::delta() const
{
	if (!nn_) vmerror("Delta of zero size vector not allowed.");
	Vector<T> tmp(nn_);
	for (size_t i = 0; i < nn_ - 1; ++i)
		tmp[i] = (*this)[i + 1] - (*this)[i];
	tmp[nn_ - 1] = (*this)[0] - (*this)[nn_ - 1];
	return tmp;
}

template <class T>
Vector<T> Vector<T>::cumsum() const
{
	if (!nn_) vmerror("Cumsum of zero size vector not allowed.");
	Vector<T> tmp = copy();
	for (size_t i = 1; i < nn_; ++i)
		tmp[i] += tmp[i - 1];
	return tmp;
}

template <class T>
inline size_t Vector<T>::size() const
{
	return nn_;
}

template <class T>
inline T* Vector<T>::data()
{
	return ptr_;
}

template <class T>
inline const T* Vector<T>::data() const
{
	return ptr_;
}

template <class T>
inline ptrdiff_t Vector<T>::stride() const
{
	return stride_;
}

template <class T>
inline Matrix<T> Vector<T>::matrix(size_t b, size_t n, size_t m, ptrdiff_t rs,
	ptrdiff_t cs) const
{
	if (b >= nn_ ||
		b + (n - 1) * rs < 0 ||
		b + (n - 1) * rs >= nn_ ||
		b + (m - 1) * cs < 0 ||
		b + (m - 1) * cs >= nn_ ||
		b + (n - 1) * rs + (m - 1) * cs < 0 ||
		b + (n - 1) * rs + (m - 1) * cs >= nn_)
		vmerror("Index out of bounds error.");
	Matrix<T> tmp;
	tmp.data_      = data_;
	tmp.refs_      = refs_;
	tmp.ptr_       = ptr_ + b * stride_;
	tmp.nn_        = n;
	tmp.mm_        = m;
	tmp.rstride_   = stride_ * rs;
	tmp.cstride_   = stride_ * cs;
	tmp.data_size_ = data_size_;
	if (data_) ++(*refs_);
	return tmp;
}

template <class T>
inline Matrix<T> Vector<T>::rowmat() const
{
	Matrix<T> tmp;
	tmp.data_      = data_;
	tmp.refs_      = refs_;
	tmp.ptr_       = ptr_;
	tmp.nn_        = 1;
	tmp.mm_        = nn_;
	tmp.rstride_   = nn_ * stride_;
	tmp.cstride_   = stride_;
	tmp.data_size_ = data_size_;
	if (data_) ++(*refs_);
	return tmp;
}

template <class T>
inline Matrix<T> Vector<T>::colmat() const
{
	Matrix<T> tmp;
	tmp.data_      = data_;
	tmp.refs_      = refs_;
	tmp.ptr_       = ptr_;
	tmp.nn_        = nn_;
	tmp.mm_        = 1;
	tmp.rstride_   = stride_;
	tmp.cstride_   = nn_ * stride_;
	tmp.data_size_ = data_size_;
	if (data_) ++(*refs_);
	return tmp;
}

template <class T>
inline void Vector<T>::sort()
{
	if (!nn_) return;
	if (stride_ == 1) std::sort(ptr_, ptr_ + nn_);
	else std::sort(begin(), end());
}

template <class T>
T Vector<T>::sum() const
{
	if (!nn_) vmerror("Sum of zero size Vector not allowed.");
	T tmp = *ptr_;
	T *p1 = ptr_;
	for (size_t i = 1; i < nn_; ++i)
	{
		p1  += stride_;
		tmp += (*p1);
	}
	return tmp;
}

template <class T>
T Vector<T>::min() const
{
	if (!nn_) vmerror("Minimum of zero size Vector not allowed.");
	T tmp = *ptr_;
	T *p1 = ptr_;
	for (size_t i = 1; i < nn_; ++i)
	{
		p1 += stride_;
		if ((*p1) < tmp) tmp = (*p1);
	}
	return tmp;
}

template <class T>
T Vector<T>::max() const
{
	if (!nn_) vmerror("Maximum of zero size Vector not allowed.");
	T tmp = *ptr_;
	T *p1 = ptr_;
	for (size_t i = 1; i < nn_; ++i)
	{
		p1 += stride_;
		if (tmp < (*p1)) tmp = (*p1);
	}
	return tmp;
}

template <class T>
inline Vector<T>::iterator Vector<T>::begin()
{
	return Vector<T>::iterator(ptr_, stride_);
}

template <class T>
inline Vector<T>::const_iterator Vector<T>::begin() const
{
	return Vector<T>::const_iterator(ptr_, stride_);
}

template <class T>
inline Vector<T>::iterator Vector<T>::end()
{
	return Vector<T>::iterator(ptr_ + nn_ * stride_, stride_);
}

template <class T>
inline Vector<T>::const_iterator Vector<T>::end() const
{
	return Vector<T>::const_iterator(ptr_ + nn_ * stride_, stride_);
}

template <class T>
void Vector<T>::write(std::ostream& out) const
{
	for (size_t i = 0; i < nn_; ++i)
	{
		out << ptr_[i * stride_];
		if (i < nn_ - 1) out << '\n';
	}
}

template <class T>
void Vector<T>::read(std::istream& in)
{
	T tmp;
	size_t buff_size = 1024, i = 0;
	reshape(buff_size);
	while (in >> tmp)
	{
		ptr_[i] = tmp;
		++i;
		if (i == buff_size)
		{
			buff_size *= 2;
			resize(buff_size);
		}
	}
	if (i < buff_size) resize(i);
}

template <class T>
void Vector<T>::read(std::istream& in, size_t ncount)
{
	T tmp;
	size_t i = 0;
	reshape(ncount);
	if (!ncount) return;
	while (in >> tmp && i < ncount)
	{
		ptr_[i] = tmp;
		++i;
	}
	if (i < ncount) resize(i);
}

#endif // VECTOR_H
