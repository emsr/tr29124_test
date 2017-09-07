#ifndef MATRIX_H
#define MATRIX_H

/******************************************************************************
 *
 *	matrix.h
 *	
 *	Header file for Matrix class declaration and definition.
 *
 *	This header file is intended to be included by vec_mat.h, and it depends
 *	on the following other header files:
 *	
 *	vec_mat.h  vm_traits.h  vector.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

// Matrix class. Views a data object as a matrix.
template <class T>
class Matrix
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
	size_t		nn_, mm_;			// Size of matrix.
	ptrdiff_t	rstride_;			// Distance between consecutive rows.
	ptrdiff_t	cstride_;			// Distance between consecutive columns.
	T			*ptr_;				// Pointer for indexing.
	T			*data_;				// Pointer to the data.
	size_t		*refs_;				// Pointer to number of references.
	size_t		data_size_;			// Size of data object.
	alloc_type  alloc_;				// Memory allocator.

public:
	// Constructors and destructor.
	Matrix();
	Matrix(size_t n, size_t m);
	Matrix(const NoInit_& noinit, size_t n, size_t m);
	Matrix(const T& a, size_t n, size_t m);
	Matrix(const T* a, size_t n, size_t m, ptrdiff_t s = 1);
	Matrix(const Matrix& rhs);
	Matrix(const Matrix<real_type>& re,	const Matrix<real_type>& im);
	Matrix(const Matrix<real_type>& re,	const real_type& im);
	~Matrix();

	// Overloaded operators.
	Matrix& operator=(const Matrix& rhs);
	Matrix& operator=(const T& a);
	Matrix& operator+=(const Matrix& rhs);
	Matrix& operator+=(const T& a);
	Matrix& operator-=(const Matrix& rhs);
	Matrix& operator-=(const T& a);
	Matrix& operator*=(const Matrix& rhs);
	Matrix& operator*=(const T& a);
	Matrix& operator/=(const Matrix& rhs);
	Matrix& operator/=(const T& a);
	Matrix& operator%=(const Matrix& rhs);
	Matrix& operator%=(const T& a);
	Matrix& operator&=(const Matrix& rhs);
	Matrix& operator&=(const T& a);
	Matrix& operator^=(const Matrix& rhs);
	Matrix& operator^=(const T& a);
	Matrix& operator|=(const Matrix& rhs);
	Matrix& operator|=(const T& a);
	Matrix& operator>>=(const Matrix& rhs);
	Matrix& operator>>=(const T& a);
	Matrix& operator<<=(const Matrix& rhs);
	Matrix& operator<<=(const T& a);
	Matrix operator+() const;
	Matrix operator-() const;
	Matrix<bool> operator!() const;
	Matrix operator~() const;
	Vector<T>::iterator operator[](size_t i);
	Vector<T>::const_iterator operator[](size_t i) const;
	T& operator()(size_t i, size_t j);
	const T& operator()(size_t i, size_t j) const;

	// Member functions.
	Matrix& apply(T (*fn)(T x));
	Matrix& apply(T (*fn)(const T& x));
	bool isdeep() const;
	void makedeep();
	void reshape(size_t n, size_t m);
	void free();
	void resize(size_t n, size_t m);
	void reference(const Matrix& rhs);
	Matrix copy() const;
	Matrix submatrix(size_t rb, size_t cb, size_t rn, size_t cn, ptrdiff_t rs, ptrdiff_t cs) const;
	Matrix transpose() const;
	Matrix<real_type> real() const;
	Matrix<real_type> imag() const;
	size_t nrows() const;
	size_t ncols() const;
	ptrdiff_t rstride() const;
	ptrdiff_t cstride() const;
	T* data();
	const T* data() const;
	Vector<T> row(size_t i) const;
	Vector<T> column(size_t j) const;
	Vector<T> unwrap() const;
	Vector<T>::iterator rbegin(size_t n);
	Vector<T>::const_iterator rbegin(size_t n) const;
	Vector<T>::iterator rend(size_t n);
	Vector<T>::const_iterator rend(size_t n) const;
	Vector<T>::iterator cbegin(size_t n);
	Vector<T>::const_iterator cbegin(size_t n) const;
	Vector<T>::iterator cend(size_t n);
	Vector<T>::const_iterator cend(size_t n) const;
	void read(std::istream& in);
	void read(std::istream& in, size_t ncount);
	void write(std::ostream& out) const;

	template <class G>
	Matrix& fill(const G& gen)
	{
		size_t i, j;
		T *p1;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) = gen();
				p1 += cstride_;
			}
		}
		return *this;
	}

	template <class G>
	Matrix& fill2(const G& fn)
	{
		size_t i, j;
		T *p1;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) = fn(*p1);
				p1 += cstride_;
			}
		}
		return *this;
	}

	template <class U>
	Matrix<U> convert(const U& a) const
	{
		Matrix<U> tmp(nn_, mm_);
		size_t i, j;
		T *p2;
		U *p1 = tmp.data();
		for (i = 0; i < nn_; ++i)
		{
			p2 = ptr_ + i * rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) = static_cast<U>(*p1);
				++p1;
				p2 += cstride_;
			}
		}
		return tmp;
	}

	template <class U>
	Matrix<U> shallow_cast(const U& a) const
	{
		if (sizeof(U) != sizeof(T))
			vmerror("Size mismatch error in Matrix::shallow_cast().");
		Matrix<U> tmp;
		Matrix<T> *tp    = reinterpret_cast<Matrix<T> *>(&tmp);
		(*tp).nn_        = nn_;
		(*tp).mm_        = mm_;
		(*tp).rstride_   = rstride_;
		(*tp).cstride_   = cstride_;
		(*tp).ptr_       = ptr_;
		(*tp).data_      = data_;
		(*tp).refs_      = refs_;
		(*tp).data_size_ = data_size_;
		if (data_) ++(*refs_);
		return tmp;
	}

	friend class Vector<T>;
};


// Constructors and destructor.
template <class T>
inline Matrix<T>::Matrix() : nn_(0), mm_(0), rstride_(0), cstride_(1), ptr_(0), data_(0),
	refs_(0), data_size_(0) {;}

template <class T>
Matrix<T>::Matrix(size_t n, size_t m) : nn_(n), mm_(m), rstride_(m), cstride_(1), ptr_(0),
	data_(0), refs_(0),	data_size_(n * m)
{
	if (data_size_)
	{
		size_t i;
		data_ = alloc_.allocate(data_size_, 0);
		if (!data_) vmerror("Memory allocation error in Matrix constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		T zero = Zero<T>();
		T *p1 = ptr_;
		for (i = 0; i < data_size_; ++i)
		{
			alloc_.construct(p1, zero);
			++p1;
		}
	}
}

template <class T>
Matrix<T>::Matrix(const NoInit_& noinit, size_t n, size_t m) : nn_(n), mm_(m), rstride_(m),
	cstride_(1), ptr_(0), data_(0), refs_(0), data_size_(n * m)
{
	if (data_size_)
	{
		size_t i;
		data_ = alloc_.allocate(data_size_, 0);
		if (!data_) vmerror("Memory allocation error in Matrix constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		if (!VMTraits<T>::is_simple)
		{
			T *p1 = ptr_;
			for (i = 0; i < data_size_; ++i)
			{
				new ((void *)p1) T();
				++p1;
			}
		}
	}
}

template <class T>
Matrix<T>::Matrix(const T& a, size_t n, size_t m) : nn_(n), mm_(m),	rstride_(m), cstride_(1),
	ptr_(0), data_(0), refs_(0), data_size_(n * m)
{
	if (data_size_)
	{
		size_t i;
		data_ = alloc_.allocate(data_size_, 0);
		if (!data_) vmerror("Memory allocation error in Matrix constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		for (i = 0; i < data_size_; ++i)
		{
			alloc_.construct(p1, a);
			++p1;
		}
	}
}

template <class T>
Matrix<T>::Matrix(const T* a, size_t n, size_t m, ptrdiff_t s) : nn_(n), mm_(m), rstride_(m),
	cstride_(1), ptr_(0), data_(0), refs_(0), data_size_(n * m)
{
	if (data_size_)
	{
		size_t i;
		data_ = alloc_.allocate(data_size_, 0);
		if (!data_) vmerror("Memory allocation error in Matrix constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		const T *p2 = a;
		for (i = 0; i < data_size_; ++i)
		{
			alloc_.construct(p1, *p2);
			++p1;
			p2 += s;
		}
	}
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& rhs) : nn_(rhs.nn_), mm_(rhs.mm_), rstride_(rhs.rstride_),
	cstride_(rhs.cstride_), ptr_(rhs.ptr_), data_(rhs.data_), refs_(rhs.refs_),
	data_size_(rhs.data_size_)
{
	if (data_) ++(*refs_);
}

template <class T>
Matrix<T>::Matrix(const Matrix<real_type>& re, const Matrix<real_type>& im) :
	nn_(re.nrows()), mm_(re.ncols()), rstride_(re.ncols()), cstride_(1), ptr_(0),
	data_(0), refs_(0), data_size_(re.nrows() * re.ncols())
{
	if (im.cols() != mm_ || im.rows() != nn_)
		vmerror("Size mismatch error in complex Matrix constructor.");
	if (data_size_)
	{
		size_t i, j;
		data_ = alloc_.allocate(data_size_);
		if (!data_) vmerror("Memory allocation error in Matrix constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		const real_type *p2, *p3;
		ptrdiff_t rs2 = re.rstride(), cs2 = re.cstride();
		ptrdiff_t rs3 = im.rstride(), cs3 = im.cstride();
		for (i = 0; i < nn_; ++i)
		{
			p2 = re.data() + i * rs2;
			p3 = im.data() + i * rs3;
			for (j = 0; j < mm_; ++j)
			{
				alloc_.construct(p1, T(*p2, *p3));
				++p1;
				p2 += cs2;
				p3 += cs3;
			}
		}
	}
}

template <class T>
Matrix<T>::Matrix(const Matrix<real_type>& re, const real_type& im) :
	nn_(re.nrows()), mm_(re.ncols()), rstride_(re.ncols()), cstride_(1), ptr_(0),
	data_(0), refs_(0), data_size_(re.nrows() * re.ncols())
{
	if (data_size_)
	{
		size_t i, j;
		data_ = alloc_.allocate(data_size_, 0);
		if (!data_) vmerror("Memory allocation error in Matrix constructor.");
		ptr_ = data_;
		refs_ = new size_t(1);
		T *p1 = ptr_;
		const real_type *p2;
		ptrdiff_t rs2 = re.rstride(), cs2 = re.cstride();
		for (i = 0; i < nn_; ++i)
		{
			p2 = re.data() + i * rs2;
			for (j = 0; j < mm_; ++j)
			{
				alloc_.construct(p1, T(*p2, im));
				++p1;
				p2 += cs2;
			}
		}
	}
}

// Destructor.
template <class T>
inline Matrix<T>::~Matrix()
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
}


// Overloaded operators.
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs)
{
	if (nn_ == rhs.nn_ && mm_ == rhs.mm_)
	{
		if (this == &rhs) return *this;
		if (data_ == rhs.data_) return operator=(rhs.copy());
		size_t i, j;
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) = (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator+=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) += (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator-=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) -= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator*=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) *= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator/=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) /= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator%=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator%=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) %= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator&=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator&=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) &= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator^=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator^=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) ^= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator|=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator|=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) |= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator<<=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator<<=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) <<= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator>>=(const Matrix<T>& rhs)
{
	size_t i, j;
	size_t n = rhs.nn_;
	size_t m = rhs.mm_;
	if (n == nn_ && m == mm_)
	{
		if (data_ == rhs.data_) return operator>>=(rhs.copy());
		T *p1, *p2;
		for (i = 0; i < nn_; ++i)
		{
			p1 = ptr_ + i * rstride_;
			p2 = rhs.ptr_ + i * rhs.rstride_;
			for (j = 0; j < mm_; ++j)
			{
				(*p1) >>= (*p2);
				p1 += cstride_;
				p2 += rhs.cstride_;
			}
		}
	}
	else
		vmerror("Size mismatch error in Matrix assignment operator.");
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const T& a)
{
	size_t i, j;
	for (i = 0; i < nn_; ++i)
	T *p1;
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) = a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) += a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) -= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) *= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) /= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator%=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) %= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator&=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) &= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator^=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) ^= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator|=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) |= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator>>=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) >>= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator<<=(const T& a)
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) <<= a;
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
inline Matrix<T> Matrix<T>::operator+() const
{
	Matrix<T> result = copy();
	size_t i, nm = nn_ * mm_;
	T *p1 = result.ptr_;
	for (i = 0; i < nm; ++i)
	{
		(*p1) = +(*p1);
		++p1;
	}
	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator-() const
{
	Matrix<T> result = copy();
	size_t i, nm = nn_ * mm_;
	T *p1 = result.ptr_;
	for (i = 0; i < nm; ++i)
	{
		(*p1) = -(*p1);
		++p1;
	}
	return result;
}

template <class T>
Matrix<bool> Matrix<T>::operator!() const
{
	Matrix<bool> result(nn_, mm_);
	size_t i, j;
	bool *p1 = result.data();
	T    *p2; 
	for (i = 0; i < nn_; ++i)
	{
		p2 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) = !(*p2);
			++p1;
			p2 += cstride_;
		}
	}
	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator~() const
{
	Matrix<T> result = copy();
	size_t i, nm = nn_ * mm_;
	T *p1 = result.ptr_;
	for (i = 0; i < nm; ++i)
	{
			(*p1) = ~(*p1);
			++p1;
	}
	return result;
}

template <class T>
inline Vector<T>::iterator Matrix<T>::operator[](size_t i)
{
#ifdef VM_BOUNDS_CHECK
	if (i >= nn_)
		vmerror("Index out of bounds error.");
#endif // VM_BOUNDS_CHECK
	return Vector<T>::iterator(ptr_ + i * rstride_, cstride_);
}

template <class T>
inline Vector<T>::const_iterator Matrix<T>::operator[](size_t i) const
{
#ifdef VM_BOUNDS_CHECK
	if (i >= nn_)
		vmerror("Index out of bounds error.");
#endif // VM_BOUNDS_CHECK
	return Vector<T>::const_iterator(ptr_ + i * rstride_, cstride_);
}

template <class T>
inline T& Matrix<T>::operator()(size_t i, size_t j)
{
#ifdef VM_BOUNDS_CHECK
	if (i >= nn_ || j >= mm_)
		vmerror("Index out of bounds error.");
#endif // VM_BOUNDS_CHECK
	return ptr_[i * rstride_ + j * cstride_];
}

template <class T>
inline const T& Matrix<T>::operator()(size_t i, size_t j) const
{
#ifdef VM_BOUNDS_CHECK
	if (i >= nn_ || j >= mm_)
		vmerror("Index out of bounds error.");
#endif // VM_BOUNDS_CHECK
	return ptr_[i * rstride_ + j * cstride_];
}


// Member functions.
template <class T>
Matrix<T>& Matrix<T>::apply(T (*fn)(T x))
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) = fn(*p1);
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
Matrix<T>& Matrix<T>::apply(T (*fn)(const T& x))
{
	size_t i, j;
	T *p1;
	for (i = 0; i < nn_; ++i)
	{
		p1 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) = fn(*p1);
			p1 += cstride_;
		}
	}
	return *this;
}

template <class T>
inline bool Matrix<T>::isdeep() const
{
	if (!data_) return true;
	if (*refs_ == 1 && rstride_ == mm_	&& cstride_ == 1 &&
		data_size_ == (nn_ * mm_))
		return true;
	else
		return false;
}

template <class T>
inline void Matrix<T>::makedeep()
{
	if (!isdeep()) reference(copy());
}

template <class T>
inline void Matrix<T>::reshape(size_t n, size_t m)
{
	free();
	if (n && m)
	{
		Matrix<T> tmp(n, m);
		reference(tmp);
	}
}

template <class T>
inline void Matrix<T>::free()
{
	Matrix<T> tmp;
	reference(tmp);
}

template <class T>
void Matrix<T>::resize(size_t n, size_t m)
{
	if (n && m)
	{
		size_t nm = n * m;
		T *new_data = alloc_.allocate(nm, 0);
		if (!new_data) vmerror("Memory allocation error in Matrix::resize().");
		T *p1 = new_data, *p2;
		ptrdiff_t rs = rstride_, cs = cstride_;
		size_t i, j;
		size_t mn = nn_ < n ? nn_ : n;
		size_t mm = mm_ < m ? mm_ : m;
		T zero = Zero<T>();
		for (i = 0; i < mn; ++i)
		{
			p2 = ptr_ + i * rs;
			for (j = 0; j < mm; ++j)
			{
				alloc_.construct(p1, (*p2));
				++p1;
				p2 += cs;
			}
			for (j = mn; j < m; ++j)
			{
				alloc_.construct(p1, zero);
				++p1;
			}
		}
		for (i = mn; i < n; ++i)
		{
			for (j = 0; j < m; ++j)
			{
				alloc_.construct(p1, zero);
				++p1;
			}
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
		mm_        = m;
		rstride_   = m;
		cstride_   = 1;
		data_size_ = nm;
	}
	else free();
}

template <class T>
Matrix<T> Matrix<T>::copy() const
{
	Matrix<T> tmp;
	if (data_)
	{
		alloc_type a2;
		T *new_data = a2.allocate(nn_ * mm_, 0);
		if (!new_data) vmerror("Memory allocation error in Matrix::copy().");
		T *p1 = new_data;
		const T *p2;
		ptrdiff_t rs = rstride_, cs = cstride_;
		size_t i, j;
		for (i = 0; i < nn_; ++i)
		{
			p2 = ptr_ + i * rs;
			for (j = 0; j < mm_; ++j)
			{
				a2.construct(p1, (*p2));
				++p1;
				p2 += cs;
			}
		}
		tmp.data_      = new_data;
		tmp.ptr_       = new_data;
		tmp.nn_        = nn_;
		tmp.mm_        = mm_;
		tmp.rstride_   = mm_;
		tmp.cstride_   = 1;
		tmp.data_size_ = nn_ * mm_;
		tmp.refs_ = new size_t(1);
	}
	return tmp;
}

template <class T>
void Matrix<T>::reference(const Matrix<T>& rhs)
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
		mm_        = rhs.mm_;
		rstride_   = rhs.rstride_;
		cstride_   = rhs.cstride_;
		ptr_       = rhs.ptr_;
		data_size_ = rhs.data_size_;
		if (data_) ++(*refs_);
	}
}

template <class T>
inline Matrix<T> Matrix<T>::submatrix(size_t rb, size_t cb, size_t rn, size_t cn,
	ptrdiff_t rs, ptrdiff_t cs) const
{
	if (rb >= nn_ ||
		cb >= mm_ ||
		rb + (rn - 1) * rs < 0 ||
		rb + (rn - 1) * rs >= nn_ ||
		cb + (cn - 1) * cs < 0 ||
		cb + (cn - 1) * cs >= mm_)
		vmerror("Index out of bounds error.");
	Matrix<T> tmp = *this;
	tmp.nn_      = rn;
	tmp.mm_      = cn;
	tmp.rstride_ = rstride_ * rs;
	tmp.cstride_ = cstride_ * cs;
	tmp.ptr_     = ptr_ + rstride_ * rb + cstride_ * cb;
	return tmp;
}

template <class T>
inline Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T> tmp = *this;
	tmp.nn_      = mm_;
	tmp.mm_      = nn_;
	tmp.rstride_ = cstride_;
	tmp.cstride_ = rstride_;
	return tmp;
}

template <class T>
Matrix<VMTraits<T>::real_type> Matrix<T>::real() const
{
	Matrix<real_type> tmp(nn_, mm_);
	size_t i, j;
	real_type *p1 = tmp.data();
	T  *p2;
	for (i = 0; i < nn_; ++i)
	{
		p2 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) = (*p2).real();
			++p1;
			p2 += cstride_;
		}
	}
	return tmp;
}

template <class T>
Matrix<VMTraits<T>::real_type> Matrix<T>::imag() const
{
	Matrix<real_type> tmp(nn_, mm_);
	size_t i, j;
	real_type *p1 = tmp.data();
	T  *p2;
	for (i = 0; i < nn_; ++i)
	{
		p2 = ptr_ + i * rstride_;
		for (j = 0; j < mm_; ++j)
		{
			(*p1) = (*p2).imag();
			++p1;
			p2 += cstride_;
		}
	}
	return tmp;
}

template <class T>
inline size_t Matrix<T>::nrows() const
{
	return nn_;
}

template <class T>
inline size_t Matrix<T>::ncols() const
{
	return mm_;
}

template <class T>
inline T* Matrix<T>::data()
{
	return ptr_;
}

template <class T>
inline const T* Matrix<T>::data() const
{
	return ptr_;
}

template <class T>
inline ptrdiff_t Matrix<T>::rstride() const
{
	return rstride_;
}

template <class T>
inline ptrdiff_t Matrix<T>::cstride() const
{
	return cstride_;
}

template <class T>
inline Vector<T> Matrix<T>::row(size_t i) const
{
	if (i >= nn_) vmerror("Index out of bounds error.");
	Vector<T> tmp;
	tmp.data_      = data_;
	tmp.refs_      = refs_;
	tmp.nn_        = mm_;
	tmp.stride_    = cstride_;
	tmp.ptr_       = ptr_ + i * rstride_;
	tmp.data_size_ = data_size_;
	if (data_) ++(*refs_);
	return tmp;
}

template <class T>
inline Vector<T> Matrix<T>::column(size_t j) const
{
	if (j >= mm_) vmerror("Index out of bounds error.");
	Vector<T> tmp;
	tmp.data_      = data_;
	tmp.refs_      = refs_;
	tmp.nn_        = nn_;
	tmp.stride_    = rstride_;
	tmp.ptr_       = ptr_ + j * cstride_;
	tmp.data_size_ = data_size_;
	if (data_) ++(*refs_);
	return tmp;
}

template <class T>
inline Vector<T> Matrix<T>::unwrap() const
{
	if (rstride_ != cstride_ * mm_)
		vmerror("Size mismatch error in Matrix::unwrap().");
	Vector<T> tmp;
	tmp.data_      = data_;
	tmp.refs_      = refs_;
	tmp.nn_        = nn_ * mm_;
	tmp.stride_    = cstride_;
	tmp.ptr_       = ptr_;
	tmp.data_size_ = data_size_;
	if (data_) ++(*refs_);
	return tmp;
}

template <class T>
inline Vector<T>::iterator Matrix<T>::rbegin(size_t n)
{
	if (n >= nn_) vmerror("Index out of bounds error.");
	return Vector<T>::iterator(ptr_ + n * rstride_, cstride_);
}

template <class T>
inline Vector<T>::const_iterator Matrix<T>::rbegin(size_t n) const
{
	if (n >= nn_) vmerror("Index out of bounds error.");
	return Vector<T>::const_iterator(ptr_ + n * rstride_, cstride_);
}

template <class T>
inline Vector<T>::iterator Matrix<T>::rend(size_t n)
{
	if (n >= nn_) vmerror("Index out of bounds error.");
	return Vector<T>::iterator(ptr_ + n * rstride_ + mm_ * cstride_, cstride_);
}

template <class T>
inline Vector<T>::const_iterator Matrix<T>::rend(size_t n) const
{
	if (n >= nn_) vmerror("Index out of bounds error.");
	return Vector<T>::const_iterator(ptr_ + n * rstride_ + mm_ * cstride_, cstride_);
}

template <class T>
inline Vector<T>::iterator Matrix<T>::cbegin(size_t n)
{
	if (n >= mm_) vmerror("Index out of bounds error.");
	return Vector<T>::iterator(ptr_ + n * cstride_, rstride_);
}

template <class T>
inline Vector<T>::const_iterator Matrix<T>::cbegin(size_t n) const
{
	if (n >= mm_) vmerror("Index out of bounds error.");
	return Vector<T>::const_iterator(ptr_ + n * cstride_, rstride_);
}

template <class T>
inline Vector<T>::iterator Matrix<T>::cend(size_t n)
{
	if (n >= mm_) vmerror("Index out of bounds error.");
	return Vector<T>::iterator(ptr_ + n * cstride_ + nn_ * rstride_, rstride_);
}

template <class T>
inline Vector<T>::const_iterator Matrix<T>::cend(size_t n) const
{
	if (n >= mm_) vmerror("Index out of bounds error.");
	return Vector<T>::const_iterator(ptr_ + n * cstride_ + nn_ * rstride_, rstride_);
}

template <class T>
void Matrix<T>::write(std::ostream& out) const
{
	size_t i, j;
	for (i = 0; i < nn_; ++i)
	{
		for (j = 0; j < mm_; ++j)
		{
			out << ptr_[i * rstride_ + j * cstride_];
			if (j == mm_ - 1)
			{
				if (i < nn_ - 1) out << '\n';
			}
			else out << '\t';
		}
	}
}

template <class T>
void Matrix<T>::read(std::istream& in)
{
	Vector<T> tmp;
	tmp.read(in);
	if (tmp.size()) reference(tmp.rowmat());
	else free();
	return;
}

template <class T>
void Matrix<T>::read(std::istream& in, size_t ncount)
{
	Vector<T> tmp;
	tmp.read(in, ncount);
	if (tmp.size()) reference(tmp.rowmat());
	else free();
	return;
}

#endif // MATRIX_H
