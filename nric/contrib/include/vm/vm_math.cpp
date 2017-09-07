#ifndef VM_MATH_CPP
#define VM_MATH_CPP

/******************************************************************************
 *
 *	vm_math.cpp
 *	
 *	Source file containing math operators and functions for Vectors and
 *	Matrices.
 *
 *	This source file is included by the header file vm_math.h, and it is
 *	dependent on the following header files:
 *
 *	vec_mat.h  vm_traits.h  vector.h  matrix.h  vm_mat.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

// Binary arithmetic operators for vectors. Vectors must be the same length.
template <class T>
Vector<T> operator+(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result += rhs;
	return result;
}

template <class T>
Vector<T> operator-(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result -= rhs;
	return result;
}

template <class T>
Vector<T> operator*(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result *= rhs;
	return result;
}

template <class T>
Vector<T> operator/(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result /= rhs;
	return result;
}

template <class T>
Vector<T> operator%(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result %= rhs;
	return result;
}

template <class T>
Vector<T> operator&(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result &= rhs;
	return result;
}

template <class T>
Vector<T> operator^(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result ^= rhs;
	return result;
}

template <class T>
Vector<T> operator|(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result |= rhs;
	return result;
}

template <class T>
Vector<T> operator<<(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result <<= rhs;
	return result;
}

template <class T>
Vector<T> operator>>(const Vector<T>& lhs, const Vector<T>& rhs)
{
	Vector<T> result = lhs.copy();
	result >>= rhs;
	return result;
}

// Binary relational operators for vectors. Vectors must be the same size.
template <class T>
Vector<bool> operator==(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] == rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator!=(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] != rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator<(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] < rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator>(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] > rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator<=(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] <= rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator>=(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] >= rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator&&(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] && rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator||(const Vector<T>& lhs, const Vector<T>& rhs)
{
	if (lhs.size() != rhs.size())
		vmerror("Size mismatch error in relational operator.");
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] || rhs[i]);
	return result;
}


// Binary arithmetic operators for vectors and scalers.
template <class T>
Vector<T> operator+(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result += a;
	return result;
}

template <class T>
Vector<T> operator+(const T& a, const Vector<T>& rhs)
{
	Vector<T> result = rhs.copy();
	result += a;
	return result;
}

template <class T>
Vector<T> operator-(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result -= a;
	return result;
}

template <class T>
Vector<T> operator-(const T& a, const Vector<T>& rhs)
{
	Vector<T> result(a, rhs.size());
	result -= rhs;
	return result;
}

template <class T>
Vector<T> operator*(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result *= a;
	return result;
}

template <class T>
Vector<T> operator*(const T& a, const Vector<T>& rhs)
{
	Vector<T> result = rhs.copy();
	result *= a;
	return result;
}

template <class T>
Vector<T> operator/(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result /= a;
	return result;
}

template <class T>
Vector<T> operator/(const T& a, const Vector<T>& rhs)
{
	Vector<T> result(a, rhs.size());
	result /= rhs;
	return result;
}

template <class T>
Vector<T> operator%(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result %= a;
	return result;
}

template <class T>
Vector<T> operator%(const T& a, const Vector<T>& rhs)
{
	Vector<T> result(a, rhs.size());
	result %= rhs;
	return result;
}

template <class T>
Vector<T> operator&(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result &= a;
	return result;
}

template <class T>
Vector<T> operator&(const T& a, const Vector<T>& rhs)
{
	Vector<T> result = rhs.copy();
	result &= a;
	return result;
}

template <class T>
Vector<T> operator^(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result ^= a;
	return result;
}

template <class T>
Vector<T> operator^(const T& a, const Vector<T>& rhs)
{
	Vector<T> result = rhs.copy();
	result ^= a;
	return result;
}

template <class T>
Vector<T> operator|(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result |= a;
	return result;
}

template <class T>
Vector<T> operator|(const T& a, const Vector<T>& rhs)
{
	Vector<T> result = rhs.copy();
	result |= a;
	return result;
}

template <class T>
Vector<T> operator<<(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result <<= a;
	return result;
}

template <class T>
Vector<T> operator<<(const T& a, const Vector<T>& rhs)
{
	Vector<T> result(a, rhs.size());
	result <<= rhs;
	return result;
}

template <class T>
Vector<T> operator>>(const Vector<T>& lhs, const T& a)
{
	Vector<T> result = lhs.copy();
	result >>= a;
	return result;
}

template <class T>
Vector<T> operator>>(const T& a, const Vector<T>& rhs)
{
	Vector<T> result(a, rhs.size());
	result >>= rhs;
	return result;
}

// Binary relational operators for vectors and scalers.
template <class T>
Vector<bool> operator==(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] == a);
	return result;
}

template <class T>
Vector<bool> operator==(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (rhs[i] == a);
	return result;
}

template <class T>
Vector<bool> operator!=(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] != a);
	return result;
}

template <class T>
Vector<bool> operator!=(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (rhs[i] != a);
	return result;
}

template <class T>
Vector<bool> operator<(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] < a);
	return result;
}

template <class T>
Vector<bool> operator<(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (a < rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator>(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] > a);
	return result;
}

template <class T>
Vector<bool> operator>(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (a > rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator<=(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] <= a);
	return result;
}

template <class T>
Vector<bool> operator<=(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (a <= rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator>=(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] >= a);
	return result;
}

template <class T>
Vector<bool> operator>=(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (a >= rhs[i]);
	return result;
}

template <class T>
Vector<bool> operator&&(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] && a);
	return result;
}

template <class T>
Vector<bool> operator&&(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (rhs[i] && a);
	return result;
}

template <class T>
Vector<bool> operator||(const Vector<T>& lhs, const T& a)
{
	size_t n = lhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (lhs[i] || a);
	return result;
}

template <class T>
Vector<bool> operator||(const T& a, const Vector<T>& rhs)
{
	size_t n = rhs.size(), i;
	Vector<bool> result(n);
	for (i = 0; i < n; ++i)
		result[i] = (rhs[i] || a);
	return result;
}


// Binary arithmetic operators for matrices. Matrices must be the same dimension.
template <class T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result += rhs;
	return result;
}

template <class T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result -= rhs;
	return result;
}

template <class T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result *= rhs;
	return result;
}

template <class T>
Matrix<T> operator%(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result %= rhs;
	return result;
}

template <class T>
Matrix<T> operator&(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result &= rhs;
	return result;
}

template <class T>
Matrix<T> operator^(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result ^= rhs;
	return result;
}

template <class T>
Matrix<T> operator|(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result |= rhs;
	return result;
}

template <class T>
Matrix<T> operator<<(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result <<= rhs;
	return result;
}

template <class T>
Matrix<T> operator>>(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result >>= rhs;
	return result;
}

template <class T>
Matrix<T> operator/(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix<T> result = lhs.copy();
	result /= rhs;
	return result;
}

// Binary relational operators for matrices. Matrices must have the same dimension.
template <class T>
Matrix<bool> operator==(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) == rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) != rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator<(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) < rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator>(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) > rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator<=(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) <= rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator>=(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) >= rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator&&(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) && rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator||(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	if (n != rhs.nrows() || m != rhs.ncols())
		vmerror("Size mismatch error in relational operator.");
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) || rhs(i, j));
	}
	return result;
}

// Binary arithmetic operators for matrices with scalers.
template <class T>
Matrix<T> operator+(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result += a;
	return result;
}

template <class T>
Matrix<T> operator+(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result = rhs.copy();
	result += a;
	return result;
}

template <class T>
Matrix<T> operator-(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result -= a;
	return result;
}

template <class T>
Matrix<T> operator-(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result(a, rhs.nrows() * rhs.ncols());
	result -= rhs;
	return result;
}

template <class T>
Matrix<T> operator*(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result *= a;
	return result;
}

template <class T>
Matrix<T> operator*(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result = rhs.copy();
	result *= a;
	return result;
}

template <class T>
Matrix<T> operator/(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result /= a;
	return result;
}

template <class T>
Matrix<T> operator/(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result(a, rhs.nrows() * rhs.ncols());
	result /= rhs;
	return result;
}

template <class T>
Matrix<T> operator%(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result %= a;
	return result;
}

template <class T>
Matrix<T> operator%(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result(a, rhs.nrows() * rhs.ncols());
	result %= rhs;
	return result;
}

template <class T>
Matrix<T> operator&(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result &= a;
	return result;
}

template <class T>
Matrix<T> operator&(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result = rhs.copy();
	result &= a;
	return result;
}

template <class T>
Matrix<T> operator^(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result ^= a;
	return result;
}

template <class T>
Matrix<T> operator^(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result = rhs.copy();
	result ^= a;
	return result;
}

template <class T>
Matrix<T> operator|(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result |= a;
	return result;
}

template <class T>
Matrix<T> operator|(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result = rhs.copy();
	result |= a;
	return result;
}

template <class T>
Matrix<T> operator<<(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result <<= a;
	return result;
}

template <class T>
Matrix<T> operator<<(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result(a, rhs.nrows() * rhs.ncols());
	result <<= rhs;
	return result;
}

template <class T>
Matrix<T> operator>>(const Matrix<T>& lhs, const T& a)
{
	Matrix<T> result = lhs.copy();
	result >>= a;
	return result;
}

template <class T>
Matrix<T> operator>>(const T& a, const Matrix<T>& rhs)
{
	Matrix<T> result(a, rhs.nrows() * rhs.ncols());
	result >>= rhs;
	return result;
}

// Binary relational operators for matrices and scalers.
template <class T>
Matrix<bool> operator==(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) == a);
	}
	return result;
}

template <class T>
Matrix<bool> operator==(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a == rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator!=(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) != a);
	}
	return result;
}

template <class T>
Matrix<bool> operator!=(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a != rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator<(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) < a);
	}
	return result;
}

template <class T>
Matrix<bool> operator<(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a < rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator>(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) > a);
	}
	return result;
}

template <class T>
Matrix<bool> operator>(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a > rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator<=(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) <= a);
	}
	return result;
}

template <class T>
Matrix<bool> operator<=(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a <= rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator>=(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) >= a);
	}
	return result;
}

template <class T>
Matrix<bool> operator>=(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a >= rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator&&(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) && a);
	}
	return result;
}

template <class T>
Matrix<bool> operator&&(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a && rhs(i, j));
	}
	return result;
}

template <class T>
Matrix<bool> operator||(const Matrix<T>& lhs, const T& a)
{
	size_t n = lhs.nrows(), m = lhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (lhs(i, j) || a);
	}
	return result;
}

template <class T>
Matrix<bool> operator||(const T& a, const Matrix<T>& rhs)
{
	size_t n = rhs.nrows(), m = rhs.ncols(), i, j;
	Matrix<bool> result(n, m);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			result(i, j) = (a || rhs(i, j));
	}
	return result;
}


// Overloaded math functions for vectors.
// Complex specific.
template <class T>
Vector<T> norm(const Vector<std::complex<T> >& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = norm(x[i]);
	return tmp;
}

template <class T>
Vector<T> arg(const Vector<std::complex<T> >& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = arg(x[i]);
	return tmp;
}

template <class T>
Vector<std::complex<T> > polar(const Vector<T>& rho, const Vector<T>& theta)
{
	using namespace std;
	if (rho.size() != theta.size())
		vmerror("Size mismatch error in polar function.");
	size_t n = rho.size();
	Vector<complex<T> > tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = polar(rho[i], theta[i]);
	return tmp;
}

template <class T>
Vector<std::complex<T> > polar(const Vector<T>& rho, const T& theta)
{
	using namespace std;
	size_t n = rho.size();
	Vector<complex<T> > tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = polar(rho[i], theta);
	return tmp;
}

template <class T>
Vector<std::complex<T> > polar(const T& rho, const Vector<T>& theta)
{
	using namespace std;
	size_t n = theta.size();
	Vector<complex<T> > tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = polar(rho, theta[i]);
	return tmp;
}

template <class T>
inline Vector<T> real(const Vector<std::complex<T> >& x)
{
	return x.real();
}

template <class T>
inline Vector<T> imag(const Vector<std::complex<T> >& x)
{
	return x.imag();
}

template <class T>
Vector<std::complex<T> > conj(const Vector<std::complex<T> >& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<complex<T> > tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = conj(x[i]);
	return tmp;
}

// General.
template<class T>
Vector<VMTraits<T>::real_type> abs(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<VMTraits<T>::real_type> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = abs(x[i]);
	return tmp;
}

template<class T>
Vector<T> acos(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = acos(x[i]);
	return tmp;
}

template<class T>
Vector<T> asin(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = asin(x[i]);
	return tmp;
}

template<class T>
Vector<T> atan(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = atan(x[i]);
	return tmp;
}

template<class T>
Vector<T> atan2(const Vector<T>& x, const Vector<T>& y)
{
	using namespace std;
	if (x.size() != y.size())
		vmerror("Size mismatch error in atan2 function.");
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = atan2(x[i], y[i]);
	return tmp;
}

template<class T>
Vector<T> atan2(const Vector<T> x, const T& y)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = atan2(x[i], y);
	return tmp;
}

template<class T>
Vector<T> atan2(const T& x, const Vector<T>& y)
{
	using namespace std;
	size_t n = y.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = atan2(x, y[i]);
	return tmp;
}

template<class T>
Vector<T> cos(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = cos(x[i]);
	return tmp;
}

template<class T>
Vector<T> cosh(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = cosh(x[i]);
	return tmp;
}

template<class T>
Vector<T> exp(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = exp(x[i]);
	return tmp;
}

template<class T>
Vector<T> log(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = log(x[i]);
	return tmp;
}

template<class T>
Vector<T> log10(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = log10(x[i]);
	return tmp;
}

template<class T>
Vector<T> pow(const Vector<T>& x, const Vector<T>& y)
{
	using namespace std;
	if (x.size() != y.size())
		vmerror("Size mismatch error in pow function.");
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = pow(x[i], y[i]);
	return tmp;
}

template<class T>
Vector<T> pow(const Vector<T> x, const T& y)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = pow(x[i], y);
	return tmp;
}

template<class T>
Vector<T> pow(const T& x, const Vector<T>& y)
{
	using namespace std;
	size_t n = y.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = pow(x, y[i]);
	return tmp;
}

template<class T>
Vector<T> sin(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = sin(x[i]);
	return tmp;
}

template<class T>
Vector<T> sinh(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = sinh(x[i]);
	return tmp;
}

template<class T>
Vector<T> sqrt(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = sqrt(x[i]);
	return tmp;
}

template<class T>
Vector<T> tan(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = tan(x[i]);
	return tmp;
}

template<class T>
Vector<T> tanh(const Vector<T>& x)
{
	using namespace std;
	size_t n = x.size();
	Vector<T> tmp(n);
	for (size_t i = 0; i < n; ++i)
		tmp[i] = tanh(x[i]);
	return tmp;
}


// Overloaded math functions for matrices.
// Complex specific.
template <class T>
Matrix<T> norm(const Matrix<std::complex<T> >& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = norm(x(i, j));
	}
	return tmp;
}

template <class T>
Matrix<T> arg(const Matrix<std::complex<T> >& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = arg(x(i, j));
	}
	return tmp;
}

template <class T>
Matrix<std::complex<T> > polar(const Matrix<T>& rho, const Matrix<T>& theta)
{
	using namespace std;
	if (rho.nrows() != theta.nrows() || rho.ncols() != theta.ncols())
		vmerror("Size mismatch error in polar function.");
	Matrix<complex<T> > tmp(rho.nrows(), rho.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = polar(rho(i, j), theta(i, j));
	}
	return tmp;
}

template <class T>
Matrix<std::complex<T> > polar(const Matrix<T>& rho, const T& theta)
{
	using namespace std;
	Matrix<complex<T> > tmp(rho.nrows(), rho.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = polar(rho(i, j), theta);
	}
	return tmp;
}

template <class T>
Matrix<std::complex<T> > polar(const T& rho, const Matrix<T>& theta)
{
	using namespace std;
	Matrix<complex<T> > tmp(theta.nrows(), theta.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = polar(rho, theta(i, j));
	}
	return tmp;
}

template <class T>
inline Matrix<T> real(const Matrix<std::complex<T> >& x)
{
	return x.real();
}

template <class T>
inline Matrix<T> imag(const Matrix<std::complex<T> >& x)
{
	return x.imag();
}

template <class T>
Matrix<std::complex<T> > conj(const Matrix<std::complex<T> >& x)
{
	using namespace std;
	Matrix<complex<T> > tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = conj(x(i, j));
	}
	return tmp;
}

// General.
template<class T>
Matrix<VMTraits<T>::real_type> abs(const Matrix<T>& x)
{
	using namespace std;
	Matrix<VMTraits<T>::real_type> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = abs(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> acos(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = acos(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> asin(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = asin(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> atan(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = atan(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> atan2(const Matrix<T>& x, const Matrix<T>& y)
{
	using namespace std;
	if (x.nrows() != y.nrows() || x.ncols() != y.ncols())
		vmerror("Size mismatch in atan2 function.");
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = atan2(x(i, j), y(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> atan2(const Matrix<T> x, const T& y)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = atan2(x(i, j), y);
	}
	return tmp;
}

template<class T>
Matrix<T> atan2(const T& x, const Matrix<T>& y)
{
	using namespace std;
	Matrix<T> tmp(y.nrows(), y.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = atan2(x, y(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> cos(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = cos(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> cosh(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = cosh(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> exp(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = exp(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> log(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = log(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> log10(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = log10(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> pow(const Matrix<T>& x, const Matrix<T>& y)
{
	using namespace std;
	if (x.nrows() != y.nrows() || x.ncols() != y.ncols())
		vmerror("Size mismatch in pow function.");
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = pow(x(i, j), y(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> pow(const Matrix<T> x, const T& y)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = pow(x(i, j), y);
	}
	return tmp;
}

template<class T>
Matrix<T> pow(const T& x, const Matrix<T>& y)
{
	using namespace std;
	Matrix<T> tmp(y.nrows(), y.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = pow(x, y(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> sin(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = sin(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> sinh(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = sinh(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> sqrt(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = sqrt(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> tan(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = tan(x(i, j));
	}
	return tmp;
}

template<class T>
Matrix<T> tanh(const Matrix<T>& x)
{
	using namespace std;
	Matrix<T> tmp(x.nrows(), x.ncols());
	size_t n = tmp.nrows();
	size_t m = tmp.ncols();
	size_t i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
			tmp(i, j) = tanh(x(i, j));
	}
	return tmp;
}


#endif // VM_MATH_CPP
