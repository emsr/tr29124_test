#ifndef VM_TOOLS_CPP
#define VM_TOOLS_CPP

/******************************************************************************
 *
 *	vm_tools.cpp
 *	
 *	Source file containing tools for Vectors and Matrices.
 *
 *	This source file is included by the header file vm_tools.h, and it is
 *	dependent on the following header files:
 *
 *	vec_mat.h  vm_traits.h  vector.h  matrix.h  vm_tools.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/


// Statistics.
template <class T>
void Stat(const Vector<T>& x, double& m, double& var)
{
	size_t    len = x.size();
	T sum(0), sum2(0);
	Vector<T>::const_iterator p = x.begin();
	for (size_t i = 0; i < len; ++i)
	{
		sum  += (*p);
		sum2 += (*p) * (*p);
		++p;
	}
	m  = sum;
	m /= len;
	double tmp = sum2;
	tmp /= len;
	var = tmp - m * m;
}

// Sorting algorithms.
template <class T>
Vector<T> Sort(const Vector<T>& x)
{
	Vector<T> cx = x.copy();
	cx.sort();
	return cx;
}

template <class T>
Vector<int> Index(const Vector<T>& x)
{
	size_t n = x.size();
	Vector<int> indx(0, 1, n);
	IndirectComp<T> comp(x);
	int *first, *last;
	first = indx.data();
	last  = first + n;
	std::sort(first, last, comp);
	return indx;
}

// Matrix multiplication.
template <class T>
Matrix<T> MatMult(const Matrix<T>& x, const Matrix<T>& y)
{
	size_t xn = x.nrows();
	size_t xm = x.ncols();
	size_t yn = y.nrows();
	size_t ym = y.ncols();
	if (!xn || !ym || xm != yn) vmerror("Size mismatch error in MatMult() function.");
	Matrix<T> tmp(xn, ym);
	Vector<T>::const_iterator px, py;
	Vector<T>::iterator       pz;
	size_t i, j, k;
	for (i = 0; i < xn; ++i)
	{
		pz = tmp.rbegin(i);
		for (j = 0; j < ym; ++j)
		{
			px = x.rbegin(i);
			py = y.cbegin(j);
			(*pz) = (*px) * (*py);
			for (k = 1; k < xm; ++k)
			{
				++px;
				++py;
				(*pz) += (*px) * (*py);
			}
			++pz;
		}
	}
	return tmp;
}

template <class T>
Vector<T> MatMult(const Vector<T>& x, const Matrix<T>& y)
{
	size_t xm = x.size();
	size_t yn = y.nrows();
	size_t ym = y.ncols();
	if (!ym || xm != yn) vmerror("Size mismatch error in MatMult() function.");
	Vector<T> tmp(ym);
	Vector<T>::const_iterator px, py;
	Vector<T>::iterator       pz = tmp.begin();
	size_t i, k;
	for (i = 0; i < ym; ++i)
	{
		px = x.begin();
		py = y.cbegin(i);
		(*pz) = (*px) * (*py);
		for (k = 1; k < xm; ++k)
		{
			++px;
			++py;
			(*pz) += (*px) * (*py);
		}
		++pz;
	}
	return tmp;
}

template <class T>
Vector<T> MatMult(const Matrix<T>& x, const Vector<T>& y)
{
	size_t xn = x.nrows();
	size_t xm = x.ncols();
	size_t yn = y.size();
	if (!xn || xm != yn) vmerror("Size mismatch error in MatMult() function.");
	Vector<T> tmp(xn);
	Vector<T>::const_iterator px, py;
	Vector<T>::iterator       pz = tmp.begin();
	size_t i, k;
	for (i = 0; i < xn; ++i)
	{
		px = x.rbegin(i);
		py = y.begin();
		(*pz) = (*px) * (*py);
		for (k = 1; k < yn; ++k)
		{
			++px;
			++py;
			(*pz) += (*px) * (*py);
		}
		++pz;
	}
	return tmp;
}

template <class T>
Matrix<T> MatMult(const Vector<T>& x, const Vector<T>& y)
{
	size_t n = x.size();
	size_t m = y.size();
	if (!n || !m) vmerror("Size mismatch error in MatMult() function.");
	Matrix<T> tmp(n, m);
	Vector<T>::const_iterator px, py;
	Vector<T>::iterator		  pz;
	size_t i, j;
	px = x.begin();
	for (i = 0; i < n; ++i)
	{
		py = y.begin();
		pz = tmp.rbegin(i);
		for (j = 0; j < m; ++j)
		{
			(*pz) = (*px) * (*py);
			++py;
			++pz;
		}
		++px;
	}
	return tmp;
}

template <class T>
T DotProduct(const Vector<T>& x, const Vector<T>& y)
{
	size_t n = x.size();
	if (!n || y.size() != n) vmerror("Size mismatch error in DotProduct() function.");
	T tmp = x[0] * y[0];
	Vector<T>::const_iterator px = x.begin(), py = y.begin();
	for (size_t i = 1; i < n; ++i)
	{
		++px;
		++py;
		tmp += (*px) * (*py);
	}
	return tmp;
}


#endif // VM_TOOLS_CPP
