#ifndef VM_MATH_H
#define VM_MATH_H

/******************************************************************************
 *
 *	vm_math.h
 *	
 *	Header file for Vector and Matrix math operators and functions.
 *
 *	This header file is intended to be included by vec_mat.h, and it depends
 *	on the following other header files:
 *	
 *	vec_mat.h  vm_traits.h  vector.h  matrix.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

// Binary arithmetic operators for vectors. Vectors must be the same size.
template <class T>
Vector<T> operator+(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator-(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator*(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator/(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator%(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator&(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator^(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator|(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator<<(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<T> operator>>(const Vector<T>& lhs, const Vector<T>& rhs);

// Binary relational operators for vectors. Vectors must be the same size.
template <class T>
Vector<bool> operator==(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator!=(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator<(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator>(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator<=(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator>=(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator&&(const Vector<T>& lhs, const Vector<T>& rhs);
template <class T>
Vector<bool> operator||(const Vector<T>& lhs, const Vector<T>& rhs);

// Binary arithemtic operators for vectors and scalers.
template <class T>
Vector<T> operator+(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator+(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator-(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator-(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator*(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator*(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator/(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator/(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator%(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator%(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator&(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator&(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator^(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator^(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator|(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator|(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator<<(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator<<(const T& a, const Vector<T>& rhs);
template <class T>
Vector<T> operator>>(const Vector<T>& lhs, const T& a);
template <class T>
Vector<T> operator>>(const T& a, const Vector<T>& rhs);

// Binary relational operators for vectors and scalers.
template <class T>
Vector<bool> operator==(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator==(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator!=(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator!=(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator<(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator<(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator>(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator>(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator<=(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator<=(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator>=(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator>=(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator&&(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator&&(const T& a, const Vector<T>& rhs);
template <class T>
Vector<bool> operator||(const Vector<T>& lhs, const T& a);
template <class T>
Vector<bool> operator||(const T& a, const Vector<T>& rhs);


// Binary arithmetic operators for matrices. Matrices must have the same dimension.
template <class T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator/(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator%(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator&(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator^(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator|(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator<<(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator>>(const Matrix<T>& lhs, const Matrix<T>& rhs);

// Binary relational operators for matrices. Matrices must have the same dimension.
template <class T>
Matrix<bool> operator==(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator<(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator>(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator<=(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator>=(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator&&(const Matrix<T>& lhs, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator||(const Matrix<T>& lhs, const Matrix<T>& rhs);

// Binary arithmetic operators for matrices with scalers.
template <class T>
Matrix<T> operator+(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator+(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator-(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator-(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator*(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator*(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator/(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator/(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator%(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator%(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator&(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator&(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator^(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator^(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator|(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator|(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator<<(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator<<(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<T> operator>>(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<T> operator>>(const T& a, const Matrix<T>& rhs);

// Binary relational operators for matrices and scalers.
template <class T>
Matrix<bool> operator==(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator==(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator!=(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator!=(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator<(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator<(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator>(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator>(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator<=(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator<=(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator>=(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator>=(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator&&(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator&&(const T& a, const Matrix<T>& rhs);
template <class T>
Matrix<bool> operator||(const Matrix<T>& lhs, const T& a);
template <class T>
Matrix<bool> operator||(const T& a, const Matrix<T>& rhs);


// Overloaded math functions for vectors.
template <class T>
Vector<typename VMTraits<T>::real_type> abs(const Vector<T>& x);
template <class T>
Vector<T> acos(const Vector<T>& x);
template <class T>
Vector<T> asin(const Vector<T>& x);
template <class T>
Vector<T> atan(const Vector<T>& x);
template <class T>
Vector<T> atan2(const Vector<T>& x, const Vector<T>& y);
template <class T>
Vector<T> atan2(const Vector<T> x, const T& y);
template <class T>
Vector<T> atan2(const T& x, const Vector<T>& y);
template <class T>
Vector<T> cos(const Vector<T>& x);
template <class T>
Vector<T> cosh(const Vector<T>& x);
template <class T>
Vector<T> exp(const Vector<T>& x);
template <class T>
Vector<T> log(const Vector<T>& x);
template <class T>
Vector<T> log10(const Vector<T>& x);
template <class T>
Vector<T> pow(const Vector<T>& x, const Vector<T>& y);
template <class T>
Vector<T> pow(const Vector<T> x, const T& y);
template <class T>
Vector<T> pow(const T& x, const Vector<T>& y);
template <class T>
Vector<T> sin(const Vector<T>& x);
template <class T>
Vector<T> sinh(const Vector<T>& x);
template <class T>
Vector<T> sqrt(const Vector<T>& x);
template <class T>
Vector<T> tan(const Vector<T>& x);
template <class T>
Vector<T> tanh(const Vector<T>& x);

// Complex specific.
template <class T>
Vector<T> arg(const Vector<std::complex<T> >& x);
template <class T>
Vector<T> norm(const Vector<std::complex<T> >& x);
template <class T>
Vector<std::complex<T> > polar(const Vector<T>& rho, const Vector<T>& theta);
template <class T>
Vector<std::complex<T> > polar(const Vector<T>& rho, const T& theta);
template <class T>
Vector<std::complex<T> > polar(const T& rho, const Vector<T>& theta);
template <class T>
Vector<T> real(const Vector<std::complex<T> >& x);
template <class T>
Vector<T> imag(const Vector<std::complex<T> >& x);
template <class T>
Vector<std::complex<T> > conj(const Vector<std::complex<T> >& x);


// Overloaded math functions for matrices.
template <class T>
Matrix<typename VMTraits<T>::real_type> abs(const Matrix<T>& x);
template <class T>
Matrix<T> acos(const Matrix<T>& x);
template <class T>
Matrix<T> asin(const Matrix<T>& x);
template <class T>
Matrix<T> atan(const Matrix<T>& x);
template <class T>
Matrix<T> atan2(const Matrix<T>& x, const Matrix<T>& y);
template <class T>
Matrix<T> atan2(const Matrix<T> x, const T& y);
template <class T>
Matrix<T> atan2(const T& x, const Matrix<T>& y);
template <class T>
Matrix<T> cos(const Matrix<T>& x);
template <class T>
Matrix<T> cosh(const Matrix<T>& x);
template <class T>
Matrix<T> exp(const Matrix<T>& x);
template <class T>
Matrix<T> log(const Matrix<T>& x);
template <class T>
Matrix<T> log10(const Matrix<T>& x);
template <class T>
Matrix<T> pow(const Matrix<T>& x, const Matrix<T>& y);
template <class T>
Matrix<T> pow(const Matrix<T> x, const T& y);
template <class T>
Matrix<T> pow(const T& x, const Matrix<T>& y);
template <class T>
Matrix<T> sin(const Matrix<T>& x);
template <class T>
Matrix<T> sinh(const Matrix<T>& x);
template <class T>
Matrix<T> sqrt(const Matrix<T>& x);
template <class T>
Matrix<T> tan(const Matrix<T>& x);
template <class T>
Matrix<T> tanh(const Matrix<T>& x);

// Complex specific.
template <class T>
Matrix<T> arg(const Matrix<std::complex<T> >& x);
template <class T>
Matrix<T> norm(const Matrix<std::complex<T> >& x);
template <class T>
Matrix<std::complex<T> > polar(const Matrix<T>& rho, const Matrix<T>& theta);
template <class T>
Matrix<std::complex<T> > polar(const Matrix<T>& rho, const T& theta);
template <class T>
Matrix<std::complex<T> > polar(const T& rho, const Matrix<T>& theta);
template <class T>
Matrix<T> real(const Matrix<std::complex<T> >& x);
template <class T>
Matrix<T> imag(const Matrix<std::complex<T> >& x);
template <class T>
Matrix<std::complex<T> > conj(const Matrix<std::complex<T> >& x);


#include "vm/vm_math.cpp"

#endif // VM_MATH_H
