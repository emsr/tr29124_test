#ifndef VM_IO_H
#define VM_IO_H

/******************************************************************************
 *
 *	vm_io.h
 *	
 *	Header file for Vector and Matrix streaming IO operators.
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

// Stream operators for vectors.
template <class T>
std::ostream& operator<<(std::ostream& out, const Vector<T>& x);
template <class T>
std::istream& operator>>(std::istream& in, Vector<T>& x);

// Stream operators for matrices.
template <class T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& x);
template <class T>
std::istream& operator>>(std::istream& in, Matrix<T>& x);

// Binary file functions for vectors.
template <class T>
void ReadBinary(const std::string filename, Vector<T>& array, size_t skip);
template <class T>
void WriteBinary(const std::string filename, const Vector<T>& array);

// Binary file functions for matrices.
template <class T>
void ReadBinary(const std::string filename, Matrix<T>& array, size_t skip);
template <class T>
void WriteBinary(const std::string filename, const Matrix<T>& array);


#include "vm/vm_io.cpp"


#endif // VM_IO_H
