#ifndef VM_TYPES_H
#define VM_TYPES_H 

/******************************************************************************
 *
 *	vm_types.h
 *	
 *	Header file containing useful typedefs for vectors and matrices.
 *
 *	This source file is intended to be included by vec_mat.h, and it is
 *	dependent on the following header files:
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

// Complex types
typedef std::complex<float>		CPLX_SP;
typedef std::complex<double>	CPLX_DP;

// Vector Types
typedef Vector<bool>			Vec_BOOL;
typedef Vector<char>			Vec_CHR;
typedef Vector<signed char>		Vec_SCHR;
typedef Vector<unsigned char>	Vec_UCHR;
typedef Vector<short>			Vec_SHRT;
typedef Vector<unsigned short>	Vec_USHRT;
typedef Vector<int>				Vec_INT;
typedef Vector<unsigned int>	Vec_UINT;
typedef Vector<long>			Vec_LNG;
typedef Vector<unsigned long>	Vec_ULNG;
typedef Vector<float>			Vec_SP;
typedef Vector<double>			Vec_DP;
typedef Vector<CPLX_SP>			Vec_CPLX_SP;
typedef Vector<CPLX_DP>			Vec_CPLX_DP;

// Matrix Types
typedef Matrix<bool>			Mat_BOOL;
typedef Matrix<char>			Mat_CHR;
typedef Matrix<signed char>		Mat_SCHR;
typedef Matrix<unsigned char>	Mat_UCHR;
typedef Matrix<short>			Mat_SHRT;
typedef Matrix<unsigned short>	Mat_USHRT;
typedef Matrix<int>				Mat_INT;
typedef Matrix<unsigned int>	Mat_UINT;
typedef Matrix<long>			Mat_LNG;
typedef Matrix<unsigned long>	Mat_ULNG;
typedef Matrix<float>			Mat_SP;
typedef Matrix<double>			Mat_DP;
typedef Matrix<CPLX_SP>			Mat_CPLX_SP;
typedef Matrix<CPLX_DP>			Mat_CPLX_DP;


#endif // VM_TYPES_H
