#ifndef _NR_TYPES_H_
#define _NR_TYPES_H_

// This is a modified version of nrtypes_nr.h for use with
// the VecMat software. Note that all typedefs included in
// the "vm_types.h" header have been removed from this header.

#include <complex>
#include <fstream>
#include "nrutil.h"
using namespace std;

typedef double DP;

// Vector Types
typedef const NRVec<bool> Vec_I_BOOL;
typedef NRVec<bool> Vec_O_BOOL, Vec_IO_BOOL;

typedef const NRVec<char> Vec_I_CHR;
typedef NRVec<char> Vec_O_CHR, Vec_IO_CHR;

typedef const NRVec<unsigned char> Vec_I_UCHR;
typedef NRVec<unsigned char> Vec_O_UCHR, Vec_IO_UCHR;

typedef const NRVec<int> Vec_I_INT;
typedef NRVec<int> Vec_O_INT, Vec_IO_INT;

typedef const NRVec<unsigned int> Vec_I_UINT;
typedef NRVec<unsigned int> Vec_O_UINT, Vec_IO_UINT;

typedef const NRVec<long> Vec_I_LNG;
typedef NRVec<long> Vec_O_LNG, Vec_IO_LNG;

typedef const NRVec<unsigned long> Vec_I_ULNG;
typedef NRVec<unsigned long> Vec_O_ULNG, Vec_IO_ULNG;

typedef const NRVec<float> Vec_I_SP;
typedef NRVec<float> Vec_O_SP, Vec_IO_SP;

typedef const NRVec<DP> Vec_I_DP;
typedef NRVec<DP> Vec_O_DP, Vec_IO_DP;

typedef const NRVec<complex<float> > Vec_I_CPLX_SP;
typedef NRVec<complex<float> > Vec_O_CPLX_SP, Vec_IO_CPLX_SP;

typedef const NRVec<complex<DP> > Vec_I_CPLX_DP;
typedef NRVec<complex<DP> > Vec_O_CPLX_DP, Vec_IO_CPLX_DP;

// Matrix Types

typedef const NRMat<bool> Mat_I_BOOL;
typedef NRMat<bool> Mat_O_BOOL, Mat_IO_BOOL;

typedef const NRMat<char> Mat_I_CHR;
typedef NRMat<char> Mat_O_CHR, Mat_IO_CHR;

typedef const NRMat<unsigned char> Mat_I_UCHR;
typedef NRMat<unsigned char> Mat_O_UCHR, Mat_IO_UCHR;

typedef const NRMat<int> Mat_I_INT;
typedef NRMat<int> Mat_O_INT, Mat_IO_INT;

typedef const NRMat<unsigned int> Mat_I_UINT;
typedef NRMat<unsigned int> Mat_O_UINT, Mat_IO_UINT;

typedef const NRMat<long> Mat_I_LNG;
typedef NRMat<long> Mat_O_LNG, Mat_IO_LNG;

typedef const NRVec<unsigned long> Mat_I_ULNG;
typedef NRMat<unsigned long> Mat_O_ULNG, Mat_IO_ULNG;

typedef const NRMat<float> Mat_I_SP;
typedef NRMat<float> Mat_O_SP, Mat_IO_SP;

typedef const NRMat<DP> Mat_I_DP;
typedef NRMat<DP> Mat_O_DP, Mat_IO_DP;

typedef const NRMat<complex<float> > Mat_I_CPLX_SP;
typedef NRMat<complex<float> > Mat_O_CPLX_SP, Mat_IO_CPLX_SP;

typedef const NRMat<complex<DP> > Mat_I_CPLX_DP;
typedef NRMat<complex<DP> > Mat_O_CPLX_DP, Mat_IO_CPLX_DP;


// 3D Matrix Types

typedef const NRMat3d<DP> Mat3D_I_DP;
typedef NRMat3d<DP> Mat3D_DP, Mat3D_O_DP, Mat3D_IO_DP;

// Miscellaneous Types

typedef NRVec<unsigned long *> Vec_ULNG_p;
typedef NRVec<NRMat<DP> *> Vec_Mat_DP_p;
typedef NRVec<fstream *> Vec_FSTREAM_p;

#endif /* _NR_TYPES_H_ */
