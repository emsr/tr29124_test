#ifndef VEC_MAT_H
#define VEC_MAT_H 

/******************************************************************************
 *
 *	vec_mat.h
 *	
 *	Main header file for the VecMat software.
 *
 *	This header file should be included by any programs using this software.
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

// To remove bounds-checking in debug mode, comment out the following
// three lines.
#ifdef _DEBUG
#define VM_BOUNDS_CHECK
#endif // _DEBUG

// To include bounds-checking in release mode, uncomment the following
// line.
//#define VM_BOUNDS_CHECK

#include <cmath>
#include <sys/stat.h>
#include <string>
#include <complex>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <memory>

inline void vmerror(const std::string& error_text)
// VecMat error handler.
{
	std::cerr << "VecMat run-time error." << std::endl;
	std::cerr << error_text << std::endl;
	exit(1);
}

// Set the allocator to be used for memory allocation.
template <class T>
struct MemAlloc
{
	typedef std::allocator<T> alloc_type;
};

// NoInit object for use with constructors.
struct NoInit_ {};
static NoInit_ NoInit;
  
#include "vm/vm_traits.h"
#include "vm/vector.h"
#include "vm/matrix.h"
#include "vm/vm_types.h"
#include "vm/vm_math.h"
#include "vm/vm_tools.h"
#include "vm/vm_io.h"

#endif // VEC_MAT_H
