#ifndef VM_TRAITS_H
#define VM_TRAITS_H

/******************************************************************************
 *
 *	vm_traits.h
 *	
 *	Header file for numeric traits used by VecMat software.
 *
 *	This header file is intended to be included by vec_mat.h, and it depends
 *  on the following other header files:
 *
 *  vec_mat.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

/* This function provides a templated zero value for general types.
   An uninitialized static variable is set to zero if it is a standard
   type, and has its data members initialized with zero otherwise. The
   default constructor is then used for types with constructors, but any
   data members not initialized by the constructor will remain zero. This
   gives us the desired behavior for a zero function that can be used
   to initialize vectors and matrices.
*/

template <class T>
inline T Zero()
{
	static T zero;
	return zero;
}

/* Traits for various types.
   real_type should be the same as T, unless the type is
   complex, in which case it should be the real equivelent.
   This is used as the return type for the abs() function,
   and is also used for construction of complex objects from
   real ones.
   is_specialized should be true for all specialized types.
   is_simple should only be true for types which have no default
   constructor or destructor, or classes whose default constructor
   and destructor don't do anything.
*/

template <class T>
struct VMTraits
{
	typedef T real_type;
	enum {is_specialized = false};
	enum {is_simple = false};
};

// Traits for standard numeric types.
template <>
struct VMTraits<bool>
{
	typedef bool real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<char>
{
	typedef char real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<signed char>
{
	typedef signed char real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<unsigned char>
{
	typedef unsigned char real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<short>
{
	typedef short real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<unsigned short>
{
	typedef unsigned short real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<int>
{
	typedef int real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<unsigned int>
{
	typedef unsigned int real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<long>
{
	typedef long real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<unsigned long>
{
	typedef unsigned long real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<float>
{
	typedef float real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<double>
{
	typedef double real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};


// Traits for complex types.
template <>
struct VMTraits<std::complex<float> >
{
	typedef float real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};

template <>
struct VMTraits<std::complex<double> >
{
	typedef double real_type;
	enum {is_specialized = true};
	enum {is_simple = true};
};


#endif // VM_TRAITS_H

