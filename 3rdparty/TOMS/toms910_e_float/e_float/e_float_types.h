
#ifndef _E_FLOAT_TYPES_2004_02_28_H_
  #define _E_FLOAT_TYPES_2004_02_28_H_

  #if !defined(_MSC_VER) && !defined(__INTEL_COMPILER) && !defined(__GNUC__)
    #error "Compiler not supported: Types can not be determined"
  #endif

  #if defined(_MSC_VER) || defined(__INTEL_COMPILER)
  
    typedef   signed __int64  INT64;
    typedef unsigned __int64 UINT64;

  #elif defined(__GNUC__)

    typedef   signed long long  INT64;
    typedef unsigned long long UINT64;

  #endif

  typedef   signed int    INT32;
  typedef unsigned int   UINT32;
  typedef   signed short  INT16;
  typedef unsigned short UINT16;

  #ifndef _WINDEF_

    typedef signed char    INT8;
    typedef unsigned char UINT8;

  #endif

#endif // _E_FLOAT_TYPES_2004_02_28_H_
