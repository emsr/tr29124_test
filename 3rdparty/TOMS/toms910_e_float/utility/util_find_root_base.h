
#ifndef _UTIL_FIND_ROOT_BASE_2009_10_31_H_
  #define _UTIL_FIND_ROOT_BASE_2009_10_31_H_

  #include <utility/util_ranged_function_operation.h>

  namespace Util
  {
    template<typename T> class FindRootBase : public RangedFunctionOperation<T>
    {
    protected:
    
      FindRootBase(const T& lo,
                   const T& hi,
                   const T& tol) : RangedFunctionOperation<T>(lo, hi, tol) { }

    public:

      virtual ~FindRootBase() { }
    };
  }

#endif // _UTIL_FIND_ROOT_BASE_2009_10_31_H_
