
#ifndef _UTIL_RANGED_FUNCTION_OPERATION_2009_10_27_H_
  #define _UTIL_RANGED_FUNCTION_OPERATION_2009_10_27_H_

  #include <utility/util_function_operation.h>

  namespace Util
  {
    template<typename T> class RangedFunctionOperation : public FunctionOperation<T>
    {
    protected:

      const T xlo;
      const T xhi;
      const T eps;

    protected:

      RangedFunctionOperation(const T& lo,
                              const T& hi,
                              const T& tol) : xlo(lo),
                                              xhi(hi),
                                              eps(tol) { }

    public:

      virtual ~RangedFunctionOperation() { }
    };
  }

#endif // _UTIL_RANGED_FUNCTION_OPERATION_2009_10_27_H_
