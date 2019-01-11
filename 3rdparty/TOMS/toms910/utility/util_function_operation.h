
#ifndef _UTIL_FUNCTION_OPERATION_2009_10_27_H_
  #define _UTIL_FUNCTION_OPERATION_2009_10_27_H_

  #include <utility/util_function.h>

  namespace Util
  {
    template<typename T> class FunctionOperation : public Function<T>
    {
    protected:

      mutable bool op_ok;

    protected:

      FunctionOperation() : op_ok(false) { }

    public:

      virtual ~FunctionOperation() { }

      bool success(void) const { return op_ok; }
      T operation(void) const { return my_operation(); }

    protected:

      virtual T my_operation(void) const = 0;
    };
  }

#endif // _UTIL_FUNCTION_OPERATION_2009_10_27_H_
