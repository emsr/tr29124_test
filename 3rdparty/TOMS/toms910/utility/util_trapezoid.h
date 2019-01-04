
#ifndef _UTIL_TRAPEZOID_2008_09_06_H_
  #define _UTIL_TRAPEZOID_2008_09_06_H_

  #include <utility/util_ranged_function_operation.h>

  namespace Util
  {
    template<typename T> class RecursiveTrapezoidRule : public RangedFunctionOperation<T>
    {
    protected:

      RecursiveTrapezoidRule(const T& lo, const T& hi, const T& tol) : RangedFunctionOperation<T>(lo, hi, tol) { }

    public:

      virtual ~RecursiveTrapezoidRule() { }

    private:

      virtual T my_operation(void) const
      {
        INT32 n = static_cast<INT32>(1);

        T a = RangedFunctionOperation<T>::xlo;
        T b = RangedFunctionOperation<T>::xhi;
        T h = (b - a) / T(n);
        
        static const T one  = T(1.0);
        static const T half = T(0.5);

        T I = (function(a) + function(b)) * (h * half);

        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(31); k++)
        {
          h *= half;

          const T I0 = I;

          T sum(0);

          for(INT32 j = static_cast<INT32>(1); j <= n; j++)
          {
            sum += function(a + (T((j * 2) - 1) * h));
          }

          I = (I0 * half) + (h * sum);

          const T ratio = I0 / I;
          const T delta = ((ratio > one) ? (ratio - one) : (one - ratio));

          if((k > static_cast<INT32>(2)) && delta < RangedFunctionOperation<T>::eps)
          {
            FunctionOperation<T>::op_ok = true;
            break;
          }

          n = n * 2;
        }

        return I;
      }
    };
  }

#endif // _UTIL_TRAPEZOID_2008_09_06_H_
