
#ifndef _UTIL_FIND_ROOT_NEWTON_RAPHSON_2009_10_27_H_
  #define _UTIL_FIND_ROOT_NEWTON_RAPHSON_2009_10_27_H_

  #include <utility/util_find_root_base.h>

  namespace Util
  {
    template<typename T> class FindRootNewtonRaphson : public FindRootBase<T>
    {
    protected:
    
      FindRootNewtonRaphson(const T& lo,
                            const T& hi,
                            const T& tol) : FindRootBase<T>(lo, hi, tol) { }

    public:

      virtual ~FindRootNewtonRaphson() { }

    public:

      T derivative(const T& x) const { return my_derivative(x); }

    private:

      virtual T my_derivative(const T& x) const = 0;

      virtual T my_operation(void) const
      {
        FunctionOperation<T>::op_ok = true;

        T x = T(0.5) * (FindRootBase<T>::xlo + FindRootBase<T>::xhi);

        for(INT32 j = static_cast<INT32>(0); j < static_cast<INT32>(128); j++)
        {
          const T dx = function(x) / derivative(x);
          
          x -= dx;

          if((FindRootBase<T>::xlo - x) * (x - FindRootBase<T>::xhi) < ef::zero())
          {
            FunctionOperation<T>::op_ok = false;
            break;
          }

          const T delta = (dx < T(0) ? -dx : dx);
          
          if(delta < FindRootBase<T>::eps)
          {
            break;
          }
        }

        return FunctionOperation<T>::op_ok ? x : T(0);
      }
    };
  }

#endif // _UTIL_FIND_ROOT_NEWTON_RAPHSON_2009_10_27_H_
