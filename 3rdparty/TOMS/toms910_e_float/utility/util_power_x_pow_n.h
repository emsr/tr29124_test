
#ifndef _UTIL_POWER_X_POW_N_2009_11_23_H_
  #define _UTIL_POWER_X_POW_N_2009_11_23_H_

  namespace Util
  {
    template<typename T> inline T x_pow_n_template(const T& t, const INT64 p)
    {
      // Compute the pure power of typename T t^p. Binary splitting of the power is
      // used. The resulting computational complexity has the order of log2[abs(p)].

      if(p < static_cast<INT64>(0))
      {
        return T(1) / x_pow_n_template(t, -p);
      }

      switch(p)
      {
        case static_cast<INT64>(0):
          return T(1);

        case static_cast<INT64>(1):
          return t;

        case static_cast<INT64>(2):
          return t * t;

        case static_cast<INT64>(3):
          return (t * t) * t;

        case static_cast<INT64>(4):
        {
          const T tsq = t * t;
          return tsq * tsq;
        }

        default:
        {
          T value(t);
          INT64 n;

          for(n = static_cast<INT64>(1); n <= static_cast<INT64>(p / static_cast<INT64>(2)); n *= static_cast<INT64>(2))
          {
            value *= value;
          }

          const INT64 p_minus_n = static_cast<INT64>(p - n);

          // Call the function recursively for computing the remaining power of n.
          return ((p_minus_n == static_cast<INT64>(0)) ? value
                                                       : value * Util::x_pow_n_template(t, p_minus_n));
        }
      }
    }
  }

#endif // _UTIL_POWER_X_POW_N_2009_11_23_H_
