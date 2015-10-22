
#include <tr1/cmath>
#include <iostream>
#include <iomanip>




    template<typename _Tp>
    _Tp
    __hurwitz_zeta_sum(const _Tp __a, const _Tp __s)
    {
      _Tp __zeta = _Tp(0);
      const unsigned int __maxit = 10000;
      for (unsigned int __i = 0; __i < __maxit; ++__i)
        {
          if (std::abs(_Tp(__i) + __a) < std::numeric_limits<_Tp>::epsilon())
            continue;
          const _Tp __term = std::pow(_Tp(__i) + __a, -__s);
          __zeta += __term;
          if (__term / __zeta < std::numeric_limits<_Tp>::epsilon())
            break;
        }
      return __zeta;
    }



    template<typename _Tp>
    _Tp
    __hurwitz_zeta_glob(const _Tp __a, const _Tp __s)
    {
      _Tp __zeta = _Tp(0);

      //  Max e exponent before overflow.
      const _Tp __max_bincoeff = std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1);

      const unsigned int __maxit = 10000;
      for (unsigned int __i = 0; __i < __maxit; ++__i)
        {
          bool __punt = false;
          _Tp __sgn = _Tp(1);
          _Tp __term = _Tp(0);
          for (unsigned int __j = 0; __j <= __i; ++__j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              _Tp __bincoeff =  std::tr1::lgamma(_Tp(1 + __i))
                              - std::tr1::lgamma(_Tp(1 + __j))
                              - std::tr1::lgamma(_Tp(1 + __i - __j));
#else
              _Tp __bincoeff =  std::tr1::__detail::__log_gamma(_Tp(1 + __i))
                              - std::tr1::__detail::__log_gamma(_Tp(1 + __j))
                              - std::tr1::__detail::__log_gamma(_Tp(1 + __i - __j));
#endif
              __bincoeff += (_Tp(1) - __s) * std::log(_Tp(__a + __j));
              if (__bincoeff > __max_bincoeff)
                {
                  //  This only gets hit for x << 0.
                  __punt = true;
                  break;
                }
              __bincoeff = std::exp(__bincoeff);
              __term += __sgn * __bincoeff;
              __sgn *= _Tp(-1);
            }
          if (__punt)
            break;
          __term /= _Tp(__i + 1);
          __zeta += __term;
 std::cout << __zeta << std::endl;
          if (std::abs(__term / __zeta) < 1000 * std::numeric_limits<_Tp>::epsilon())
            break;
        }

      __zeta /= __s - _Tp(1);

      return __zeta;
    }


int main(int, char**)
{

  std::cout.precision(16);

  for (unsigned int i = 1; i <= 100; ++i)
    {
      double x = i * 1.0;
      std::cout << std::setw(24) << x;
      std::cout << std::setw(24) << std::tr1::riemann_zeta(x);
      for (unsigned int n = 0; n < 3; ++n)
        {
          std::cout << std::setw(24) << __hurwitz_zeta_sum(double(n), x);
        }
      std::cout << std::endl;
    }

return 1;

  for (unsigned int i = 1; i <= 100; ++i)
    {
      double x = i * 1.0;
      for (unsigned int n = 0; n < 3; ++n)
        {
          std::cout << std::setw(24) << __hurwitz_zeta_glob(double(n), x);
        }
      std::cout << std::endl;
    }

  return 0;
}

