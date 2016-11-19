/*
$HOME/bin/bin/g++ -std=c++17 -Wall -Wextra -Wno-psabi -g -I. -o zeta_glob_deathmatch zeta_glob_deathmatch.cpp -lquadmath -L. -lwgsl
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./zeta_glob_deathmatch > zeta_glob_deathmatch.txt

$HOME/bin/bin/g++ -std=c++17 -Wall -Wextra -g -I. -o zeta_glob_deathmatch zeta_glob_deathmatch.cpp -lquadmath -L. -lwgsl
./zeta_glob_deathmatch > zeta_glob_deathmatch.txt
*/

#include<ext/cmath>
#include<limits>
#include<iostream>
#include<iomanip>

#include "wrap_gsl.h"

#define _GLIBCXX_MATH_NS std

    template<typename _Tp>
    _Tp
    __riemann_zeta_glob_old(_Tp __s)
    {
      _Tp __zeta = _Tp(0);

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp __max_bincoeff = std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1);

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (__s < _Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(__s,_Tp(2)) == _Tp(0))
            return _Tp(0);
          else
#endif
            {
              _Tp __zeta = __riemann_zeta_glob_old(_Tp(1) - __s);
              __zeta *= std::pow(_Tp(2)
                     * __gnu_cxx::__math_constants<_Tp>::__pi, __s)
                     * std::sin(__gnu_cxx::__math_constants<_Tp>::__pi_half * __s)
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(_Tp(1) - __s))
#else
                     * std::exp(__log_gamma(_Tp(1) - __s))
#endif
                     / __gnu_cxx::__math_constants<_Tp>::__pi;
              return __zeta;
            }
        }

      _Tp __num = _Tp(0.5L);
      const unsigned int __maxit = 10000;
      for (unsigned int __i = 0; __i < __maxit; ++__i)
        {
          bool __punt = false;
          _Tp __sgn = _Tp(1);
          _Tp __term = _Tp(0);
          for (unsigned int __j = 0; __j <= __i; ++__j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              _Tp __bincoeff =  _GLIBCXX_MATH_NS::lgamma(_Tp(1 + __i))
                              - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + __j))
                              - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + __i - __j));
#else
              _Tp __bincoeff =  __log_gamma(_Tp(1 + __i))
                              - __log_gamma(_Tp(1 + __j))
                              - __log_gamma(_Tp(1 + __i - __j));
#endif
              if (__bincoeff > __max_bincoeff)
                {
                  //  This only gets hit for x << 0.
                  __punt = true;
                  break;
                }
              __bincoeff = std::exp(__bincoeff);
              __term += __sgn * __bincoeff * std::pow(_Tp(1 + __j), -__s);
              __sgn *= _Tp(-1);
            }
          if (__punt)
            break;
          __term *= __num;
          __zeta += __term;
          if (std::abs(__term/__zeta) < __eps)
            break;
          __num *= _Tp(0.5L);
        }

      __zeta /= _Tp(1) - std::pow(_Tp(2), _Tp(1) - __s);

      return __zeta;
    }


    template<typename _Tp>
    _Tp
    __hurwitz_zeta_glob_old(_Tp __s, _Tp __a)
    {
      _Tp __zeta = _Tp(0);

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
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
              _Tp __bincoeff =  _GLIBCXX_MATH_NS::lgamma(_Tp(1 + __i))
                              - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + __j))
                              - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + __i - __j));
#else
              _Tp __bincoeff =  __log_gamma(_Tp(1 + __i))
                              - __log_gamma(_Tp(1 + __j))
                              - __log_gamma(_Tp(1 + __i - __j));
#endif
              if (__bincoeff > __max_bincoeff)
                {
                  //  This only gets hit for x << 0.
                  __punt = true;
                  break;
                }
              __bincoeff = std::exp(__bincoeff);
              __term += __sgn * __bincoeff * std::pow(_Tp(__j + __a), -__s);
              __sgn *= _Tp(-1);
            }
          if (__punt)
            break;
          __term /= _Tp(__i + 1);
          if (std::abs(__term / __zeta) < __eps)
            break;
          __zeta += __term;
        }

      __zeta /= __s - _Tp(1);

      return __zeta;
    }


    template<typename _Tp>
    _Tp
    __riemann_zeta_glob_new(_Tp __s)
    {
      _Tp __zeta = _Tp(0);

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp __max_bincoeff = std::exp(std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1));

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (__s < _Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(__s, _Tp(2)) == _Tp(0))
            return _Tp(0);
          else
#endif
            {
              _Tp __zeta = __riemann_zeta_glob_new(_Tp(1) - __s);
              __zeta *= std::pow(_Tp(2)
                     * __gnu_cxx::__math_constants<_Tp>::__pi, __s)
                     * std::sin(__gnu_cxx::__math_constants<_Tp>::__pi_half * __s)
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(_Tp(1) - __s))
#else
                     * std::exp(__log_gamma(_Tp(1) - __s))
#endif
                     / __gnu_cxx::__math_constants<_Tp>::__pi;
              return __zeta;
            }
        }

      _Tp __num = _Tp(0.25L);
      const unsigned int __maxit = 10000;
      __zeta = _Tp(0.5L);
      // This for loop starts at 1 because we already calculated the
      // value of the zeroeth order in __zeta above
      for (unsigned int __i = 1; __i < __maxit; ++__i)
        {
          bool __punt = false;
          _Tp __term = _Tp(1.0L);
          _Tp __bincoeff = _Tp(1.0L);
          // This for loop starts at 1 because we already calculated the value
          // of the zeroeth order in __term above.
          for (unsigned int __j = 1; __j <= __i; ++__j)
            {
              _Tp incr = _Tp(__i - __j + 1) / _Tp(__j);
              __bincoeff *= -incr;
              if(std::abs(__bincoeff) > __max_bincoeff )
                {
                  // This only gets hit for x << 0.
                  __punt = true;
                  break;
                }
              __term += __bincoeff * std::pow(_Tp(1 + __j), -__s);
            }
          if (__punt)
            break;
          __term *= __num;
          __zeta += __term;
          if (std::abs(__term / __zeta) < __eps)
            break;
          __num *= _Tp(0.5L);
        }

      __zeta /= _Tp(1) - std::pow(_Tp(2), _Tp(1) - __s);

      return __zeta;
    }


    template<typename _Tp>
    _Tp
    __hurwitz_zeta_glob_new(_Tp __s, _Tp __a)
    {
      _Tp __zeta = _Tp(0);

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp __max_bincoeff = std::exp(std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1));

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (__s < _Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(__s, _Tp(2)) == _Tp(0))
            return _Tp(0);
          else
#endif
            {
              _Tp __zeta = __hurwitz_zeta_glob_new(_Tp(1) - __s, __a);
              __zeta *= std::pow(_Tp(2)
                     * __gnu_cxx::__math_constants<_Tp>::__pi, __s)
                     * std::sin(__gnu_cxx::__math_constants<_Tp>::__pi_half * __s)
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(_Tp(1) - __s))
#else
                     * std::exp(__log_gamma(_Tp(1) - __s))
#endif
                     / __gnu_cxx::__math_constants<_Tp>::__pi;
              return __zeta;
            }
        }

      _Tp __num = _Tp(0.25L);
      const unsigned int __maxit = 10000;
      __zeta = _Tp(0.5L * std::pow(__a, -__s));
      // This for loop starts at 1 because we already calculated the
      // value of the zeroeth order in __zeta above
      for (unsigned int __i = 1; __i < __maxit; ++__i)
        {
          bool __punt = false;
          _Tp __term = std::pow(__a, -__s);
          _Tp __bincoeff = _Tp(1.0L);
          // This for loop starts at 1 because we already calculated the value
          // of the zeroeth order in __term above.
          for (unsigned int __j = 1; __j <= __i; ++__j)
            {
              _Tp incr = _Tp(__i - __j + 1) / _Tp(__j);
              __bincoeff *= -incr;
              if(std::abs(__bincoeff) > __max_bincoeff )
                {
                  // This only gets hit for x << 0.
                  __punt = true;
                  break;
                }
              __term += __bincoeff * std::pow(_Tp(__a + __j), -__s);
            }
          if (__punt)
            break;
          __term *= __num;
          __zeta += __term;
          if (std::abs(__term / __zeta) < __eps)
            break;
          __num *= _Tp(0.5L);
        }

      __zeta /= _Tp(1) - std::pow(__a + _Tp(1), _Tp(1) - __s);

      return __zeta;
    }

int
main()
{
  using _Tp = long double;

  std::cout.precision(std::numeric_limits<_Tp>::max_digits10);
  auto width = std::numeric_limits<_Tp>::max_digits10 + 6;

  // Test a Bernoulli thing for the regular zeta function.
  std::cout << '\n';
  for (auto is = -9; is < 100; ++is)
    {
      _Tp s = 0.1L * is;
      if (s == _Tp{1})
	continue;
      auto ozeta = __riemann_zeta_glob_old(s);
      auto nzeta = __riemann_zeta_glob_new(s);
      auto gzeta = gsl::riemann_zeta(s);
      std::cout << ' ' << std::setw(width) << s
		<< ' ' << std::setw(width) << gzeta
		<< ' ' << std::setw(width) << ozeta
		<< ' ' << std::setw(width) << nzeta
		<< ' ' << std::setw(width) << nzeta - ozeta
		<< ' ' << std::setw(width) << nzeta - gzeta
		<< std::endl;
    }

  std::cout << '\n';
  for (auto is = -9; is < 100; ++is)
    {
      _Tp s = 0.1L * is;
      if (s == _Tp{1})
	continue;
      //auto ozeta = __hurwitz_zeta_glob_old(s, _Tp{1});
      auto nzeta = __riemann_zeta_glob_new(s);
      auto hzeta = __hurwitz_zeta_glob_new(s, _Tp{1});
      std::cout << ' ' << std::setw(width) << s
		//<< ' ' << std::setw(width) << ozeta
		<< ' ' << std::setw(width) << nzeta
		<< ' ' << std::setw(width) << hzeta
		//<< ' ' << std::setw(width) << nzeta - ozeta
		<< ' ' << std::setw(width) << hzeta - nzeta
		<< std::endl;
    }

  std::cout << '\n';
  _Tp a = 3.0L;
  for (auto is = -9; is < 100; ++is)
    {
      _Tp s = 0.1L * is;
      if (s <= _Tp{1})
	continue;
      auto ozeta = __hurwitz_zeta_glob_old(s, a);
      auto hzeta = __hurwitz_zeta_glob_new(s, a);
      auto gzeta = gsl::hurwitz_zeta(s, a);
      std::cout << ' ' << std::setw(width) << s
		<< ' ' << std::setw(width) << gzeta
		<< ' ' << std::setw(width) << ozeta
		<< ' ' << std::setw(width) << hzeta
		<< ' ' << std::setw(width) << ozeta - gzeta
		<< ' ' << std::setw(width) << hzeta - ozeta
		<< ' ' << std::setw(width) << hzeta - gzeta
		<< std::endl;
    }
}
