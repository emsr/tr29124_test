/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <wrap_gsl.h>
#include <emsr/math_constants.h>

#define _GLIBCXX_MATH_NS std

    template<typename _Tp>
    _Tp
    riemann_zeta_glob_old(_Tp s)
    {
      _Tp zeta = _Tp(0);

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp max_binom = std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1);

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (s < _Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(s,_Tp(2)) == _Tp(0))
            return _Tp(0);
          else
#endif
            {
              _Tp zeta = riemann_zeta_glob_old(_Tp(1) - s);
              zeta *= std::pow(_Tp(2)
                     * emsr::pi_v<_Tp>, s)
                     * std::sin(emsr::pi_v<_Tp> * s / _Tp{2})
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(_Tp(1) - s))
#else
                     * std::exp(log_gamma(_Tp(1) - s))
#endif
                     / emsr::pi_v<_Tp>;
              return zeta;
            }
        }

      _Tp num = _Tp(0.5L);
      const unsigned int maxit = 10000;
      for (unsigned int i = 0; i < maxit; ++i)
        {
          bool punt = false;
          _Tp sgn = _Tp(1);
          _Tp term = _Tp(0);
          for (unsigned int j = 0; j <= i; ++j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              _Tp binom = _GLIBCXX_MATH_NS::lgamma(_Tp(1 + i))
                          - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + j))
                          - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + i - j));
#else
              _Tp binom = log_gamma(_Tp(1 + i))
                          - log_gamma(_Tp(1 + j))
                          - log_gamma(_Tp(1 + i - j));
#endif
              if (binom > max_binom)
                {
                  //  This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              binom = std::exp(binom);
              term += sgn * binom * std::pow(_Tp(1 + j), -s);
              sgn *= _Tp(-1);
            }
          if (punt)
            break;
          term *= num;
          zeta += term;
          if (std::abs(term/zeta) < eps)
            break;
          num *= _Tp(0.5L);
        }

      zeta /= _Tp(1) - std::pow(_Tp(2), _Tp(1) - s);

      return zeta;
    }


    template<typename _Tp>
    _Tp
    hurwitz_zeta_glob_old(_Tp s, _Tp a)
    {
      _Tp zeta = _Tp(0);

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp max_binom = std::numeric_limits<_Tp>::max_exponent10
                            * std::log(_Tp(10)) - _Tp(1);

      const unsigned int maxit = 10000;
      for (unsigned int i = 0; i < maxit; ++i)
        {
          bool punt = false;
          _Tp sgn = _Tp(1);
          _Tp term = _Tp(0);
          for (unsigned int j = 0; j <= i; ++j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              _Tp binom =  _GLIBCXX_MATH_NS::lgamma(_Tp(1 + i))
                          - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + j))
                          - _GLIBCXX_MATH_NS::lgamma(_Tp(1 + i - j));
#else
              _Tp binom =  log_gamma(_Tp(1 + i))
                          - log_gamma(_Tp(1 + j))
                          - log_gamma(_Tp(1 + i - j));
#endif
              if (binom > max_binom)
                {
                  //  This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              binom = std::exp(binom);
              term += sgn * binom * std::pow(_Tp(j + a), -s);
              sgn *= _Tp(-1);
            }
          if (punt)
            break;
          term /= _Tp(i + 1);
          if (std::abs(term / zeta) < eps)
            break;
          zeta += term;
        }

      zeta /= s - _Tp(1);

      return zeta;
    }


    template<typename _Tp>
    _Tp
    riemann_zeta_glob_new(_Tp s)
    {
      _Tp zeta = _Tp(0);

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp max_binom = std::exp(std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1));

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (s < _Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(s, _Tp(2)) == _Tp(0))
            return _Tp(0);
          else
#endif
            {
              _Tp zeta = riemann_zeta_glob_new(_Tp(1) - s);
              zeta *= std::pow(_Tp(2)
                     * emsr::pi_v<_Tp>, s)
                     * std::sin(emsr::pi_v<_Tp> * s / _Tp{2})
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(_Tp(1) - s))
#else
                     * std::exp(log_gamma(_Tp(1) - s))
#endif
                     / emsr::pi_v<_Tp>;
              return zeta;
            }
        }

      _Tp num = _Tp(0.25L);
      const unsigned int maxit = 10000;
      zeta = _Tp(0.5L);
      // This for loop starts at 1 because we already calculated the
      // value of the zeroeth order in zeta above
      for (unsigned int i = 1; i < maxit; ++i)
        {
          bool punt = false;
          _Tp term = _Tp(1.0L);
          _Tp binom = _Tp(1.0L);
          // This for loop starts at 1 because we already calculated the value
          // of the zeroeth order in term above.
          for (unsigned int j = 1; j <= i; ++j)
            {
              _Tp incr = _Tp(i - j + 1) / _Tp(j);
              binom *= -incr;
              if(std::abs(binom) > max_binom )
                {
                  // This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              term += binom * std::pow(_Tp(1 + j), -s);
            }
          if (punt)
            break;
          term *= num;
          zeta += term;
          if (std::abs(term / zeta) < eps)
            break;
          num *= _Tp(0.5L);
        }

      zeta /= _Tp(1) - std::pow(_Tp(2), _Tp(1) - s);

      return zeta;
    }


    template<typename _Tp>
    _Tp
    hurwitz_zeta_glob_new(_Tp s, _Tp a)
    {
      _Tp zeta = _Tp(0);

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      //  Max e exponent before overflow.
      const _Tp max_binom = std::exp(std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1));

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (s < _Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(s, _Tp(2)) == _Tp(0))
            return _Tp(0);
          else
#endif
            {
              _Tp zeta = hurwitz_zeta_glob_new(_Tp(1) - s, a);
              zeta *= std::pow(_Tp(2)
                     * emsr::pi_v<_Tp>, s)
                     * std::sin(emsr::pi_v<_Tp> * s / _Tp{2})
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(_Tp(1) - s))
#else
                     * std::exp(log_gamma(_Tp(1) - s))
#endif
                     / emsr::pi_v<_Tp>;
              return zeta;
            }
        }

      _Tp num = _Tp(0.25L);
      const unsigned int maxit = 10000;
      zeta = _Tp(0.5L * std::pow(a, -s));
      // This for loop starts at 1 because we already calculated the
      // value of the zeroeth order in zeta above
      for (unsigned int i = 1; i < maxit; ++i)
        {
          bool punt = false;
          _Tp term = std::pow(a, -s);
          _Tp binom = _Tp(1.0L);
          // This for loop starts at 1 because we already calculated the value
          // of the zeroeth order in term above.
          for (unsigned int j = 1; j <= i; ++j)
            {
              _Tp incr = _Tp(i - j + 1) / _Tp(j);
              binom *= -incr;
              if(std::abs(binom) > max_binom )
                {
                  // This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              term += binom * std::pow(_Tp(a + j), -s);
            }
          if (punt)
            break;
          term *= num;
          zeta += term;
          if (std::abs(term / zeta) < eps)
            break;
          num *= _Tp(0.5L);
        }

      zeta /= _Tp(1) - std::pow(a + _Tp(1), _Tp(1) - s);

      return zeta;
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
      auto ozeta = riemann_zeta_glob_old(s);
      auto nzeta = riemann_zeta_glob_new(s);
      auto gzeta = gsl::riemann_zeta(s);
      std::cout << ' ' << std::setw(width) << s
		<< ' ' << std::setw(width) << gzeta
		<< ' ' << std::setw(width) << ozeta
		<< ' ' << std::setw(width) << nzeta
		<< ' ' << std::setw(width) << nzeta - ozeta
		<< ' ' << std::setw(width) << nzeta - gzeta
		<< '\n';
    }

  std::cout << '\n';
  for (auto is = -9; is < 100; ++is)
    {
      _Tp s = 0.1L * is;
      if (s == _Tp{1})
	continue;
      //auto ozeta = hurwitz_zeta_glob_old(s, _Tp{1});
      auto nzeta = riemann_zeta_glob_new(s);
      auto hzeta = hurwitz_zeta_glob_new(s, _Tp{1});
      std::cout << ' ' << std::setw(width) << s
		//<< ' ' << std::setw(width) << ozeta
		<< ' ' << std::setw(width) << nzeta
		<< ' ' << std::setw(width) << hzeta
		//<< ' ' << std::setw(width) << nzeta - ozeta
		<< ' ' << std::setw(width) << hzeta - nzeta
		<< '\n';
    }

  std::cout << '\n';
  _Tp a = 3.0L;
  for (auto is = -9; is < 100; ++is)
    {
      _Tp s = 0.1L * is;
      if (s <= _Tp{1})
	continue;
      auto ozeta = hurwitz_zeta_glob_old(s, a);
      auto hzeta = hurwitz_zeta_glob_new(s, a);
      auto gzeta = gsl::hurwitz_zeta(s, a);
      std::cout << ' ' << std::setw(width) << s
		<< ' ' << std::setw(width) << gzeta
		<< ' ' << std::setw(width) << ozeta
		<< ' ' << std::setw(width) << hzeta
		<< ' ' << std::setw(width) << ozeta - gzeta
		<< ' ' << std::setw(width) << hzeta - ozeta
		<< ' ' << std::setw(width) << hzeta - gzeta
		<< '\n';
    }
}
