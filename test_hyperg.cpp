/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hyperg test_hyperg.cpp wrap_boost.cpp -lquadmath
./test_hyperg > test_hyperg.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hyperg test_hyperg.cpp wrap_boost.cpp -lquadmath
./test_hyperg > test_hyperg.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_hyperg test_hyperg.cpp wrap_boost.cpp -lquadmath
./test_hyperg > test_hyperg.txt
*/

#include <bits/specfun.h>
#include <bits/float128_io.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

  template<typename _Tp>
    bool
    __isnan(_Tp __x)
    { return std::isnan(__x); }

  /**
   *   @brief Return the hypergeometric function @f$ _2F_1(a,b;c;x) @f$
   *   by series expansion.
   *
   *   The hypergeometric function is defined by
   *   @f[
   *     _2F_1(a,b;c;x) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}
   *                      \sum_{n=0}^{\infty}
   *                      \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)}
   *                      \frac{x^n}{n!}
   *   @f]
   *
   *   This works and it's pretty fast.
   *
   *   @param  __a  The first @a numerator parameter.
   *   @param  __b  The second @a numerator parameter.
   *   @param  __c  The @a denominator parameter.
   *   @param  __x  The argument of the confluent hypergeometric function.
   *   @return  The confluent hypergeometric function.
   */
  template<typename _Tp>
    _Tp
    __hyperg_series(_Tp __a, _Tp __b, _Tp __c, _Tp __x)
    {
      const auto __eps = __gnu_cxx::__epsilon(__x);
      const unsigned int _S_max_iter = 100000;
      auto __aint = __gnu_cxx::__fp_is_integer(__a);
      auto __bint = __gnu_cxx::__fp_is_integer(__b);
      auto __max_iter = _S_max_iter;
      if (__aint && __aint() < 0 && _S_max_iter > -__aint())
	__max_iter = -__aint();
      if (__bint && __bint() < 0 && _S_max_iter > -__bint())
	__max_iter = -__bint();

      auto __term = _Tp{1};
      auto __Fabc = _Tp{1};
      unsigned int __i;
      for (__i = 0; __i < __max_iter; ++__i)
	{
	  __term *= (__a + _Tp(__i)) * (__b + _Tp(__i)) * __x
		  / ((__c + _Tp(__i)) * _Tp(1 + __i));
	  if (std::abs(__term) < __eps)
	    break;
	  __Fabc += __term;
	}
      if (__i == __max_iter)
	std::__throw_runtime_error(__N("Series failed to converge "
				       "in __hyperg_series."));

      return __Fabc;
    }

  /**
   * @brief Return the hypergeometric function @f$ _2F_1(a,b;c;x) @f$.
   */
  template<typename _Tp>
    _Tp
    __hyperg(_Tp __a, _Tp __b, _Tp __c, _Tp __x, _Tp __rho = _Tp{0.5Q})
    {
      auto __aint = __gnu_cxx::__fp_is_integer(__a);
      auto __bint = __gnu_cxx::__fp_is_integer(__b);
      auto __cint = __gnu_cxx::__fp_is_integer(__c);
      auto __d = __c - __a - __b;
      auto __dint = __gnu_cxx::__fp_is_integer(__d);
      const auto __toler = _Tp{1000} * __gnu_cxx::__epsilon(__x);

      if (__isnan(__a) || __isnan(__b) || __isnan(__c) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__cint && __cint() <= 0)
	return __gnu_cxx::__infinity(__x);
      else if (std::abs(__c - __b) < __toler
	    || std::abs(__c - __a) < __toler)
	return std::pow(_Tp{1} - __x, __d);
      else if (std::abs(__x) <= __rho)
        return __hyperg_series(__a, __b, __c, __x);
      else if (std::abs(_Tp{1} - __x) <= __rho)
        {
	  if (__dint)
            {
	      auto __m = __dint();
	      auto __gammam = std::__detail::__gamma(_Tp(__m));
	      auto __gammaabm = std::__detail::__gamma(__a + __b + _Tp(__m));
	      auto __gammaam = std::__detail::__gamma(__a + _Tp(__m));
	      auto __gammabm = std::__detail::__gamma(__b + _Tp(__m));
	      auto __sum = _Tp{1};
	      auto __term = _Tp{1};
	      for (int __k = 0; __k < __m; ++__k)
	        {
		  __term *= (__a + _Tp(__m + __k)) * (__b + _Tp(__m + __k))
			  / _Tp(1 - __m + __k) / _Tp(__k) * (_Tp{1} - __x);
		  __sum += __term;
		}
	      return __sum * __gammam * __gammaabm / __gammaam / __gammabm;
	    }
	  else
            {
	      // This is where Buhring's gamma ratio might come in handy.
	      auto __gammaa = std::__detail::__gamma(__a);
	      auto __gammab = std::__detail::__gamma(__b);
	      auto __gammac = std::__detail::__gamma(__c);
	      auto __gammad = std::__detail::__gamma(__d);
	      auto __gammamd = std::__detail::__gamma(-__d);
	      auto __gammacma = std::__detail::__gamma(__c - __a);
	      auto __gammacmb = std::__detail::__gamma(__c - __b);
	      return __gammac * __gammad
		   * __hyperg_series(__a, __b, _Tp{1} - __d, _Tp{1} - __x)
		   / __gammacma / __gammacmb
		  + __gammac * __gammamd
		   * __hyperg_series(__c - __a, __c - __b, _Tp{1} + __d, _Tp{1} - __x)
		   / __gammaa / __gammab;
	    }
	}
      else
	return __gnu_cxx::__quiet_NaN(__x);
    }

/**
 * Test harness.
 */
template<typename _Tp>
  void
  test_hyperg(_Tp proto = _Tp{})
  {
    //using _Val = _Tp;
    //using _Real = std::__detail::__num_traits_t<_Val>;
    std::vector<_Tp> parm{_Tp{0.25}, _Tp{0.5}, _Tp{1}, _Tp{2}, _Tp{5}};

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "b"
	      << ' ' << std::setw(width) << "c"
	      << ' ' << std::setw(width) << "z"
	      << ' ' << std::setw(width) << "hyperg0"
	      << ' ' << std::setw(width) << "hyperg"
	      << '\n';
    int i_min = -99;
    for (auto a : parm)
      for (auto b : parm)
	for (auto c : parm)
	  for (int i = i_min; i <= +99; ++i)
	    {
	      auto z = _Tp{0.01Q} * i;
	      auto hyperg0 = __gnu_cxx::hyperg(a, b, c, z);
	      auto hyperg = __hyperg(a, b, c, z);
	      std::cout << ' ' << std::setw(width) << a
			<< ' ' << std::setw(width) << b
			<< ' ' << std::setw(width) << c
			<< ' ' << std::setw(width) << z
			<< ' ' << std::setw(width) << hyperg0
			<< ' ' << std::setw(width) << hyperg
			<< ' ' << std::setw(width) << hyperg - hyperg0
			<< '\n';
	    }
  }

int
main()
{
  test_hyperg(1.0F);

  test_hyperg(1.0);

  test_hyperg(1.0L);
}
