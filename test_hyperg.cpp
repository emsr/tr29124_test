/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hyperg test_hyperg.cpp wrap_boost.cpp -lquadmath
./test_hyperg > test_hyperg.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hyperg test_hyperg.cpp wrap_boost.cpp -lquadmath
./test_hyperg > test_hyperg.txt

g++ -std=gnu++17 -g -Wall -Wextra -DNO_LOGBQ -I. -o test_hyperg test_hyperg.cpp wrap_boost.cpp -lquadmath
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

      auto __term = _Tp{1};
      auto __Fabc = _Tp{1};
      const unsigned int __max_iter = 100000;
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
   *   @brief Return the hypergeometric function @f$ _2F_1(a,b;c;x) @f$.
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
    }

/**
 * Test harness.
 */
template<typename _Tp>
  void
  test_hyperg(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "b"
	      << ' ' << std::setw(width) << "c"
	      << ' ' << std::setw(width) << "z"
	      << ' ' << std::setw(width) << "hyperg"
	      //<< ' ' << std::setw(width) << "delta_boost"
	      << '\n';
    int i_min = -200;
    for (auto a : {0.0, 0.25, 0.5, 1.0, 2.0, 5.0})
      for (auto b : {0.0, 0.25, 0.5, 1.0, 2.0, 5.0})
	for (auto c : {0.0, 0.25, 0.5, 1.0, 2.0, 5.0})
	  for (int i = i_min; i <= +500; ++i)
	    {
	      auto z = _Tp{0.10Q} * i;
	      auto hyperg0 = __gnu_cxx::hyperg(a, b, c, z);
	      std::cout << ' ' << std::setw(width) << a
			<< ' ' << std::setw(width) << b
			<< ' ' << std::setw(width) << c
			<< ' ' << std::setw(width) << z
			<< ' ' << std::setw(width) << hyperg0
//			<< ' ' << std::setw(width) << (glgam - blgam) / std::abs(blgam)
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
