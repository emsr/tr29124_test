/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_gamma_ratio test_gamma_ratio.cpp wrap_boost.cpp -lquadmath
./test_gamma_ratio > test_gamma_ratio.txt

*/

#include <ext/cmath>
#include <bits/float128_io.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <bits/summation.h>

#include "wrap_boost.h"
  /**
   * Return the four-gamma ratio.
   * @f[
   *    \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)\Gamma(a+b-c+n)}
   *      = 1 + sum_{m=1}^{M}\frac{(c-a)_m(c-b)_m}{m!(1+c-a-b-n)_m}
   *          + O(n^{-M-1})
   *      = 1 + sum_{m=1}^{M}\frac{(a-c)_m(b-c)_m}{m!(1-c-n)_m}
   *          + O(n^{-M-1})
   * @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_ratio_buhring(int __n, _Tp __a, _Tp __b, _Tp __c)
    {
      const auto _S_eps = 10 * __gnu_cxx::__epsilon(__a);
      if (std::real(1 - __c - __n) < _Tp{0})
	{
	  const auto _S_M = 20;
	  __gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> __sum;
	  auto __term = _Tp{1};
	  __sum += _Tp{__term};
	  auto __ca = __c - __a;
	  auto __cb = __c - __b;
	  auto __cabn = _Tp{1} + __c - __a - __b - __n;
	  for (int __m = 1; __m <= _S_M; ++__m)
	    {
	      auto __prev = __term;
	      __term *= __ca / __m * __cb / __cabn;
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * __sum())
		break;
	      if (std::abs(__term) > std::abs(__prev))
		break;
	      ++__ca;
	      ++__cb;
	      ++__cabn;
	    }
	  return __sum();
	}
      else if (std::real(1 + __c - __a - __b - __n) < _Tp{0})
	{
	  const auto _S_M = 20;
	  __gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> __sum;
	  auto __term = _Tp{1};
	  __sum += _Tp{__term};
	  auto __ac = __a - __c;
	  auto __bc = __b - __c;
	  auto __cn = _Tp{1} - __c - __n;
	  for (int __m = 1; __m <= _S_M; ++__m)
	    {
	      auto __prev = __term;
	      __term *= __ac / __m * __bc / __cn;
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * __sum())
		break;
	      if (std::abs(__term) > std::abs(__prev))
		break;
	      ++__ac;
	      ++__bc;
	      ++__cn;
	    }
	  return __sum();
	}
      else
	return __gnu_cxx::__infinity(__a);
    }

/**
 * Return the Stirling numbers of the first kind.
 *
 */
template<typename _Tp>
  _Tp
  stirling_number_1(int __n, int __k)
  {
    return _Tp{0};
  }

/**
 * Return the Stirling numbers of the first kind.
 *
 */
template<typename _Tp>
  _Tp
  stirling_number_1(int __n, _Tp __kappa)
  {
    return _Tp{0};
  }

/**
 * Return the ratio
 * @f[
 *   \frac{\Gamma(a+\alpha)}{\Gamma(a)}
 * @f]
 * for large @f$ a @f$.
 */
template<typename _Tp>
  _Tp
  __gamma_ratio_asymp()
  {
  }

template<typename _Tp>
  void
  test_gamma_ratio(_Tp proto = _Tp{})
  {
    //using _Val = _Tp;
    //using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();
  }

template<typename _Tp>
  void
  test_gamma_ratio_buhring(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    _Tp a = -11.7;
    _Tp b = -11.2;
    _Tp c = -11.4;
    auto b10 = __gamma_ratio_buhring(10, a, b, c);
    std::cout << b10 << '\n';
    auto b20 = __gamma_ratio_buhring(20, a, b, c);
    std::cout << b20 << '\n';
  }

int
main()
{
  test_gamma_ratio_buhring(1.0);

  test_gamma_ratio(1.0f);
  
  test_gamma_ratio(1.0);

  test_gamma_ratio(1.0l);

  test_gamma_ratio(1.0q);
}
