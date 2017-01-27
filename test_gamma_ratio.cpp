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

  enum __buhring_mode
  {
    automatic,
    equation2p7,
    equation3p1
  };


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
    __gamma_ratio_buhring(int __n, _Tp __a, _Tp __b, _Tp __c,
			  int _S_M = 20, __buhring_mode mode = automatic)
    {
      const auto _S_eps = 10 * __gnu_cxx::__epsilon(__a);
      if (mode == equation2p7
	|| (mode == automatic && std::real(1 - __c - __n) < _Tp{0}))
	{
	  //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> __sum;
	  __gnu_cxx::_BasicSum<_Tp> __sum;
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
	      if (mode == automatic
		&& (std::abs(__term) < _S_eps * __sum()
		 || std::abs(__term) > std::abs(__prev)))
		break;
	      ++__ca;
	      ++__cb;
	      ++__cabn;
	    }
	  return __sum();
	}
      else if (mode == equation3p1
	|| (mode == automatic && std::real(1 + __c - __a - __b - __n) < _Tp{0}))
	{
	  //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> __sum;
	  __gnu_cxx::_BasicSum<_Tp> __sum;
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
	      if (mode == automatic
		&& (std::abs(__term) < _S_eps * __sum()
		 || std::abs(__term) > std::abs(__prev)))
		break;
	      if (mode == automatic && std::abs(__term) > std::abs(__prev))
		break;
	      ++__ac;
	      ++__bc;
	      ++__cn;
	    }
	  return __sum();
	}
      else
	{
	  // Average 2.7 and 3.1 sums.
	  __gnu_cxx::_BasicSum<_Tp> __sum2p7;
	  auto __term2p7 = _Tp{1};
	  __sum2p7 += _Tp{__term2p7};
	  auto __ca = __c - __a;
	  auto __cb = __c - __b;
	  auto __cabn = _Tp{1} + __c - __a - __b - __n;

	  __gnu_cxx::_BasicSum<_Tp> __sum3p1;
	  auto __term3p1 = _Tp{1};
	  __sum3p1 += _Tp{__term3p1};
	  auto __ac = __a - __c;
	  auto __bc = __b - __c;
	  auto __cn = _Tp{1} - __c - __n;

	  bool __conv2p7 = false;
	  bool __conv3p1 = false;
	  for (int __m = 1; __m <= _S_M; ++__m)
	    {
	      auto __prev2p7 = __term2p7;
	      __term2p7 *= __ca / __m * __cb / __cabn;
	      __sum2p7 += __term2p7;
	      ++__ca;
	      ++__cb;
	      ++__cabn;
	      if (std::abs(__term2p7) < _S_eps * __sum2p7()
		|| std::abs(__term2p7) > std::abs(__prev2p7))
		__conv2p7 = true;

	      auto __prev3p1 = __term3p1;
	      __term3p1 *= __ac / __m * __bc / __cn;
	      __sum3p1 += __term3p1;
	      ++__ac;
	      ++__bc;
	      ++__cn;
	      if (std::abs(__term3p1) < _S_eps * __sum3p1()
		|| std::abs(__term3p1) > std::abs(__prev3p1))
		__conv3p1 = true;

	      if (__conv2p7 && __conv3p1)
		break;
	    }
	  return (__sum2p7() + __sum3p1()) / _Tp{2}
		* __gnu_cxx::sin_pi(__c + __n) * __gnu_cxx::sin_pi(__a + __b - __c + __n)
		/ __gnu_cxx::sin_pi(__a + __n) * __gnu_cxx::sin_pi(__b + __n);
	}
    }

  /**
   * Juggle gamma special cases - from hyperg.
   */
  template<typename _Tp>
    _Tp
    __hyperg_reflect(_Tp __a, _Tp __b, _Tp __c)
    {
      const auto __d = __c - __a - __b;
      const auto __intd = __gnu_cxx::__fp_is_integer(__d);
      auto __F1 = _Tp{0};
      if (__intd)
	{
	  if (__intd() <= 0)
	    {
	      __F1 = _Tp{0};
	    }
	  else
	    {
	      __F1 = _Tp{0};
	    }
	}
      else
	{
	  __F1 = _Tp{0};
	}
    }

/**
 * Return the gamma ratio by asymptotic series:
 * @f[
 *    \frac{\Gamma(z+\alpha)}{\Gamma(z+\beta)}
 *  ~ z^{\alpha-\beta}\sum_{n=0}^{\infty}\frac{(\beta-\alpha)_n}{n!}
 *         {}_2F_0(-n, z+\beta;;1/z)
 * @f]
 * The hypergeometric function @f$ {}_2F_0(a,b;;z) @f$
 * with nonpositive numeratorial parameter n
 * is a polynomial with a simple TTRR:
 * @f[
 *  z {}_2F_0(-n-1,z+\beta;;1/z) = -(n+\beta) {}_2F_0(-n,z+\beta;;1/z)
 *      + n {}_2F_0(-n+1,z+\beta;;1/z)
 * @f]
 * where @f$ {}_2F_0(-0,z+\beta;;1/z) = 1 @f$
 * and @f$ {}_2F_0(-1,z+\beta;;1/z) = -\beta/z @f$.
 */
template<typename _Tp>
  _Tp
  gamma_ratio_2f0(_Tp __alpha, _Tp __beta, _Tp __z)
  {
    return _Tp{0};
  }

/**
 * Return the gamma ratio by asymptotic series:
 * @f[
 *    \frac{\Gamma(z+\alpha)}{\Gamma(z+\beta)}
 *  ~ z^{\alpha-\beta}\sum_{n=0}^{\infty}\frac{(\alpha-\beta)_n}{n!}
 *         {}_1F_1(-n; z+\alpha; z)
 * @f]
 * The hypergeometric function @f$ {}_1F_1(a;b;z) @f$
 * with nonpositive numeratorial parameter n
 * is a polynomial with a simple TTRR:
 * @f[
 *  (z+\alpha+n) {}_1F_1(-n-1;z+\alpha;z) = -(2n+\alpha) {}_1F_1(-n;z+\alpha;z)
 *      - n {}_1F_1(-n+1;z+\alpha;z)
 * @f]
 * where @f$ {}_1F_1(-0;z+\alpha;z) = 1 @f$
 * and @f$ {}_1F_1(-1;z+\alpha;z) = \alpha/(z+\alpha) @f$.
 */
template<typename _Tp>
  _Tp
  gamma_ratio_1f1(_Tp __alpha, _Tp __beta, _Tp __z)
  {
    return _Tp{0};
  }

/**
 * Return the gamma ratio by asymptotic series:
 * @f[
 *    \frac{\Gamma(z+\alpha)}{\Gamma(z+\beta)}
 *  ~ z^{\alpha-\beta}\sum_{n=0}^{\infty}
 *      \frac{(\alpha-\beta+1-n)_n B_n^{\alpha-\beta+1}(\alpha)}{n!}z^{-n}
 * @f]
 * where @f$@f$ is the Norlund or generalized Bernoulli polynomial.
 */
  _Tp
  gamma_ratio_erdelyi_tricomi(_Tp __alpha, _Tp __beta, _Tp __z)
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
    int i_min = -200;
    for (auto a : parm)
      for (auto b : parm)
	for (auto c : parm)
	  for (int i = i_min; i <= +200; ++i)
	    {
	      auto n = _Tp{0.1Q} * i;
	      //auto gamrat0 = __gamma_ratio_log(a, b, c, z);
	      auto gamrat = __gamma_ratio_buhring(n, a, b, c);
	      std::cout << ' ' << std::setw(width) << n
			<< ' ' << std::setw(width) << a
			<< ' ' << std::setw(width) << b
			<< ' ' << std::setw(width) << c
			//<< ' ' << std::setw(width) << gamrat0
			<< ' ' << std::setw(width) << gamrat
			//<< ' ' << std::setw(width) << gamrat - gamrat0
			<< '\n';
	    }
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

    std::cout << '\n';

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b10_2p7 = __gamma_ratio_buhring(10, a, b, c, M, equation2p7);
	auto b10_3p1 = __gamma_ratio_buhring(10, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b10_2p7
		  << std::setw(width) << b10_3p1
		  << '\n';
      }

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b20_2p7 = __gamma_ratio_buhring(20, a, b, c, M, equation2p7);
	auto b20_3p1 = __gamma_ratio_buhring(20, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b20_2p7
		  << std::setw(width) << b20_3p1
		  << '\n';
      }

    a = 11.7;
    b = 11.2;
    c = 11.4;

    std::cout << '\n';

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b10_2p7 = __gamma_ratio_buhring(-15, a, b, c, M, equation2p7);
	auto b10_3p1 = __gamma_ratio_buhring(-15, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b10_2p7
		  << std::setw(width) << b10_3p1
		  << '\n';
      }

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b20_2p7 = __gamma_ratio_buhring(-5, a, b, c, M, equation2p7);
	auto b20_3p1 = __gamma_ratio_buhring(-5, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b20_2p7
		  << std::setw(width) << b20_3p1
		  << '\n';
      }
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
