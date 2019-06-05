/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <ext/math_const.h>

/**
 * Compute the inverse Gudermann function relating a circular argument
 * to a hyperbolic argument using the functional definition:
 * @f[
 *    invgd(\phi) = asinh(tan(\phi))
 * @f]
 */
template<typename _Tp>
  _Tp
  invgd_trig(_Tp phi)
  {
    return std::asinh(std::tan(phi));
  }

/**
 * Compute the inverse Gudermann function relating a circular argument
 * to a hyperbolic argument using the series:
 * @f[
 *    invgd(\phi) = \sum_{k=0}^{\infty}\frac{(-1)^kE_{2k}}{(2k+1)!}
 * @f]
 * wnere @f$ E_n @f$ are the Euler numbers.
 */
template<typename _Tp>
  _Tp
  invgd_series(_Tp x)
  {
    const auto xx = -x * x;
    auto term = x;
    auto sum = term;
    for (int i = 1; i < 20; ++i)
      {
	term *= xx
	     / _Tp(2 * i) / _Tp(2 * i + 1);
	sum += term * __gnu_cxx::euler<_Tp>(2 * i);
      }
    return sum;
  }

/**
 * 
 * @f[
 *    invgd(\phi) = F(1,\phi)
 * @f]
 */
template<typename _Tp>
  _Tp
  invgd_ellint_2(_Tp x)
  {
    return std::ellint_2(_Tp{1}, x);
  }

/**
 * Compute the Gudermann function relating a hyperbolic argument
 * to a circular argument using the functional definition:
 * @f[
 *    gd(x) = atan(sinh(x))
 * @f]
 */
template<typename _Tp>
  _Tp
  gd_trig(_Tp x)
  {
    return std::atan(std::sinh(x));
  }

/**
 * Compute the inverse Gudermann function relating a circular argument
 * to a hyperbolic argument using the series:
 * @f[
 *    gd(x) = \sum_{k=0}^{\infty}\frac{E_{2k}}{(2k+1)!}
 * @f]
 * wnere @f$ E_n @f$ are the Euler numbers.
 */
template<typename _Tp>
  _Tp
  gd_series(_Tp x)
  {
    const auto xx = x * x;
    auto term = x;
    auto sum = term;
    for (int i = 1; i < 20; ++i)
      {
	term *= xx
	     / _Tp(2 * i) / _Tp(2 * i + 1);
	sum += term * __gnu_cxx::euler<_Tp>(2 * i);
      }
    return sum;
  }

template<typename _Tp>
  _Tp
  gd(_Tp x)
  {
    const auto _S_pi = __const_2_pi(x);

    if (std::isnan(x))
      return x;
    else if (x == -std::numeric_limits<_Tp>::infinity())
      return -_S_pi / _Tp{2};
    else if (x == +std::numeric_limits<_Tp>::infinity())
      return +_S_pi / _Tp{2};
    else
      return gd_trig(x);
  }

template<typename _Tp>
  void
  test_gudermannian()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    const auto _S_pi = __gnu_cxx::__const_2_pi<_Tp>();
    const auto del = _S_pi / 400;
    for (int i = -200; i <= 200; ++i)
      {
	auto x = i * del;
	auto gd_fun = gd_trig(x);
	auto gd_ser = gd_series(x);
	auto gd_diff = gd_ser - gd_fun;
	auto igd_fun = invgd_trig(x);
	auto igd_ser = invgd_series(x);
	auto igd_F   = invgd_ellint_2(x);
	auto test_gd = invgd_series(gd_fun) - x;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << gd_ser
		  << ' ' << std::setw(w) << gd_fun
		  << ' ' << std::setw(w) << gd_diff
		  << ' ' << std::setw(w) << igd_ser
		  << ' ' << std::setw(w) << igd_fun
		  << ' ' << std::setw(w) << igd_F
		  << ' ' << std::setw(w) << test_gd
		  << '\n';
      }
  }

int
main()
{
  test_gudermannian<double>();
}
