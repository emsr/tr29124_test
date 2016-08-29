/*
$HOME/bin_tr29124/bin/g++ -std=gnu++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -g -Wall -Wextra -Wno-compare-reals -o test_struve test_struve.cpp wrap_burkhardt.cpp burkhardt/special_functions.f90 -lgfortran
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_struve > test_struve.new

g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -Wall -Wextra -o test_struve test_struve.cpp wrap_burkhardt.cpp burkhardt/special_functions.f90 -lgfortran
./test_struve > test_struve.txt
*/

#include <cassert>
#include <cmath> // There are issues with <complex> inclusion if this isn't up here!
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <string>
#include <ext/math_const.h>
#include <bits/numeric_limits.h>
#include <bits/specfun_util.h>
#include <bits/complex_util.h>
#include <bits/summation.h>

#include "wrap_burkhardt.h"

namespace std
{
namespace __detail
{

  /**
   * An enum to dispatch Struve function summation.
   */
  enum  _StruveType
  : int
  {
    _StruveH,
    _StruveK,
    _StruveL,
    _StruveM
  };

  /**
   * Return either the Struve function of the first kind
   * @f$ \boldmath{H}_\nu(x) @f$ or the modified Struve function
   * of the first kind @f$ \boldmath{L}_\nu(x) @f$
   * depending on whether @c sign is -1 or +1 respectively.
   */
  template<_StruveType _Type, typename _Tp>
    _Tp
    __struve_series(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;

      using _BasicSum = __gnu_cxx::_BasicSum<_Val>;
      using _WenigerBasSum = __gnu_cxx::_WenigerDeltaSum<_BasicSum>;
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Val>;
      using _WenigerWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;
      using _WenigerSum = std::conditional_t<_Type == _StruveH,
					     _WenigerWijnSum, _WenigerBasSum>;
      int __sign = (_Type == _StruveH ? -1 : _Type == _StruveL ? +1 : 0);
      assert(__sign != 0);

      constexpr int _S_max_iter = 1000;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Real>::__root_pi;

      auto __x2 = __x / _Val{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Val{1};
      auto __struve = _WenigerSum(_Val{1});
      __struve += __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
      	  __term *= __xx4 / _Val(__k + 0.5L) / (__nu + _Val(__k + 0.5L));
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      auto __factor = _Val{2} * std::pow(__x2, __nu + _Val{1})
		    / std::__detail::__gamma(__nu + _Val{1.5L}) / _S_sqrt_pi;

      return __factor * __struve();
    }

  /**
   * Return either the Struve function of the second kind
   * @f$ \boldmath{K}_\nu(x) @f$ or the modified Struve function
   * of the second kind @f$ \boldmath{M}_\nu(x) @f$
   * depending on whether @c sign is +1 or -1 respectively.
   */
  template<_StruveType _Type, typename _Tp>
    _Tp
    __struve_asymp(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;

      using _BasicSum = __gnu_cxx::_BasicSum<_Val>;
      using _WenigerBasSum = __gnu_cxx::_WenigerDeltaSum<_BasicSum>;
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Val>;
      using _WenigerWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;
      using _WenigerSum = std::conditional_t<_Type == _StruveM,
					     _WenigerWijnSum, _WenigerBasSum>;

      int __sign = (_Type == _StruveK ? +1 : _Type == _StruveM ? -1 : 0);
      assert(__sign != 0);

      constexpr int _S_max_iter = 1000;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Real>::__root_pi;

      auto __x2 = __x / _Val{2};
      auto __xx4 = _Val(__sign) * __x2 * __x2;
      auto __term = _Val{1};
      auto __struve = _WenigerSum(_Val{1});
      __struve += __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
      	  __term *= _Val(__k - 0.5L) / (_Val(-__k - 0.5L) + __nu) / __xx4;
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      auto __fact = _Val(__sign) * std::pow(__x2, __nu - _Val{1})
		  / std::__detail::__gamma(__nu + _Val{0.5L}) / _S_sqrt_pi;

      return __fact * __struve();
    }

  /**
   * Return the Struve function of the first kind
   * @f$ \boldmath{H}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_h(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_max = _Real{__gnu_cxx::__digits10<_Real>()};

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_h: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_max)
	return __struve_series<_StruveH>(__nu, __x);
      else
	{
	  auto _Nnu = __cyl_neumann_n(__nu, __x);
	  auto _Knu = __struve_asymp<_StruveK>(__nu, __x);
	  return _Knu + _Nnu;
	}
    }

  /**
   * Return the Struve function of the second kind
   * @f$ \boldmath{K}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_k(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_max = _Real{__gnu_cxx::__digits10<_Real>()};

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_k: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) >= _S_max)
	return __struve_asymp<_StruveK>(__nu, __x);
      else
	{
	  auto _Nnu = __cyl_neumann_n(__nu, __x);
	  auto _Hnu = __struve_series<_StruveH>(__nu, __x);
	  return _Hnu - _Nnu;
	}
    }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_l(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_max = _Real{__gnu_cxx::__digits10<_Real>()};

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_l: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_max)
	return __struve_series<_StruveL>(__nu, __x);
      else
	{
	  auto _Inu = __cyl_bessel_i(__nu, __x);
	  auto _Mnu = __struve_asymp<_StruveM>(__nu, __x);
	  return _Mnu + _Inu;
	}
    }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_m(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_max = _Real{__gnu_cxx::__digits10<_Real>()};

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_k: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) >= _S_max)
	return __struve_asymp<_StruveM>(__nu, __x);
      else
	{
	  auto _Inu = __cyl_bessel_i(__nu, __x);
	  auto _Lnu = __struve_series<_StruveL>(__nu, __x);
	  return _Lnu - _Inu;
	}
    }

} // namespace __detail
} // namespace std

namespace __gnu_cxx
{

  // Struve functions (of the first kind)

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline float
  struve_hf(float __nu, float __x)
  { return std::__detail::__struve_h<float>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline long double
  struve_hl(long double __nu, long double __x)
  { return std::__detail::__struve_h<long double>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The Struve function is:
   * @f[
   *    \boldmath{H}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline __promote_fp_t<_Tpnu, _Tp>
    struve_h(_Tpnu __nu, _Tp __x)
    {
      using __type = __promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__struve_h<__type>(__nu, __x);
    }

  // Struve functions (of the second kind)

  /**
   * Return the Struve function of the first kind @f$ \boldmath{K}_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline float
  struve_kf(float __nu, float __x)
  { return std::__detail::__struve_k<float>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{K}_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline long double
  struve_kl(long double __nu, long double __x)
  { return std::__detail::__struve_k<long double>(__nu, __x); }

  /**
   * Return the Struve function of the second kind @f$ \boldmath{K}_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The Struve function of the second kind is:
   * @f[
   *    \boldmath{K}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline __promote_fp_t<_Tpnu, _Tp>
    struve_k(_Tpnu __nu, _Tp __x)
    {
      using __type = __promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__struve_k<__type>(__nu, __x);
    }

  // Modified Struve functions (of the first kind)

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ for @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline float
  struve_lf(float __nu, float __x)
  { return std::__detail::__struve_l<float>(__nu, __x); }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ for <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline long double
  struve_ll(long double __nu, long double __x)
  { return std::__detail::__struve_l<long double>(__nu, __x); }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The modified Struve function is:
   * @f[
   *    \boldmath{L}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline __promote_fp_t<_Tpnu, _Tp>
    struve_l(_Tpnu __nu, _Tp __x)
    {
      using __type = __promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__struve_l<__type>(__nu, __x);
    }

  // Modified Struve functions of the second kind

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ for @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline float
  struve_mf(float __nu, float __x)
  { return std::__detail::__struve_m<float>(__nu, __x); }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ for <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline long double
  struve_ml(long double __nu, long double __x)
  { return std::__detail::__struve_m<long double>(__nu, __x); }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The modified Struve function is:
   * @f[
   *    \boldmath{M}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline __promote_fp_t<_Tpnu, _Tp>
    struve_m(_Tpnu __nu, _Tp __x)
    {
      using __type = __promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__struve_m<__type>(__nu, __x);
    }

} // namespace __gnu_cxx


/**
 * Take a hard look at the series/asymptotic transition.
 */
template<typename _Tp>
  void
  test_struve_transition()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    for (int i = 500; i <= +5500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	std::cout << std::setw(width) << t;
	for (int n = 0; n <= 5; ++n)
	  {
	    auto nu = _Tp(1.0Q * n);
	    auto series = std::__detail::__struve_series<std::__detail::_StruveH>(nu, t);
	    auto asymp = std::__detail::__struve_asymp<std::__detail::_StruveK>(nu, t)
		       + std::__detail::__cyl_neumann_n(nu, t);
	    std::cout << '\t'
		      << std::setw(width) << series
		      << std::setw(width) << asymp
		      << std::setw(width) << asymp - series;
	  }
	std::cout << '\n';
      }
    std::cout << "\n\n";
  }


/**
 * Plot the Struve functions.
 */
template<typename _Tp>
  void
  plot_struve(std::string filename)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "H"
	 << std::setw(width) << "L"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = 0; i <= +3000; ++i)
      {
	auto t = _Tp(0.01Q * i);
	data << std::setw(width) << t;
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = _Tp(1.0Q * n);
	    data << '\t'
		 << std::setw(width) << __gnu_cxx::struve_h(nu, t)
		 << std::setw(width) << __gnu_cxx::struve_l(nu, t);
	  }
	data << '\n';
      }
    data << "\n\n";

    for (int i = 0; i <= +3000; ++i)
      {
	auto t = _Tp(0.01Q * i);
	data << std::setw(width) << t;
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = _Tp(1.0Q * n);
	    data << '\t'
		 << std::setw(width) << __gnu_cxx::struve_k(nu, t)
		 << std::setw(width) << __gnu_cxx::struve_m(nu, t);
	  }
	data << '\n';
      }
    data << "\n\n";
  }

void
test_struve()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  for (int i = 0; i <= +3000; ++i)
    {
      auto t = double(0.01Q * i);
      std::cout << std::setw(width) << t;
      for (int n = 0; n <= 20; ++n)
	{
	  auto nu = double(1.0Q * n);
	  auto h = __gnu_cxx::struve_h(nu, t);
	  auto l = __gnu_cxx::struve_l(nu, t);
	  auto hb = burkhardt::struve_h(nu, t);
	  auto lb = burkhardt::struve_l(nu, t);
	  std::cout << '\t'
	       << std::setw(width) << h
	       << std::setw(width) << hb
	       << std::setw(width) << (h - hb) / std::abs(hb)
	       << std::setw(width) << l
	       << std::setw(width) << lb
	       << std::setw(width) << (l - lb) / std::abs(lb);
	}
      std::cout << '\n';
    }
  std::cout << "\n\n";
}

int
main()
{
  test_struve_transition<float>();
  test_struve_transition<double>();
  test_struve_transition<long double>();

  //using cmplx = std::complex<double>;
  plot_struve<float>("plot/struve_float.txt");
  plot_struve<double>("plot/struve_double.txt");
  plot_struve<long double>("plot/struve_long_double.txt");

  test_struve();

  return 0;
}
