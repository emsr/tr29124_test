// $HOME/bin_specfun/bin/g++ -std=gnu++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -g -Wall -Wextra -o test_struve_new test_struve_new.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_struve_new > test_struve_new.new

// g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -Wall -Wextra -o test_struve_new test_struve_new.cpp

// ./test_struve_new > test_struve_new.txt

#include <cmath> // There are issues with <complex> inclusion if this isn't uphere!
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

namespace std
{
namespace __detail
{

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __struve_series(_Tp __nu, _Tp __x, int __sign)
    {
      using _Val = std::__detail::__num_traits_t<_Tp>;

      constexpr auto _S_eps = std::numeric_limits<_Val>::epsilon();
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Val>::__root_pi;

      auto __x2 = __x / _Tp{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Tp{1};
      auto __struve = __term;
      for (int __k = 1; __k < 100; ++__k)
	{
      	  __term *= __xx4 / _Tp(__k - 1 + 1.5L) / (__nu + _Tp(__k - 1 + 1.5L));
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      __struve *= _Tp{2} * std::pow(__x2, __nu + _Tp{1})
		/ std::__detail::__gamma(__nu + _Tp{1.5L}) / _S_sqrt_pi;

      return __struve;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __struve_asymp(_Tp __nu, _Tp __x, int __sign)
    {
      using _Val = std::__detail::__num_traits_t<_Tp>;

      constexpr auto _S_eps = std::numeric_limits<_Val>::epsilon();
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Val>::__root_pi;

      auto __x2 = __x / _Tp{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Tp{1};
      auto __struve = __term;
      for (int __k = 1; __k < 100; ++__k)
	{
      	  __term *= _Tp(__k - 1 + 0.5L) / (_Tp(-__k - 1 + 0.5L) + __nu) / __xx4;
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      __struve *= _Tp(__sign) * std::pow(__x2, __nu - _Tp{1})
		/ std::__detail::__gamma(__nu + _Tp{0.5L}) / _S_sqrt_pi;

      return __struve;
    }

  template<typename _Tp>
    _Tp
    __struve_h(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      constexpr auto _S_max = _Tp{20};
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__struve_h: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_max)
	return __struve_series(__nu, __x, -1);
      else
	{
	  auto _Nnu = __cyl_neumann_n(__nu, __x);
	  auto _Knu = __struve_asymp(__nu, __x, +1);
	  return _Knu + _Nnu;
	}
    }

  template<typename _Tp>
    _Tp
    __struve_k(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      constexpr auto _S_max = _Tp{20};
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__struve_k: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) >= _S_max)
	return __struve_asymp(__nu, __x, +1);
      else
	{
	  auto _Nnu = __cyl_neumann_n(__nu, __x);
	  auto _Hnu = __struve_series(__nu, __x, -1);
	  return _Hnu - _Nnu;
	}
    }

  template<typename _Tp>
    _Tp
    __struve_l(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      constexpr auto _S_max = _Tp{20};
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__struve_l: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_max)
	return __struve_series(__nu, __x, +1);
      else
	{
	  auto _Inu = __cyl_bessel_i(__nu, __x);
	  auto _Mnu = __struve_asymp(__nu, __x, -1);
	  return _Mnu + _Inu;
	}
    }

  template<typename _Tp>
    _Tp
    __struve_m(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      constexpr auto _S_max = _Tp{20};
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__struve_k: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return _S_nan;
      else if (std::abs(__x) >= _S_max)
	return __struve_asymp(__nu, __x, -1);
      else
	{
	  auto _Inu = __cyl_bessel_i(__nu, __x);
	  auto _Lnu = __struve_series(__nu, __x, +1);
	  return _Lnu - _Inu;
	}
    }

} // namespace __detail
} // namespace std

namespace __gnu_cxx
{

  // Struve functions (of the first kind)

  /**
   * Return the Struve function of the first kind @f$ H_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline float
  struve_hf(float __nu, float __x)
  { return std::__detail::__struve_h<float>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ H_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline long double
  struve_hl(long double __nu, long double __x)
  { return std::__detail::__struve_h<long double>(__nu, __x); }

  /**
   * Return the Struve function @f$ H_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The Struve function is:
   * @f[
   *    H_{\nu}(x) = \sum_{k=0}^{\infty}
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
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    struve_h(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return std::__detail::__struve_h<__type>(__nu, __x);
    }

  // Struve functions (of the second kind)

  /**
   * Return the Struve function of the first kind @f$ K_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline float
  struve_kf(float __nu, float __x)
  { return std::__detail::__struve_k<float>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ K_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline long double
  struve_kl(long double __nu, long double __x)
  { return std::__detail::__struve_k<long double>(__nu, __x); }

  /**
   * Return the Struve function @f$ K_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The Struve function is:
   * @f[
   *    K_{\nu}(x) = \sum_{k=0}^{\infty}
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
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    struve_k(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return std::__detail::__struve_k<__type>(__nu, __x);
    }

  // Modified Struve functions (of the first kind)

  /**
   * Return the modified Struve function of the first kind @f$ L_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline float
  struve_lf(float __nu, float __x)
  { return std::__detail::__struve_l<float>(__nu, __x); }

  /**
   * Return the modified Struve function of the first kind @f$ L_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline long double
  struve_ll(long double __nu, long double __x)
  { return std::__detail::__struve_l<long double>(__nu, __x); }

  /**
   * Return the modified Struve function @f$ L_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The modified Struve function is:
   * @f[
   *    L_{\nu}(x) = \sum_{k=0}^{\infty}
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
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    struve_l(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return std::__detail::__struve_l<__type>(__nu, __x);
    }

  // Modified Struve functions (of the second kind)

  /**
   * Return the modified Struve function of the first kind @f$ M_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline float
  struve_mf(float __nu, float __x)
  { return std::__detail::__struve_m<float>(__nu, __x); }

  /**
   * Return the modified Struve function of the first kind @f$ M_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline long double
  struve_ml(long double __nu, long double __x)
  { return std::__detail::__struve_m<long double>(__nu, __x); }

  /**
   * Return the modified Struve function @f$ M_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The modified Struve function is:
   * @f[
   *    M_{\nu}(x) = \sum_{k=0}^{\infty}
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
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    struve_m(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return std::__detail::__struve_m<__type>(__nu, __x);
    }

} // namespace __gnu_cxx

/**
 * 
 */
template<typename _Tp>
  void
  plot_struve(std::string filename)
  {
    using _Val = std::__detail::__num_traits_t<_Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Val>::digits10);
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
  }

int
main()
{
  //using cmplx = std::complex<double>;
  plot_struve<double>("plot/struve_double.txt");
  plot_struve<long double>("plot/struve_long_double.txt");

  return 0;
}
