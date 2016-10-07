// Special functions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/sf_theta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_THETA_TCC
#define _GLIBCXX_BITS_SF_THETA_TCC 1

#pragma GCC system_header

#include <vector>
#include <tuple>
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  Compute and return the @f$ \theta_1 @f$ function by series expansion.
   */
  template<typename _Tp>
    _Tp
    __theta_2_sum(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Val>();
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;
      auto __sum = std::exp(-__nu * __nu / __x);
      auto __sign = _Tp{-1};
      for (auto __k = 1; __k < 20; ++__k)
	{
	  auto __nup = __nu + _Tp(__k);
	  auto __termp = __sign * std::exp(-__nup * __nup / __x);
	  auto __num = __nu - _Tp(__k);
	  auto __termm = __sign * std::exp(-__num * __num / __x);
	  __sum += __termp + __termm;
	  __sign = -__sign;
	  if (std::abs(__termp) < _S_eps * std::abs(__sum)
	   && std::abs(__termm) < _S_eps * std::abs(__sum))
	    break;
	}
      return __sum / std::sqrt(_S_pi * __x);
    }

  /**
   *  Compute and return the @f$ \theta_3 @f$ function by series expansion.
   */
  template<typename _Tp>
    _Tp
    __theta_3_sum(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Val>();
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;
      auto __sum = std::exp(-__nu * __nu / __x);
      for (auto __k = 1; __k < 20; ++__k)
	{
	  auto __nup = __nu + _Tp(__k);
	  auto __termp = std::exp(-__nup * __nup / __x);
	  auto __num = __nu - _Tp(__k);
	  auto __termm = std::exp(-__num * __num / __x);
	  __sum += __termp + __termm;
	  if (std::abs(__termp) < _S_eps * std::abs(__sum)
	   && std::abs(__termm) < _S_eps * std::abs(__sum))
	    break;
	}
      return __sum / std::sqrt(_S_pi * __x);
    }

  /**
   *  Compute and return the @f$ \theta_2 @f$ function by series expansion.
   */
  template<typename _Tp>
    _Tp
    __theta_2_asymp(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Val>();
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;
      auto __sum = _Tp{0};
      for (auto __k = 0; __k < 20; ++__k)
	{
	  auto __thing = _Tp(2 * __k + 1) * _S_pi;
	  auto __cosarg = __nu * __thing;
	  auto __exparg = __thing * __thing * __x / _Tp{4};
	  auto __term = std::exp(-__exparg) * std::cos(__cosarg);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Tp{2} * __sum;
    }

  /**
   *  Compute and return the @f$ \theta_3 @f$ function by asymptotic series expansion.
   */
  template<typename _Tp>
    _Tp
    __theta_3_asymp(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Val>();
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;
      auto __sum = _Tp{0};
      for (auto __k = 1; __k < 20; ++__k)
	{
	  auto __thing = _Tp(2 * __k) * _S_pi;
	  auto __cosarg = __nu * __thing;
	  auto __exparg = __thing * __thing * __x / _Tp{4};
	  auto __term = std::exp(-__exparg) * std::cos(__cosarg);
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}
      return _Tp{1} + _Tp{2} * __sum;
    }

  /**
   * Return the exponential theta-2 function of period @c nu and argument @c x.
   *
   * The exponential theta-2 function is defined by
   * @f[
   *    \theta_2(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_2(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__x) <= _Tp{1} / _S_pi)
	return __theta_2_sum(__nu, __x);
      else
	return __theta_2_asymp(__nu, __x);
    }

  /**
   * Return the exponential theta-1 function of period @c nu and argument @c x.
   *
   * The Neville theta-1 function is defined by
   * @f[
   *    \theta_1(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j - 1/2)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_1(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else
	return __theta_2(__nu - _Tp{0.5L}, __x);
    }

  /**
   * Return the exponential theta-3 function of period @c nu and argument @c x.
   *
   * The exponential theta-3 function is defined by
   * @f[
   *    \theta_3(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    \exp\left( \frac{-(\nu+j)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 1) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_3(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__x) <= _Tp{1} / _S_pi)
	return __theta_3_sum(__nu, __x);
      else
	return __theta_3_asymp(__nu, __x);
    }

  /**
   * Return the exponential theta-2 function of period @c nu and argument @c x.
   *
   * The exponential theta-2 function is defined by
   * @f[
   *    \theta_2(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __theta_4(_Tp __nu, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Val>::__pi;

      if (__isnan(__nu) || __isnan(__x))
	return _S_NaN;
      else
	return __theta_3(__nu + _Tp{0.5L}, __x);
    }

  /**
   *  Use MacLaurin series to calculate the elliptic nome
   *  given the , k.
   */
  template<typename _Tp>
    _Tp
    __ellnome_series(_Tp __k)
    {
      auto __m = __k * __k; 
      return __m * ((_Tp{1} / _Tp{16})
	   + __m * ((_Tp{1} / _Tp{32})
	   + __m * ((_Tp{21} / _Tp{1024})
	   + __m * ((_Tp{31} / _Tp{2048})
	   + __m * (_Tp{6257} / _Tp{524288})))));
    }

  /**
   *  Use the arithmetic-geometric mean to calculate the elliptic nome
   *  given the , k.
   */
  template<typename _Tp>
    _Tp
    __ellnome_k(_Tp __k)
    {
      constexpr auto _S_pi = _Tp{3.1415926535897932384626433832795029L};
      auto __kp = std::sqrt((_Tp{1} - __k) * (_Tp{1} + __k));
      auto __K = __comp_ellint_1(__k);
      auto __Kp = __comp_ellint_1(__kp);
      return std::exp(-_S_pi * __Kp / __K);
    }

  /**
   * Return the elliptic nome given the modulus @c k.
   */
  template<typename _Tp>
    _Tp
    __ellnome(_Tp __k)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      if (__isnan(__k))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__ellnome:"
				      " argument k out of range"));
      else if (__k < std::pow(_Tp{67} * _S_eps, _Tp{0.125L}))
	return __ellnome_series(__k);
      else
	return __ellnome_k(__k);
    }

  /**
   * Return the Neville @f$ \theta_s @f$ function
   */
  template<typename _Tp>
    _Tp
    __theta_s(_Tp __k, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Val>::__pi_half;

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_s:"
				      " argument k out of range"));
      else
	{
	  auto __kc = std::sqrt((_Tp{1} - __k) * (_Tp{1} + __k));
	  auto _Kk = __comp_ellint_1(__k);
	  auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / (__k * __kc * _Kk))
	       * __theta_1(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_c @f$ function
   */
  template<typename _Tp>
    _Tp
    __theta_c(_Tp __k, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Val>::__pi_half;

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_c:"
				      " argument k out of range"));
      else
	{
	  auto _Kk = __comp_ellint_1(__k);
	  auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / (__k * _Kk))
	       * __theta_2(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_d @f$ function
   */
  template<typename _Tp>
    _Tp
    __theta_d(_Tp __k, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Val>::__pi_half;

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_d:"
				      " argument k out of range"));
      else
	{
	  auto _Kk = __comp_ellint_1(__k);
	  auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / _Kk)
	       * __theta_3(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return the Neville @f$ \theta_n @f$ function
   */
  template<typename _Tp>
    _Tp
    __theta_n(_Tp __k, _Tp __x)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Val>::__pi_half;

      if (__isnan(__k) || __isnan(__x))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__theta_n:"
				      " argument k out of range"));
      else
	{
	  auto __kc = std::sqrt((_Tp{1} - __k) * (_Tp{1} + __k));
	  auto _Kk = __comp_ellint_1(__k);
	  auto __q = __ellnome(__k);
	  return std::sqrt(_S_pi_2 / (__kc * _Kk))
	       * __theta_4(__q, _S_pi_2 * __x / _Kk);
	}
    }

  /**
   * Return a tuple of the three primary Jacobi elliptic functions:
   * @f$ sn(k, u), cn(k, u), dn(k, u) @f$.
   */
  template<typename _Tp>
    std::tuple<_Tp, _Tp, _Tp>
    __jacobi_sncndn(_Tp __k, _Tp __u)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Val>();
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Val>();

      if (__isnan(__k) || __isnan(__u))
	return std::make_tuple(_S_NaN, _S_NaN, _S_NaN);
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__jacobi_sncndn:"
				      " argument k out of range"));
      else if (std::abs(_Tp{1} - __k) < _Tp{2} * _S_eps)
	{
	  auto __sn = std::tanh(__u);
	  auto __cn = _Tp{1} / std::cosh(__u);
	  auto __dn = __cn;
	  return std::make_tuple(__sn, __cn, __dn);
	}
      else if (std::abs(__k) < _Tp{2} * _S_eps)
	{
	  auto __sn = std::sin(__u);
	  auto __cn = std::cos(__u);
	  auto __dn = _Tp{1};
	  return std::make_tuple(__sn, __cn, __dn);
	}
      else
	{
	  constexpr auto _S_CA = std::sqrt(_S_eps);
	  constexpr auto _S_N = 100;
	  std::vector<_Tp> __m;
	  std::vector<_Tp> __n;
	  __m.reserve(20);
	  __n.reserve(20);
	  _Tp __c, __d;
	  auto __mc = _Tp{1} - __k * __k;
	  bool __bo = (__mc < _Tp{0});
	  if (__bo)
	    {
	      __d = _Tp{1} - __mc;
	      __mc /= -_Tp{1} / __d;
	      __u *= (__d = std::sqrt(__d));
	    }
	  auto __a = _Tp{1};
	  auto __dn = _Tp{1};
	  auto __l = _S_N;
	  for (auto __i = 0; __i < _S_N; ++__i)
	    {
	      __l = __i;
	      __m.push_back(__a);
	      __n.push_back(__mc = std::sqrt(__mc));
	      __c = 0.5 * (__a + __mc);
	      if (std::abs(__a - __mc) <= _S_CA * __a)
		break;
	      __mc *= __a;
	      __a = __c;
	    }
	  __u *= __c;
	  auto __sn = std::sin(__u);
	  auto __cn = std::cos(__u);
	  if (__sn != _Tp{0})
	    {
	      __a = __cn / __sn;
	      __c *= __a;
	      for (auto __ii = __l; __ii + 1 >= 1; --__ii)
		{
		  _Tp __b = __m[__ii];
		  __a *= __c;
		  __c *= (__dn);
		  __dn = (__n[__ii] + __a) / (__b + __a);
		  __a = __c / __b;
		}
	      __a = _Tp{1} / std::hypot(_Tp{1}, __c);
	      __sn = std::copysign(__a, __sn);
	      __cn = __c * __sn;
	    }
	  if (__bo)
	    {
	      __a = __dn;
	      __dn = __cn;
	      __cn = __a;
	      __sn /= __d;
	    }
	  return std::make_tuple(__sn, __cn, __dn);
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_THETA_TCC
