// TR29124 math special functions -*- C++ -*-

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

/** @file bits/sf_trig.tcc
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SF_TRIG_TCC
#define _GLIBCXX_BITS_SF_TRIG_TCC 1

#pragma GCC system_header

#include <bits/complex_util.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A struct to store a cosine and a sine value.
   */
  template<typename _Tp>
    struct __sincos_t
    {
      _Tp sin_value;
      _Tp cos_value;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the reperiodized sine of argument x:
   * @f[
   *   \sin_\pi(x) = \sin(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __sin_pi(_Tp __x)
    {
      constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return -__sin_pi(-__x);
      else if (__x < _Tp{0.5L})
	return std::sin(__x * _S_pi);
      else if (__x < _Tp{1})
	return std::sin((_Tp{1} - __x) * _S_pi);
      else
	{
	  auto __nu = std::floor(__x);
	  auto __arg = __x - __nu;
	  auto __sign = (int(__nu) & 1) == 1 ? -1 : +1;
	  auto __sinval = (__arg < _Tp{0.5L})
			? __sin_pi(__arg)
			: __sin_pi(_Tp{1} - __arg);
	  return __sign * __sinval;
	}
    }

  /**
   * Return the reperiodized hyperbolic sine of argument x:
   * @f[
   *   \sinh_\pi(x) = \sinh(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __sinh_pi(_Tp __x)
    {
      constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return -__sinh_pi(-__x);
      else
	std::sinh(_S_pi * __x);
    }

  /**
   * Return the reperiodized cosine of argument x:
   * @f[
   *   \cos_\pi(x) = \cos(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __cos_pi(_Tp __x)
    {
      constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return __cos_pi(-__x);
      else if (__x < _Tp{0.5L})
	return std::cos(__x * _S_pi);
      else if (__x < _Tp{1})
	return -std::cos((_Tp{1} - __x) * _S_pi);
      else
	{
	  auto __nu = std::floor(__x);
	  auto __arg = __x - __nu;
	  auto __sign = (int(__nu) & 1) == 1 ? -1 : +1;
	  return __sign * __cos_pi(__arg);
	}
    }

  /**
   * Return the reperiodized hyperbolic cosine of argument x:
   * @f[
   *   \cosh_\pi(x) = \cosh(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __cosh_pi(_Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return __cosh_pi(-__x);
      else
	return std::cosh(_S_pi * __x);
    }

  /**
   * Return the reperiodized tangent of argument x:
   * @f[
   *   \tan_pi(x) = \tan(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __tan_pi(_Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      return std::tan(_S_pi * (__x - std::floor(__x)));
    }

  /**
   * Return the reperiodized hyperbolic tangent of argument x:
   * @f[
   *   \tanh_\pi(x) = \tanh(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __tanh_pi(_Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      return std::tanh(_S_pi * __x);
    }

  /**
   * Return the reperiodized sine of complex argument z:
   * @f[
   *   \sin_\pi(z) = \sin(\pi z)
   *     = \sin_\pi(x) \cosh_\pi(y) + i \cos_\pi(x) \sinh_\pi(y)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sin_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return __sin_pi(__x) * std::cosh(_S_pi * __y)
	   + _S_i * __cos_pi(__x) * std::sinh(_S_pi * __y);
    }

  /**
   * Return the reperiodized hyperbolic sine of complex argument z:
   * @f[
   *   \sinh_\pi(z) = \sinh(\pi z)
   *     = \sinh(\pi x) \cos_\pi(y) + i \cosh(\pi x) \sin_\pi(y)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sinh_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return std::sinh(_S_pi * __x) * __cos_pi(__y)
	+ _S_i * std::cosh(_S_pi * __x) * __sin_pi(__y);
    }

  /**
   * Return the reperiodized cosine of complex argument z:
   * \cos_\pi(z) = \cos(\pi z)
   *    = \cos_\pi(x) \cosh_\pi(y) - i \sin_\pi(x) \sinh_\pi(y)
   */
  template<typename _Tp>
    std::complex<_Tp>
    __cos_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return __cos_pi(__x) * std::cosh(_S_pi * __y)
	   - _S_i * __sin_pi(__x) * std::sinh(_S_pi * __y);
    }

  /**
   * Return the reperiodized hyperbolic cosine of complex argument z:
   * \cosh_\pi(z) = \cosh_\pi(z)
   *    = \cosh_\pi(x) \cos_\pi(y) + i \sinh_\pi(x) \sin_\pi(y)
   */
  template<typename _Tp>
    std::complex<_Tp>
    __cosh_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return std::cosh(_S_pi * __x) * __cos_pi(__y)
	   + _S_i * std::sinh(_S_pi * __x) * __sin_pi(__y);
    }

  /**
   * Return the reperiodized tangent of complex argument z:
   * @f[
   *   \tan_\pi(z) = \tan(\pi z)
   *     = \frac{\tan_\pi(x) + i \tanh_\pi(y)}{1 - i \tan_\pi(x) \tanh_\pi(y)}
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __tan_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      auto __tan = __tan_pi(__x);
      auto __tanh = std::tanh(_S_pi * __y);
      return (__tan + _S_i * __tanh) / (1 - _S_i * __tan * __tanh);
    }

  /**
   * Return the reperiodized hyperbolic tangent of complex argument z:
   * @f[
   *   \tanh_\pi(z) = \tanh(\pi z)
   *     = \frac{\tanh_\pi(x) + i \tan_\pi(y)}{1 + i \tanh_\pi(x) \tan_\pi(y)}
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __tanh_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      auto __tanh = std::tanh(_S_pi * __x);
      auto __tan = __tan_pi(__y);
      return (__tanh + _S_i * __tan) / (1 + _S_i * __tanh * __tan);
    }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__sincos_t<_Tp>
    __sincos(_Tp __x)
    { return __gnu_cxx::__sincos_t<_Tp>{std::sin(__x), std::cos(__x)}; }

  /**
   * 
   */
  template<>
    inline __gnu_cxx::__sincos_t<float>
    __sincos(float __x)
    {
      float __sin, __cos;
      __builtin_sincosf(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<float>{__sin, __cos};
    }

  /**
   * 
   */
  template<>
    inline __gnu_cxx::__sincos_t<double>
    __sincos(double __x)
    {
      double __sin, __cos;
      __builtin_sincos(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<double>{__sin, __cos};
    }

  /**
   * 
   */
  template<>
    inline __gnu_cxx::__sincos_t<long double>
    __sincos(long double __x)
    {
      long double __sin, __cos;
      __builtin_sincosl(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<long double>{__sin, __cos};
    }

  /**
   * Reperiodized sincos.
   */
  template<typename _Tp>
    __gnu_cxx::__sincos_t<_Tp>
    __sincos_pi(_Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
      if (std::isnan(__x))
	return __gnu_cxx::__sincos_t<_Tp>{_S_NaN, _S_NaN};
      else if (__x < _Tp{0})
	{
	  __gnu_cxx::__sincos_t<_Tp> __tempsc = __sincos_pi(-__x);
	  return __gnu_cxx::__sincos_t<_Tp>{-__tempsc.sin_value,
					     __tempsc.cos_value};
	}
      else if (__x < _Tp{0.5L})
	return __sincos(_S_pi * __x);
      else if (__x < _Tp{1})
	{
	  __gnu_cxx::__sincos_t<_Tp>
	    __tempsc = __sincos(_S_pi * (_Tp{1} - __x));
	  return __gnu_cxx::__sincos_t<_Tp>{__tempsc.sin_value,
					   -__tempsc.cos_value};
	}
      else
	{
	  auto __nu = std::floor(__x);
	  auto __arg = __x - __nu;
	  auto __sign = (int(__nu) & 1) == 1 ? _Tp{-1} : _Tp{+1};

	  auto __sinval = (__arg < _Tp{0.5L})
			? std::sin(_S_pi * __arg)
			: std::sin(_S_pi * (_Tp{1} - __arg));
	  auto __cosval = std::cos(_S_pi * __arg);
	  return __gnu_cxx::__sincos_t<_Tp>{__sign * __sinval,
					    __sign * __cosval};
	}
    }

  /**
   * Reperiodized complex constructor.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polar_pi(_Tp __rho, _Tp __phi_pi)
    {
      __gnu_cxx::__sincos_t<_Tp> __sc = __sincos_pi(__phi_pi);
      return std::complex<_Tp>(__rho * __sc.cos_value, __rho * __sc.sin_value);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_TRIG_TCC
