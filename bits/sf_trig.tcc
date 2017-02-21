// Special functions -*- C++ -*-

// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the reperiodized sine of argument x:
   * @f[
   *   \mathrm{sin_\pi}(x) = \sin(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __sin_pi(_Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
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
   *   \mathrm{sinh_\pi}(x) = \sinh(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __sinh_pi(_Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__x < _Tp{0})
	return -__sinh_pi(-__x);
      else
	std::sinh(_S_pi * __x);
    }

  /**
   * Return the reperiodized cosine of argument x:
   * @f[
   *   \mathrm{cos_\pi}(x) = \cos(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __cos_pi(_Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
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
   *   \mathrm{cosh_\pi}(x) = \cosh(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __cosh_pi(_Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
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
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      return std::tan(_S_pi * (__x - std::floor(__x)));
    }

  /**
   * Return the reperiodized hyperbolic tangent of argument x:
   * @f[
   *   \mathrm{tanh_\pi}(x) = \tanh(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __tanh_pi(_Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      return std::tanh(_S_pi * __x);
    }

  /**
   * Return the reperiodized sine of complex argument z:
   * @f[
   *   \mathrm{sin_\pi}(z) = \sin(\pi z)
   *     = \mathrm{sin_\pi}(x) \mathrm{cosh_\pi}(y)
   *   + i \mathrm{cos_\pi}(x) \mathrm{sinh_\pi}(y)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sin_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return __sin_pi(__x) * std::cosh(_S_pi * __y)
	   + _S_i * __cos_pi(__x) * std::sinh(_S_pi * __y);
    }

  /**
   * Return the reperiodized hyperbolic sine of complex argument z:
   * @f[
   *   \mathrm{sinh_\pi}(z) = \sinh(\pi z)
   *     = \mathrm{\sinh_\pi}(x) \mathrm{cos_\pi}(y)
   *   + i \mathrm{\cosh_\pi}(x) \mathrm{sin_\pi}(y)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sinh_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return std::sinh(_S_pi * __x) * __cos_pi(__y)
	+ _S_i * std::cosh(_S_pi * __x) * __sin_pi(__y);
    }

  /**
   * Return the reperiodized cosine of complex argument z:
   * @f[
   *    \mathrm{cos_\pi}(z) = \cos(\pi z)
   *       = \mathrm{cos_\pi}(x) \mathrm{cosh_\pi}(y)
   *     - i \mathrm{sin_\pi}(x) \mathrm{sinh_\pi}(y)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __cos_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return __cos_pi(__x) * std::cosh(_S_pi * __y)
	   - _S_i * __sin_pi(__x) * std::sinh(_S_pi * __y);
    }

  /**
   * Return the reperiodized hyperbolic cosine of complex argument z:
   * @f[
   *    \mathrm{cosh_\pi}(z) = \mathrm{cosh_\pi}(z)
   *       = \mathrm{cosh_\pi}(x) \mathrm{cos_\pi}(y)
   *     + i \mathrm{sinh_\pi}(x) \mathrm{sin_\pi}(y)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __cosh_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      return std::cosh(_S_pi * __x) * __cos_pi(__y)
	   + _S_i * std::sinh(_S_pi * __x) * __sin_pi(__y);
    }

  /**
   * Return the reperiodized tangent of complex argument z:
   * @f[
   *   \mathrm{tan_\pi}(z) = \tan(\pi z)
   *     = \frac{\mathrm{tan_\pi}(x) + i \mathrm{tanh_\pi}(y)}
   *            {1 - i \mathrm{tan_\pi}(x) \mathrm{tanh_\pi}(y)}
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __tan_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __x = std::real(__z);
      auto __y = std::imag(__z);
      auto __tan = __tan_pi(__x);
      auto __tanh = std::tanh(_S_pi * __y);
      return (__tan + _S_i * __tanh) / (1 - _S_i * __tan * __tanh);
    }

  /**
   * Return the reperiodized hyperbolic tangent of complex argument z:
   * @f[
   *   \mathrm{tanh_\pi}(z) = \tanh(\pi z)
   *     = \frac{\mathrm{tanh_\pi}(x) + i \mathrm{tan_\pi}(y)}
   *            {1 + i \mathrm{tanh_\pi}(x) \mathrm{tan_\pi}(y)}
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __tanh_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i = std::complex<_Tp>{0, 1};
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
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__x);
      if (std::isnan(__x))
	return __gnu_cxx::__sincos_t<_Tp>{_S_NaN, _S_NaN};
      else if (__x < _Tp{0})
	{
	  __gnu_cxx::__sincos_t<_Tp> __tempsc = __sincos_pi(-__x);
	  return __gnu_cxx::__sincos_t<_Tp>{-__tempsc.__sin_v,
					     __tempsc.__cos_v};
	}
      else if (__x < _Tp{0.5L})
	return __sincos(_S_pi * __x);
      else if (__x < _Tp{1})
	{
	  __gnu_cxx::__sincos_t<_Tp>
	    __tempsc = __sincos(_S_pi * (_Tp{1} - __x));
	  return __gnu_cxx::__sincos_t<_Tp>{__tempsc.__sin_v,
					   -__tempsc.__cos_v};
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
      return std::complex<_Tp>(__rho * __sc.__cos_v,
			       __rho * __sc.__sin_v);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_TRIG_TCC
