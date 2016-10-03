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

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the reperiodized sine of argument x:
   * @f[
   *   \sin_pi(x) = \sin(\pi x)
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
			? std::sin(__arg * _S_pi)
			: std::sin((_Tp{1} - __arg) * _S_pi);
	  return __sign * __sinval;
	}
    }

  // FIXME: Reperiodize the real part.
  template<typename _Tp>
    std::complex<_Tp>
    __sin_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr _Real _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      return std::sin(_S_pi * __z);
    }

  /**
   * Return the reperiodized cosine of argument x:
   * @f[
   *   \cos_pi(x) = \cos(\pi x)
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

  // FIXME: Reperiodize the real part.
  template<typename _Tp>
    std::complex<_Tp>
    __cos_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr _Real _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      return std::cos(_S_pi * __z);
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
      constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return -__tan_pi(-__x);
      else if (__x < _Tp{0.5L})
	return std::tan(__x * _S_pi);
      else
	return __tan_pi(__x - std::floor(__x));
    }

  // FIXME: Reperiodize the real part.
  template<typename _Tp>
    std::complex<_Tp>
    __tan_pi(std::complex<_Tp> __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr _Real _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      return std::tan(_S_pi * __z);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_TRIG_TCC
