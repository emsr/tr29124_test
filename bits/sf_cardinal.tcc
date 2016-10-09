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

/** @file bits/sf_airy.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SF_CARDINAL_TCC
#define _GLIBCXX_BITS_SF_CARDINAL_TCC 1

#pragma GCC system_header

#include <bits/complex_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief Return the sinus cardinal function
   * @f[
   *   sinc(x) = \frac{\sin(x)}{x}
   * @f]
   */
  template<typename _Tp>
    __gnu_cxx::__promote_fp_t<_Tp>
    __sinc(_Tp __x)
    {
      if (__isnan(__x))
        return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (std::abs(__x) == __gnu_cxx::__infinity<_Tp>())
	return _Tp{0};
      else if (std::abs(__x) < __gnu_cxx::__sqrt_min<_Tp>())
        return _Tp{1} - __x * __x / _Tp{6};
      else
	return std::sin(__x) / __x;
    }

  /**
   * @brief Return the reperiodized sinus cardinal function
   * @f[
   *   sinc_\pi(x) = \frac{\sin(\pi x)}{\pi x}
   * @f]
   */
  template<typename _Tp>
    __gnu_cxx::__promote_fp_t<_Tp>
    __sinc_pi(_Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__isnan(__x))
        return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (std::abs(__x) == __gnu_cxx::__infinity<_Tp>())
	return _Tp{0};
      else
	{
	  auto __arg = _S_pi * __x;
	  if (std::abs(__arg) < _Tp{4} * __gnu_cxx::__sqrt_min<_Tp>())
	    return _Tp{1} - __arg * __arg / _Tp{6};
	  else
	    return __sin_pi(__x) / __arg;
	}
    }

  /**
   * @brief Return the hyperbolic sinus cardinal function
   * @f[
   *   sinhc(x) = \frac{\sinh(x)}{x}
   * @f]
   */
  template<typename _Tp>
    __gnu_cxx::__promote_fp_t<_Tp>
    __sinhc(_Tp __x)
    {
      if (__isnan(__x))
        return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (std::abs(__x) < _Tp{4} * __gnu_cxx::__sqrt_min<_Tp>())
        return _Tp{1} + __x * __x / _Tp{6};
      else
	return std::sinh(__x) / __x;
    }

  /**
   * @brief Return the reperiodized hyperbolic sinus cardinal function
   * @f[
   *   sinhc_\pi(x) = \frac{\sinh(\pi x)}{\pi x}
   * @f]
   */
  template<typename _Tp>
    __gnu_cxx::__promote_fp_t<_Tp>
    __sinhc_pi(_Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__isnan(__x))
        return __gnu_cxx::__quiet_NaN<_Tp>();
      else
	{
	  auto __arg = _S_pi * __x;
	  if (std::abs(__arg) < _Tp{4} * __gnu_cxx::__sqrt_min<_Tp>())
	    return _Tp{1} + __arg * __arg / _Tp{6};
	  else
	    return __sinh_pi(__x) / __arg;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_CARDINAL_TCC
