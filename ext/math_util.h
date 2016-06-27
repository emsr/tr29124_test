// Math extensions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file ext/math_util.h
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_MATH_UTIL_H
#define _EXT_MATH_UTIL_H 1

#pragma GCC system_header

#if __cplusplus < 201103L
# include <bits/c++0x_warning.h>
#else

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A function to reliably compare two floating point numbers.
   *
   * @param __a The left hand side
   * @param __b The right hand side
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return returns true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    bool
    __fpequal(const _Tp& __a, const _Tp& __b, const _Tp __mul = _Tp{5})
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tol = __mul * _S_eps;
      bool __retval = true;
      if ((__a != _Tp{0}) || (__b != _Tp{0}))
	// Looks mean, but is necessary that the next line has sense.
	__retval = (std::abs(__a - __b) < std::max(std::abs(__a),
						   std::abs(__b)) * _S_tol);
      return __retval;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // C++11

#endif // _EXT_MATH_UTIL_H
