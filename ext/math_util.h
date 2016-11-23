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
   * A function to return the max of the absolute values of two numbers
   * ... so we won't include everything.
   */
  template<typename _Tp>
    inline _Tp
    __fpmaxabs(_Tp __a, _Tp __b)
    { std::abs(__a) < std::abs(__b) ? std::abs(__b) : std::abs(__a); }

  /**
   * A function to reliably compare two floating point numbers.
   *
   * @param __a The left hand side
   * @param __b The right hand side
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    inline bool
    __fp_is_equal(_Tp __a, _Tp __b, _Tp __mul = _Tp{1})
    {
      const auto _S_tol = __mul * std::numeric_limits<_Tp>::epsilon();
      bool __retval = true;
      if ((__a != _Tp{0}) || (__b != _Tp{0}))
	// Looks mean, but is necessary that the next line has sense.
	__retval = (std::abs(__a - __b) < __fpmaxabs(__a, __b) * _S_tol);
      return __retval;
    }

  /**
   * A function to reliably compare a floating point number with zero.
   *
   * @param __a The floating point number
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    inline bool
    __fp_is_zero(_Tp __a, _Tp __mul = _Tp{1})
    {
      const auto _S_tol = __mul * std::numeric_limits<_Tp>::epsilon();
      if (__a != _Tp{0})
	return (std::abs(__a) < _S_tol);
      else
        return true;
    }

  /**
   * A struct returned by floating point integer queries.
   */
  struct __fp_is_integer_t
  {
    // A flag indicating whether the floating point number is integralish.
    bool __is_integral;

    // An integer related to the floating point integral value.
    int __value;

    // Return __is_integral in a boolean context.
    operator bool() const
    { return this->__is_integral; }

    // Return __value.
    int
    operator()() const
    { return this->__value; }
  };

  /**
   * A function to reliably detect if a floating point number is an integer.
   *
   * @param __a The floating point number
   * @return @c true if a is an integer within mul * epsilon.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_integer(_Tp __a, _Tp __mul = _Tp{1})
    {
      auto __n = static_cast<int>(std::nearbyint(__a));
      return __fp_is_integer_t{__fp_is_equal(__a, _Tp(__n), __mul), __n};
    }

  /**
   * A function to reliably detect if a floating point number is a half-integer.
   *
   * @param __a The floating point number
   * @return @c true if 2*a is an integer within mul * epsilon.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_half_integer(_Tp __a, _Tp __mul = _Tp{1})
    {
      auto __n = static_cast<int>(std::nearbyint(_Tp{2} * __a));
      return __fp_is_integer_t{__fp_is_equal(_Tp{2} * __a, _Tp(__n), __mul), __n / 2};
    }

  /**
   * A function to reliably detect if a floating point number is a
   * half-odd-integer.
   *
   * @param __a The floating point number
   * @return @c true if 2*a is an odd integer within mul * epsilon.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_half_odd_integer(_Tp __a, _Tp __mul = _Tp{1})
    {
      auto __n = static_cast<int>(std::nearbyint(_Tp{2} * __a));
      bool __halfodd = (__n & 1 == 1)
		      && __fp_is_equal(_Tp{2} * __a, _Tp(__n), __mul);
      return __fp_is_integer_t{__halfodd, (__n - 1) / 2};
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // C++11

#endif // _EXT_MATH_UTIL_H
