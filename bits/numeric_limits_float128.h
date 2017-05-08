// -*- C++ -*- header.

// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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

/** @file bits/numeric_limits_float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_NUMERIC_LIMITS_FLOAT128_H
#define _GLIBCXX_BITS_NUMERIC_LIMITS_FLOAT128_H 1

#pragma GCC system_header

#define _GLIBCXX_HAVE_FLOAT128_MATH 0
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
#  if __has_include(<quadmath.h>)
#    include <quadmath.h>
#    define _GLIBCXX_HAVE_FLOAT128_MATH 1
#  endif
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#if _GLIBCXX_HAVE_FLOAT128_MATH

#include <bits/float128_math.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief An extension/wrapper of std::numeric_limits.
   *
   *  The idea is that types with, say, non-constexpr or even dynamic epsilon()
   *  can use this generic API.
   *  I think variable templates could be specialized with non-constexpr types
   *  but I need something to work in C++11 and variable templates won't allow
   *  extraction of variable max, epsilon, etc. from a multi-precision number.
   */

  // Constexpr function template versions of std::numeric_limits.

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_specialized<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_specialized; }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::min(); }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::max(); }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __lowest<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::lowest(); }

  template<>
    _GLIBCXX_CONSTEXPR int
    __digits<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::digits; }

  template<>
    _GLIBCXX_CONSTEXPR int
    __digits10<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::digits10; }

  template<>
    _GLIBCXX_CONSTEXPR int
    __max_digits10<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::max_digits10; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_signed<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_signed; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_integer<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_integer; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_exact<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_exact; }

  template<>
    _GLIBCXX_CONSTEXPR int
    __radix<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::radix; }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __epsilon<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::epsilon(); }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __round_error<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::round_error(); }

  template<>
    _GLIBCXX_CONSTEXPR int
    __min_exponent<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::min_exponent; }

  template<>
    _GLIBCXX_CONSTEXPR int
    __min_exponent10<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::min_exponent10; }

  template<>
    _GLIBCXX_CONSTEXPR int
    __max_exponent<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::max_exponent; }

  template<>
    _GLIBCXX_CONSTEXPR int
    __max_exponent10<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::max_exponent10; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __has_infinity<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::has_infinity; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __has_quiet_NaN<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::has_quiet_NaN; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __has_signaling_NaN<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::has_signaling_NaN; }

  template<>
    _GLIBCXX_CONSTEXPR std::float_denorm_style
    __has_denorm<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::has_denorm; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __has_denorm_loss<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::has_denorm_loss; }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __infinity<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::infinity(); }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __quiet_NaN<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::quiet_NaN(); }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __signaling_NaN<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::signaling_NaN(); }

  template<>
    _GLIBCXX_CONSTEXPR __float128
    __denorm_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::denorm_min(); }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_iec559<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_iec559; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_bounded<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_bounded; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __is_modulo<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::is_modulo; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __traps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::traps; }

  template<>
    _GLIBCXX_CONSTEXPR bool
    __tinyness_before<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::tinyness_before; }

  template<>
    _GLIBCXX_CONSTEXPR std::float_round_style
    __round_style<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<__float128>::round_style; }

  // Extra bits to help with numerics...
  // These depend math functions which aren't constexpr for __float128.
  // These are specializations of the functions in bits/numeric_limits.h

  template<>
    inline __float128
    __max_integer(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::ldexp(1, __digits(__float128{})); }

  template<>
    inline __float128
    __sqrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__max(__float128{})); }

  template<>
    inline __float128
    __cbrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__max(__float128{})); }

  template<>
    inline __float128
    __root_max(__float128 __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(__float128{}), 1 / __root); }

  template<>
    inline __float128
    __log_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__max(__float128{})); }

  template<>
    inline __float128
    __log10_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__max(__float128{})); }


  template<>
    inline __float128
    __sqrt_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__min(__float128{})); }

  template<>
    inline __float128
    __cbrt_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__min(__float128{})); }

  template<>
    inline __float128
    __root_min(__float128 __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(__float128{}), 1 / __root); }

  template<>
    inline __float128
    __log_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__min(__float128{})); }

  template<>
    inline __float128
    __log10_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__min(__float128{})); }

  template<>
    inline __float128
    __sqrt_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__epsilon(__float128{})); }

  template<>
    inline __float128
    __cbrt_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__epsilon(__float128{})); }

  template<>
    inline __float128
    __root_eps(__float128 __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(__float128{}), 1 / __root); }

  template<>
    inline __float128
    __log_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__epsilon(__float128{})); }

  template<>
    inline __float128
    __log10_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__epsilon(__float128{})); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_HAVE_FLOAT128_MATH

#endif // _GLIBCXX_BITS_NUMERIC_LIMITS_FLOAT128_H
