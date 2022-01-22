// -*- C++ -*- header.

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file ext/numeric_limits_float128.h
 */

#ifndef NUMERIC_LIMITS_FLOAT128_H
#define NUMERIC_LIMITS_FLOAT128_H 1

#define _GLIBCXX_USE_FLOAT128 0
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
#  if __has_include(<quadmath.h>)
#    include <quadmath.h>
#    define _GLIBCXX_USE_FLOAT128 1
#  endif
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#ifdef _GLIBCXX_USE_FLOAT128

#include <ext/float128_limits.h>
#include <ext/float128_math.h>
#include <emsr/numeric_limits.h>

namespace emsr
{

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
    constexpr bool
    is_specialized<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_specialized; }

  template<>
    constexpr __float128
    lim_min<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::min(); }

  template<>
    constexpr __float128
    lim_max<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::max(); }

  template<>
    constexpr __float128
    lowest<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::lowest(); }

  template<>
    constexpr int
    digits<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::digits; }

  template<>
    constexpr int
    digits10<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::digits10; }

  template<>
    constexpr int
    max_digits10<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::max_digits10; }

  template<>
    constexpr bool
    is_signed<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_signed; }

  template<>
    constexpr bool
    is_integer<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_integer; }

  template<>
    constexpr bool
    is_exact<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_exact; }

  template<>
    constexpr int
    radix<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::radix; }

  template<>
    constexpr __float128
    epsilon<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::epsilon(); }

  template<>
    constexpr __float128
    round_error<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::round_error(); }

  template<>
    constexpr int
    min_exponent<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::min_exponent; }

  template<>
    constexpr int
    min_exponent10<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::min_exponent10; }

  template<>
    constexpr int
    max_exponent<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::max_exponent; }

  template<>
    constexpr int
    max_exponent10<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::max_exponent10; }

  template<>
    constexpr bool
    has_infinity<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::has_infinity; }

  template<>
    constexpr bool
    has_quiet_NaN<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::has_quiet_NaN; }

  template<>
    constexpr bool
    has_signaling_NaN<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::has_signaling_NaN; }

  template<>
    constexpr std::float_denorm_style
    has_denorm<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::has_denorm; }

  template<>
    constexpr bool
    has_denorm_loss<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::has_denorm_loss; }

  template<>
    constexpr __float128
    infinity<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::infinity(); }

  template<>
    constexpr __float128
    quiet_NaN<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::quiet_NaN(); }

  template<>
    constexpr __float128
    signaling_NaN<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::signaling_NaN(); }

  template<>
    constexpr __float128
    denorm_min<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::denorm_min(); }

  template<>
    constexpr bool
    is_iec559<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_iec559; }

  template<>
    constexpr bool
    is_bounded<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_bounded; }

  template<>
    constexpr bool
    is_modulo<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::is_modulo; }

  template<>
    constexpr bool
    traps<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::traps; }

  template<>
    constexpr bool
    tinyness_before<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::tinyness_before; }

  template<>
    constexpr std::float_round_style
    round_style<__float128>(__float128) noexcept
    { return std::numeric_limits<__float128>::round_style; }

  // Extra bits to help with numerics...
  // These depend math functions which aren't constexpr for __float128.
  // These are specializations of the functions in ext/numeric_limits.h

  template<>
    inline __float128
    max_integer(__float128) noexcept
    { return std::ldexp(1, digits(__float128{})); }

  template<>
    inline __float128
    sqrt_max<__float128>(__float128) noexcept
    { return std::sqrt(lim_max(__float128{})); }

  template<>
    inline __float128
    cbrt_max<__float128>(__float128) noexcept
    { return std::cbrt(lim_max(__float128{})); }

  template<>
    inline __float128
    root_max(__float128 __root) noexcept
    { return std::pow(lim_max(__float128{}), 1 / __root); }

  template<>
    inline __float128
    log_max<__float128>(__float128) noexcept
    { return std::log(lim_max(__float128{})); }

  template<>
    inline __float128
    log10_max<__float128>(__float128) noexcept
    { return std::log10(lim_max(__float128{})); }


  template<>
    inline __float128
    sqrt_min<__float128>(__float128) noexcept
    { return std::sqrt(lim_min(__float128{})); }

  template<>
    inline __float128
    cbrt_min<__float128>(__float128) noexcept
    { return std::cbrt(lim_min(__float128{})); }

  template<>
    inline __float128
    root_min(__float128 __root) noexcept
    { return std::pow(lim_min(__float128{}), 1 / __root); }

  template<>
    inline __float128
    log_min<__float128>(__float128) noexcept
    { return std::log(lim_min(__float128{})); }

  template<>
    inline __float128
    log10_min<__float128>(__float128) noexcept
    { return std::log10(lim_min(__float128{})); }

  template<>
    inline __float128
    sqrt_eps<__float128>(__float128) noexcept
    { return std::sqrt(epsilon(__float128{})); }

  template<>
    inline __float128
    cbrt_eps<__float128>(__float128) noexcept
    { return std::cbrt(epsilon(__float128{})); }

  template<>
    inline __float128
    root_eps(__float128 __root) noexcept
    { return std::pow(epsilon(__float128{}), 1 / __root); }

  template<>
    inline __float128
    log_eps<__float128>(__float128) noexcept
    { return std::log(epsilon(__float128{})); }

  template<>
    inline __float128
    log10_eps<__float128>(__float128) noexcept
    { return std::log10(epsilon(__float128{})); }

} // namespace emsr

#endif // _GLIBCXX_USE_FLOAT128

#endif // NUMERIC_LIMITS_FLOAT128_H
