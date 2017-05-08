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

/** @file bits/limits_float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{limits}
 */

#ifndef _GLIBCXX_BITS_FLOAT128_LIMITS_H
#define _GLIBCXX_BITS_FLOAT128_LIMITS_H 1

#pragma GCC system_header

//#include <quadmath.h>
#define FLT128_MAX 1.18973149535723176508575932662800702e4932Q
#define FLT128_MIN 3.36210314311209350626267781732175260e-4932Q
#define FLT128_EPSILON 1.92592994438723585305597794258492732e-34Q
#define FLT128_DENORM_MIN 6.475175119438025110924438958227646552e-4966Q
#define FLT128_MANT_DIG 113
#define FLT128_MIN_EXP (-16381)
#define FLT128_MAX_EXP 16384
#define FLT128_DIG 33
#define FLT128_MIN_10_EXP (-4931)
#define FLT128_MAX_10_EXP 4932

#if _GLIBCXX_HAVE_FLOAT128_MATH

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /// numeric_limits<__float128> specialization.
  template<>
    struct numeric_limits<__float128>
    {
      static _GLIBCXX_USE_CONSTEXPR bool is_specialized = true;

      static _GLIBCXX_CONSTEXPR __float128 
      min() _GLIBCXX_USE_NOEXCEPT { return FLT128_MIN; }

      static _GLIBCXX_CONSTEXPR __float128 
      max() _GLIBCXX_USE_NOEXCEPT { return FLT128_MAX; }

#if __cplusplus >= 201103L
      static _GLIBCXX_CONSTEXPR __float128 
      lowest() _GLIBCXX_USE_NOEXCEPT { return -FLT128_MAX; }
#endif

      static _GLIBCXX_USE_CONSTEXPR int digits = FLT128_MANT_DIG;
      static _GLIBCXX_USE_CONSTEXPR int digits10 = FLT128_DIG;
#if __cplusplus >= 201103L
      static _GLIBCXX_USE_CONSTEXPR int max_digits10
	 = __glibcxx_max_digits10 (FLT128_MANT_DIG);
#endif
      static _GLIBCXX_USE_CONSTEXPR bool is_signed = true;
      static _GLIBCXX_USE_CONSTEXPR bool is_integer = false;
      static _GLIBCXX_USE_CONSTEXPR bool is_exact = false;
      static _GLIBCXX_USE_CONSTEXPR int radix = __FLT_RADIX__;

      static _GLIBCXX_CONSTEXPR __float128 
      epsilon() _GLIBCXX_USE_NOEXCEPT { return FLT128_EPSILON; }

      static _GLIBCXX_CONSTEXPR __float128 
      round_error() _GLIBCXX_USE_NOEXCEPT { return 0.5Q; }

      static _GLIBCXX_USE_CONSTEXPR int min_exponent = FLT128_MIN_EXP;
      static _GLIBCXX_USE_CONSTEXPR int min_exponent10 = FLT128_MIN_10_EXP;
      static _GLIBCXX_USE_CONSTEXPR int max_exponent = FLT128_MAX_EXP;
      static _GLIBCXX_USE_CONSTEXPR int max_exponent10 = FLT128_MAX_10_EXP;

      static _GLIBCXX_USE_CONSTEXPR bool has_infinity = true;
      static _GLIBCXX_USE_CONSTEXPR bool has_quiet_NaN = true;
      static _GLIBCXX_USE_CONSTEXPR bool has_signaling_NaN = has_quiet_NaN;
      static _GLIBCXX_USE_CONSTEXPR float_denorm_style has_denorm
	= denorm_present;
      static _GLIBCXX_USE_CONSTEXPR bool has_denorm_loss
	= __glibcxx_float128_has_denorm_loss;

      static _GLIBCXX_CONSTEXPR __float128
      infinity() _GLIBCXX_USE_NOEXCEPT { return __builtin_infq(); }

      static _GLIBCXX_CONSTEXPR __float128 
      quiet_NaN() _GLIBCXX_USE_NOEXCEPT { return __builtin_nanq(""); }

      static _GLIBCXX_CONSTEXPR __float128 
      signaling_NaN() _GLIBCXX_USE_NOEXCEPT { return __builtin_nansq(""); }

      static _GLIBCXX_CONSTEXPR __float128 
      denorm_min() _GLIBCXX_USE_NOEXCEPT { return FLT128_DENORM_MIN; }

      static _GLIBCXX_USE_CONSTEXPR bool is_iec559
	= has_infinity && has_quiet_NaN && has_denorm == denorm_present;
      static _GLIBCXX_USE_CONSTEXPR bool is_bounded = true;
      static _GLIBCXX_USE_CONSTEXPR bool is_modulo = false;

      static _GLIBCXX_USE_CONSTEXPR bool traps = __glibcxx_float128_traps;
      static _GLIBCXX_USE_CONSTEXPR bool tinyness_before
	= __glibcxx_float128_tinyness_before;
      static _GLIBCXX_USE_CONSTEXPR float_round_style round_style
	= round_to_nearest;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_HAVE_FLOAT128_MATH

#endif // _GLIBCXX_BITS_FLOAT128_LIMITS_H
