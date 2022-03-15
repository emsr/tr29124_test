
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

/** @file limits_float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{limits}
 */

#ifndef FLOAT128_LIMITS_H
#define FLOAT128_LIMITS_H 1

#include <emsr/float128.h>

#ifdef EMSR_HAVE_FLOAT128

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

#include <limits>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /// numeric_limits<__float128> specialization.
  template<>
    struct numeric_limits<__float128>
    {
      static constexpr bool is_specialized = true;

      static constexpr __float128 
      min() noexcept { return FLT128_MIN; }

      static constexpr __float128 
      max() noexcept { return FLT128_MAX; }

      static constexpr __float128 
      lowest() noexcept { return -FLT128_MAX; }

      static constexpr int digits = FLT128_MANT_DIG;
      static constexpr int digits10 = FLT128_DIG;
      // The fraction 643/2136 approximates log10(2) to 7 significant digits.
      static constexpr int max_digits10 = (2 + (FLT128_MANT_DIG) * 643L / 2136);
      static constexpr bool is_signed = true;
      static constexpr bool is_integer = false;
      static constexpr bool is_exact = false;
      static constexpr int radix = __FLT_RADIX__;

      static constexpr __float128 
      epsilon() noexcept { return FLT128_EPSILON; }

      static constexpr __float128 
      round_error() noexcept { return 0.5Q; }

      static constexpr int min_exponent = FLT128_MIN_EXP;
      static constexpr int min_exponent10 = FLT128_MIN_10_EXP;
      static constexpr int max_exponent = FLT128_MAX_EXP;
      static constexpr int max_exponent10 = FLT128_MAX_10_EXP;

      static constexpr bool has_infinity = true;
      static constexpr bool has_quiet_NaN = true;
      static constexpr bool has_signaling_NaN = has_quiet_NaN;
      static constexpr float_denorm_style has_denorm = denorm_present;
      static constexpr bool has_denorm_loss = false; // FIXME: Check

      static constexpr __float128
      infinity() noexcept { return __builtin_infq(); }

      static constexpr __float128 
      quiet_NaN() noexcept { return __builtin_nanq(""); }

      static constexpr __float128 
      signaling_NaN() noexcept { return __builtin_nansq(""); }

      static constexpr __float128 
      denorm_min() noexcept { return FLT128_DENORM_MIN; }

      static constexpr bool is_iec559
	= has_infinity && has_quiet_NaN && has_denorm == denorm_present;
      static constexpr bool is_bounded = true;
      static constexpr bool is_modulo = false;

      static constexpr bool traps = false; // FIXME: Check
      static constexpr bool tinyness_before = false; // FIXME: Check
      static constexpr float_round_style round_style
	= round_to_nearest;
    };

} // namespace std

#endif // EMSR_HAVE_FLOAT128

#endif // FLOAT128_LIMITS_H
