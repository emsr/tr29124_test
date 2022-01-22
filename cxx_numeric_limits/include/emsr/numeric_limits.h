// math special functions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file bits/numeric_limits.h
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef NUMERIC_LIMITS_H
#define NUMERIC_LIMITS_H 1

#define _GLIBCXX_NO_SPECFUN 1
#include <cmath>
#undef _GLIBCXX_NO_SPECFUN

#include <limits>

namespace emsr
{

  /**
   *  @brief Part of std::numeric_limits.
   *  The idea is that types with, say, non-constexpr or even dynamic epsilon()
   *  can participate in this.
   *  I think variable templates could be specialized with non-constexpr types
   *  but I need something to work in C++11 and variable templates won't allow
   *  extraction of variable max from a mp number.
   */

  // Constexpr function template versions of std::numeric_limits.

  template<typename Tp>
    constexpr bool
    is_specialized(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_specialized; }

  template<typename Tp>
    constexpr Tp
    lim_min(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::min(); }

  template<typename Tp>
    constexpr Tp
    lim_max(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::max(); }

  template<typename Tp>
    constexpr Tp
    lowest(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::lowest(); }

  template<typename Tp>
    constexpr int
    digits(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::digits; }

  template<typename Tp>
    constexpr int
    digits10(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::digits10; }

  template<typename Tp>
    constexpr int
    max_digits10(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::max_digits10; }

  template<typename Tp>
    constexpr bool
    is_signed(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_signed; }

  template<typename Tp>
    constexpr bool
    is_integer(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_integer; }

  template<typename Tp>
    constexpr bool
    is_exact(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_exact; }

  template<typename Tp>
    constexpr int
    radix(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::radix; }

  template<typename Tp>
    constexpr Tp
    epsilon(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::epsilon(); }

  template<typename Tp>
    constexpr Tp
    round_error(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::round_error(); }

  template<typename Tp>
    constexpr int
    min_exponent(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::min_exponent; }

  template<typename Tp>
    constexpr int
    min_exponent10(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::min_exponent10; }

  template<typename Tp>
    constexpr int
    max_exponent(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::max_exponent; }

  template<typename Tp>
    constexpr int
    max_exponent10(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::max_exponent10; }

  template<typename Tp>
    constexpr bool
    has_infinity(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::has_infinity; }

  template<typename Tp>
    constexpr bool
    has_quiet_NaN(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::has_quiet_NaN; }

  template<typename Tp>
    constexpr bool
    has_signaling_NaN(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::has_signaling_NaN; }

  template<typename Tp>
    constexpr std::float_denorm_style
    has_denorm(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::has_denorm; }

  template<typename Tp>
    constexpr bool
    has_denorm_loss(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::has_denorm_loss; }

  template<typename Tp>
    constexpr Tp
    infinity(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::infinity(); }

  template<typename Tp>
    constexpr Tp
    quiet_NaN(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::quiet_NaN(); }

  template<typename Tp>
    constexpr Tp
    signaling_NaN(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::signaling_NaN(); }

  template<typename Tp>
    constexpr Tp
    denorm_min(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::denorm_min(); }

  template<typename Tp>
    constexpr bool
    is_iec559(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_iec559; }

  template<typename Tp>
    constexpr bool
    is_bounded(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_bounded; }

  template<typename Tp>
    constexpr bool
    is_modulo(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::is_modulo; }

  template<typename Tp>
    constexpr bool
    traps(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::traps; }

  template<typename Tp>
    constexpr bool
    tinyness_before(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::tinyness_before; }

  template<typename Tp>
    constexpr std::float_round_style
    round_style(Tp = Tp{}) noexcept
    { return std::numeric_limits<Tp>::round_style; }

  // Extra bits to help with numerics...
  // These depend on constexpr math functions.

  template<typename Tp>
    inline Tp
    max_integer(Tp = Tp{}) noexcept
    { return std::ldexp(Tp{1}, digits(Tp{})); }

  template<typename Tp>
    inline Tp
    sqrt_max(Tp = Tp{}) noexcept
    { return std::sqrt(lim_max(Tp{})); }

  template<typename Tp>
    inline Tp
    cbrt_max(Tp = Tp{}) noexcept
    { return std::cbrt(lim_max(Tp{})); }

  template<typename Tp>
    inline Tp
    root_max(Tp root) noexcept
    { return std::pow(lim_max(Tp{}), 1 / root); }

  template<typename Tp>
    inline Tp
    log_max(Tp = Tp{}) noexcept
    { return std::log(lim_max(Tp{})); }

  template<typename Tp>
    inline Tp
    log10_max(Tp = Tp{}) noexcept
    { return std::log10(lim_max(Tp{})); }


  template<typename Tp>
    inline Tp
    sqrt_min(Tp = Tp{}) noexcept
    { return std::sqrt(lim_min(Tp{})); }

  template<typename Tp>
    inline Tp
    cbrt_min(Tp = Tp{}) noexcept
    { return std::cbrt(lim_min(Tp{})); }

  template<typename Tp>
    inline Tp
    root_min(Tp root) noexcept
    { return std::pow(lim_min(Tp{}), 1 / root); }

  template<typename Tp>
    inline Tp
    log_min(Tp = Tp{}) noexcept
    { return std::log(lim_min(Tp{})); }

  template<typename Tp>
    inline Tp
    log10_min(Tp = Tp{}) noexcept
    { return std::log10(lim_min(Tp{})); }

  template<typename Tp>
    inline Tp
    sqrt_eps(Tp = Tp{}) noexcept
    { return std::sqrt(epsilon(Tp{})); }

  template<typename Tp>
    inline Tp
    cbrt_eps(Tp = Tp{}) noexcept
    { return std::cbrt(epsilon(Tp{})); }

  template<typename Tp>
    inline Tp
    root_eps(Tp root) noexcept
    { return std::pow(epsilon(Tp{}), 1 / root); }

  template<typename Tp>
    inline Tp
    log_eps(Tp = Tp{}) noexcept
    { return std::log(epsilon(Tp{})); }

  template<typename Tp>
    inline Tp
    log10_eps(Tp = Tp{}) noexcept
    { return std::log10(epsilon(Tp{})); }

} // namespace emsr

#endif // NUMERIC_LIMITS_H
