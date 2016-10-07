// math special functions -*- C++ -*-

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

/** @file bits/numeric_limits.h
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_NUMERIC_LIMITS_H
#define _GLIBCXX_BITS_NUMERIC_LIMITS_H 1

#pragma GCC system_header

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief Part of std::numeric_limits.
   *  The idea is that types with, say, non-constexpr or even dynamic epsilon()
   *  can participate in this.
   *  I think variable templates could be specialized with non-constexpr types
   *  but I need something to work in C++11 and variable templates won't allow
   *  extraction of variable max from a mp number.
   */

  // Constexpr function template versions of std::numeric_limits.

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_specialized(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_specialized; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::min(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __max(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::max(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __lowest(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::lowest(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __digits(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::digits; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __digits10(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::digits10; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __max_digits10(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::max_digits10; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_signed(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_signed; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_integer(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_integer; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_exact(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_exact; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __radix(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::radix; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __epsilon(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::epsilon(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __round_error(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::round_error(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __min_exponent(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::min_exponent; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __min_exponent10(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::min_exponent10; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __max_exponent(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::max_exponent; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR int
    __max_exponent10(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::max_exponent10; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __has_infinity(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::has_infinity; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __has_quiet_NaN(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::has_quiet_NaN; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __has_signaling_NaN(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::has_signaling_NaN; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR std::float_denorm_style
    __has_denorm(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::has_denorm; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __has_denorm_loss(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::has_denorm_loss; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __infinity(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::infinity(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __quiet_NaN(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::quiet_NaN(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __signaling_NaN(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::signaling_NaN(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR _Tp
    __denorm_min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::denorm_min(); }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_iec559(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_iec559; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_bounded(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_bounded; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __is_modulo(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::is_modulo; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __traps(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::traps; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR bool
    __tinyness_before(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::tinyness_before; }

  template<typename _Tp>
    _GLIBCXX_CONSTEXPR std::float_round_style
    __round_style(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<_Tp>::round_style; }

  // Extra bits to help with numerics...
  // These depend on constexpr math functions.

  template<typename _Tp>
    _Tp
    __sqrt_max(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__max(_Tp{})); }

#ifdef NO_CBRT
  template<typename _Tp>
    _Tp
    __cbrt_max(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(_Tp{}), 1 / _Tp{3}); }
#else
  template<typename _Tp>
    _Tp
    __cbrt_max(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__max(_Tp{})); }
#endif

  template<typename _Tp>
    _Tp
    __root_max(_Tp __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(_Tp{}), 1 / __root); }

  template<typename _Tp>
    _Tp
    __log_max(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__max(_Tp{})); }

  template<typename _Tp>
    _Tp
    __log10_max(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__max(_Tp{})); }


  template<typename _Tp>
    _Tp
    __sqrt_min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__min(_Tp{})); }

#ifdef NO_CBRT
  template<typename _Tp>
    _Tp
    __cbrt_min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(_Tp{}), 1 / _Tp{3}); }
#else
  template<typename _Tp>
    _Tp
    __cbrt_min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__min(_Tp{})); }
#endif

  template<typename _Tp>
    _Tp
    __root_min(_Tp __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(_Tp{}), 1 / __root); }

  template<typename _Tp>
    _Tp
    __log_min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__min(_Tp{})); }

  template<typename _Tp>
    _Tp
    __log10_min(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__min(_Tp{})); }

  template<typename _Tp>
    _Tp
    __sqrt_eps(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__epsilon(_Tp{})); }

#ifdef NO_CBRT
  template<typename _Tp>
    _Tp
    __cbrt_eps(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(_Tp{}), 1 / _Tp{3}); }
#else
  template<typename _Tp>
    _Tp
    __cbrt_eps(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__epsilon(_Tp{})); }
#endif

  template<typename _Tp>
    _Tp
    __root_eps(_Tp __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(_Tp{}), 1 / __root); }

  template<typename _Tp>
    _Tp
    __log_eps(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__epsilon(_Tp{})); }

  template<typename _Tp>
    _Tp
    __log10_eps(_Tp = _Tp{}) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__epsilon(_Tp{})); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_NUMERIC_LIMITS_H
