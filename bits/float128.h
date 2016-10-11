// -*- C++ -*- header.

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

/** @file bits/float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_FLOAT128_H
#define _GLIBCXX_BITS_FLOAT128_H 1

#pragma GCC system_header

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

#include <limits>
#include <iosfwd>
#include <quadmath.h>

// From <limits>
#define __glibcxx_max_digits10(T) \
  (2 + (T) * 643L / 2136)

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  inline __float128
  acos(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return acosq(__x); }

  inline __float128
  asin(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return asinq(__x); }

  inline __float128
  atan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return atanq(__x); }

  inline __float128
  atan2(__float128 __y, __float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return atan2q(__y, __x); }

  inline __float128
  cbrt(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return cbrtq(__x); }

  inline __float128
  ceil(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ceilq(__x); }

  inline __float128
  copysign(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return copysignq(__x, __y); }

  inline __float128
  cos(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return cosq(__x); }

  inline __float128
  cosh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return coshq(__x); }

  inline __float128
  exp(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return expq(__x); }

  inline __float128
  erf(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return erfq(__x); }

  inline __float128
  erfc(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return erfcq(__x); }

  inline __float128
  expm1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return expm1q(__x); }

  inline __float128
  fabs(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return fabsq(__x); }

  inline __float128
  fdim(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fdimq(__x, __y); }

  inline __float128
  floor(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return floorq(__x); }

  inline __float128
  fma(__float128 __m, __float128 __x, __float128 __b) _GLIBCXX_USE_NOEXCEPT
  { return fmaq(__m, __x, __b); }

  inline __float128
  fmax(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fmaxq(__x, __y); }

  inline __float128
  fmin(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fminq(__x, __y); }

  inline __float128
  fmod(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fmodq(__x, __y); }

  inline __float128
  frexp(__float128 __x, int* __exp) _GLIBCXX_USE_NOEXCEPT
  { return frexpq(__x, __exp); }

  inline __float128
  hypot(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return hypotq(__x, __y); }

  inline int
  isinf(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return isinfq(__x); }

  inline int
  ilogb(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ilogbq(__x); }

  inline int
  isnan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return isnanq(__x); }

  inline __float128
  j0(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return j0q(__x); }

  inline __float128
  j1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return j1q(__x); }

  inline __float128
  jn(int __n, __float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return jnq(__n, __x); }

  inline __float128
  ldexp(__float128 __x, int __exp) _GLIBCXX_USE_NOEXCEPT
  { return ldexpq(__x, __exp); }

  inline __float128
  lgamma(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return lgammaq(__x); }

  inline long long int
  llrint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return llrintq(__x); }

  inline long long int
  llround(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return llroundq(__x); }

#ifndef NO_LOGBQ
  inline __float128
  logb(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return logbq(__x); }
#endif

  inline __float128
  log(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return logq(__x); }

  inline __float128
  log10(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return log10q(__x); }

  inline __float128
  log2(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return log2q(__x); }

  inline __float128
  log1p(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return log1pq(__x); }

  inline long int
  lrint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return lrintq(__x); }

  inline long int
  lround(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return lroundq(__x); }

  inline __float128
  modf(__float128 __x, __float128* __iptr) _GLIBCXX_USE_NOEXCEPT
  { return modfq(__x, __iptr); }

  inline __float128
  nanq(const char* __str) _GLIBCXX_USE_NOEXCEPT
  { return __builtin_nanq(__str); }

  inline __float128
  nearbyint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return nearbyintq(__x); }

  inline __float128
  nextafter(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return nextafterq(__x, __y); }

  inline __float128
  pow(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return powq(__x, __y); }

  inline __float128
  remainder(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return remainderq(__x, __y); }

  inline __float128
  remquo(__float128 __x, __float128 __y, int* __n) _GLIBCXX_USE_NOEXCEPT
  { return remquoq(__x, __y, __n); }

  inline __float128
  rint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return rintq(__x); }

  inline __float128
  round(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return roundq(__x); }

  inline __float128
  scalbln(__float128 __x, long int __n) _GLIBCXX_USE_NOEXCEPT
  { return scalblnq(__x, __n); }

  inline __float128
  scalbn(__float128 __x, int __n) _GLIBCXX_USE_NOEXCEPT
  { return scalbnq(__x, __n); }

  inline int
  signbit(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return signbitq(__x); }

  inline void
  sincos(__float128 __x, __float128 * __sin, __float128 * __cos)
  _GLIBCXX_USE_NOEXCEPT
  { return sincosq(__x, __sin, __cos); }

  inline __float128
  sin(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return sinq(__x); }

  inline __float128
  sinh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return sinhq(__x); }

  inline __float128
  sqrt(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return sqrtq(__x); }

  inline __float128
  tan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return tanq(__x); }

  inline __float128
  tanh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return tanhq(__x); }

  inline __float128
  tgamma(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return tgammaq(__x); }

  inline __float128
  trunc(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return truncq(__x); }

  inline __float128
  y0(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return y0q(__x); }

  inline __float128
  y1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return y1q(__x); }

  inline __float128
  yn(int __n, __float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ynq(__n, __x); }


  inline __float128
  mod(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return __x * __x; }

  inline __float128
  arg(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return __x < 0.0Q ? M_PIq : 0.0Q; }

  inline _GLIBCXX_USE_CONSTEXPR __float128
  imag(__float128) _GLIBCXX_USE_NOEXCEPT
  { return 0.0Q; }

  inline _GLIBCXX_USE_CONSTEXPR __float128
  real(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return __x; }


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
      static _GLIBCXX_USE_CONSTEXPR bool has_denorm_loss = true;

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

      static _GLIBCXX_USE_CONSTEXPR bool traps = false;//???
      static _GLIBCXX_USE_CONSTEXPR bool tinyness_before = 
					 false;//???
      static _GLIBCXX_USE_CONSTEXPR float_round_style round_style = 
						      round_to_nearest;
    };

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       __float128 __x);

  template<typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& __is, __float128& __x);

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

// From <limits>
#undef __glibcxx_max_digits10

#include <bits/numeric_limits.h>

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
    __float128
    __sqrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__max(__float128{})); }

#ifdef NO_CBRT
  template<>
    __float128
    __cbrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(__float128{}), 1 / 3.0Q); }
#else
  template<>
    __float128
    __cbrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__max(__float128{})); }
#endif

  template<>
    __float128
    __root_max(__float128 __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(__float128{}), 1 / __root); }

  template<>
    __float128
    __log_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__max(__float128{})); }

  template<>
    __float128
    __log10_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__max(__float128{})); }


  template<>
    __float128
    __sqrt_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__min(__float128{})); }

#ifdef NO_CBRT
  template<>
    __float128
    __cbrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(__float128{}), 1 / 3.0Q); }
#else
  template<>
    __float128
    __cbrt_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__min(__float128{})); }
#endif

  template<>
    __float128
    __root_min(__float128 __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(__float128{}), 1 / __root); }

  template<>
    __float128
    __log_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__min(__float128{})); }

  template<>
    __float128
    __log10_min<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__min(__float128{})); }

  template<>
    __float128
    __sqrt_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__epsilon(__float128{})); }

#ifdef NO_CBRT
  template<>
    __float128
    __cbrt_max<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(__float128{}), 1 / 3.0Q); }
#else
  template<>
    __float128
    __cbrt_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__epsilon(__float128{})); }
#endif

  template<>
    __float128
    __root_eps(__float128 __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(__float128{}), 1 / __root); }

  template<>
    __float128
    __log_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__epsilon(__float128{})); }

  template<>
    __float128
    __log10_eps<__float128>(__float128) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__epsilon(__float128{})); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#include <bits/float128.tcc>

#endif // _GLIBCXX_BITS_FLOAT128_H
