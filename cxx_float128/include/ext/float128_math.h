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

/** @file ext/float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_EXT_FLOAT128_MATH_H
#define _GLIBCXX_EXT_FLOAT128_MATH_H 1

#pragma GCC system_header

#ifdef _GLIBCXX_USE_FLOAT128

#include <quadmath.h>

extern "C++"
{
namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  inline __float128
  acos(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return acosq(__x); }

  inline __float128
  acosh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return acoshq(__x); }

  inline __float128
  asin(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return asinq(__x); }

  inline __float128
  asinh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return asinhq(__x); }

  inline __float128
  atan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return atanq(__x); }

  inline __float128
  atanh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return atanhq(__x); }

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
  erf(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return erfq(__x); }

  inline __float128
  erfc(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return erfcq(__x); }

  inline __float128
  exp(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return expq(__x); }

  inline __float128
  expm1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return expm1q(__x); }

  // libquadmath should get a real exp2!
  inline __float128
  exp2(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return powq(2, __x); }

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
  ilogb(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ilogbq(__x); }

  inline bool
  isinf(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return isinfq(__x); }

  inline bool
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

  inline __float128
  logb(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return logbq(__x); }

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

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std
} // extern "C++"

#endif // _GLIBCXX_USE_FLOAT128

#endif // _GLIBCXX_EXT_FLOAT128_MATH_H
