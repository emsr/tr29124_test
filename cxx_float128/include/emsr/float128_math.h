
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

/** @file emsr/float128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef FLOAT128_MATH_H
#define FLOAT128_MATH_H 1

#include <emsr/float128.h>

#ifdef EMSR_HAVE_FLOAT128

#include <quadmath.h>

extern "C++"
{
namespace std
{

  inline __float128
  acos(__float128 __x) noexcept
  { return acosq(__x); }

  inline __float128
  acosh(__float128 __x) noexcept
  { return acoshq(__x); }

  inline __float128
  asin(__float128 __x) noexcept
  { return asinq(__x); }

  inline __float128
  asinh(__float128 __x) noexcept
  { return asinhq(__x); }

  inline __float128
  atan(__float128 __x) noexcept
  { return atanq(__x); }

  inline __float128
  atanh(__float128 __x) noexcept
  { return atanhq(__x); }

  inline __float128
  atan2(__float128 __y, __float128 __x) noexcept
  { return atan2q(__y, __x); }

  inline __float128
  cbrt(__float128 __x) noexcept
  { return cbrtq(__x); }

  inline __float128
  ceil(__float128 __x) noexcept
  { return ceilq(__x); }

  inline __float128
  copysign(__float128 __x, __float128 __y) noexcept
  { return copysignq(__x, __y); }

  inline __float128
  cos(__float128 __x) noexcept
  { return cosq(__x); }

  inline __float128
  cosh(__float128 __x) noexcept
  { return coshq(__x); }

  inline __float128
  erf(__float128 __x) noexcept
  { return erfq(__x); }

  inline __float128
  erfc(__float128 __x) noexcept
  { return erfcq(__x); }

  inline __float128
  exp(__float128 __x) noexcept
  { return expq(__x); }

  inline __float128
  expm1(__float128 __x) noexcept
  { return expm1q(__x); }

  // libquadmath should get a real exp2!
  inline __float128
  exp2(__float128 __x) noexcept
  { return powq(2, __x); }

  inline __float128
  fabs(__float128 __x) noexcept
  { return fabsq(__x); }

  inline __float128
  fdim(__float128 __x, __float128 __y) noexcept
  { return fdimq(__x, __y); }

  inline __float128
  floor(__float128 __x) noexcept
  { return floorq(__x); }

  inline __float128
  fma(__float128 __m, __float128 __x, __float128 __b) noexcept
  { return fmaq(__m, __x, __b); }

  inline __float128
  fmax(__float128 __x, __float128 __y) noexcept
  { return fmaxq(__x, __y); }

  inline __float128
  fmin(__float128 __x, __float128 __y) noexcept
  { return fminq(__x, __y); }

  inline __float128
  fmod(__float128 __x, __float128 __y) noexcept
  { return fmodq(__x, __y); }

  inline __float128
  frexp(__float128 __x, int* __exp) noexcept
  { return frexpq(__x, __exp); }

  inline __float128
  hypot(__float128 __x, __float128 __y) noexcept
  { return hypotq(__x, __y); }

  inline __float128
  hypot(__float128 __x, __float128 __y, __float128 __z) noexcept
  { // FIXME
    return hypotq(hypotq(__x, __y), __z);
  }

  inline int
  ilogb(__float128 __x) noexcept
  { return ilogbq(__x); }

  inline bool
  isinf(__float128 __x) noexcept
  { return isinfq(__x); }

  inline bool
  isnan(__float128 __x) noexcept
  { return isnanq(__x); }

  inline __float128
  j0(__float128 __x) noexcept
  { return j0q(__x); }

  inline __float128
  j1(__float128 __x) noexcept
  { return j1q(__x); }

  inline __float128
  jn(int __n, __float128 __x) noexcept
  { return jnq(__n, __x); }

  inline __float128
  ldexp(__float128 __x, int __exp) noexcept
  { return ldexpq(__x, __exp); }

  inline __float128
  lgamma(__float128 __x) noexcept
  { return lgammaq(__x); }

  inline long long int
  llrint(__float128 __x) noexcept
  { return llrintq(__x); }

  inline long long int
  llround(__float128 __x) noexcept
  { return llroundq(__x); }

  inline __float128
  logb(__float128 __x) noexcept
  { return logbq(__x); }

  inline __float128
  log(__float128 __x) noexcept
  { return logq(__x); }

  inline __float128
  log10(__float128 __x) noexcept
  { return log10q(__x); }

  inline __float128
  log2(__float128 __x) noexcept
  { return log2q(__x); }

  inline __float128
  log1p(__float128 __x) noexcept
  { return log1pq(__x); }

  inline long int
  lrint(__float128 __x) noexcept
  { return lrintq(__x); }

  inline long int
  lround(__float128 __x) noexcept
  { return lroundq(__x); }

  inline __float128
  modf(__float128 __x, __float128* __iptr) noexcept
  { return modfq(__x, __iptr); }

  inline __float128
  nanq(const char* __str) noexcept
  { return __builtin_nanq(__str); }

  inline __float128
  nearbyint(__float128 __x) noexcept
  { return nearbyintq(__x); }

  inline __float128
  nextafter(__float128 __x, __float128 __y) noexcept
  { return nextafterq(__x, __y); }

  inline __float128
  pow(__float128 __x, __float128 __y) noexcept
  { return powq(__x, __y); }

  inline __float128
  remainder(__float128 __x, __float128 __y) noexcept
  { return remainderq(__x, __y); }

  inline __float128
  remquo(__float128 __x, __float128 __y, int* __n) noexcept
  { return remquoq(__x, __y, __n); }

  inline __float128
  rint(__float128 __x) noexcept
  { return rintq(__x); }

  inline __float128
  round(__float128 __x) noexcept
  { return roundq(__x); }

  inline __float128
  scalbln(__float128 __x, long int __n) noexcept
  { return scalblnq(__x, __n); }

  inline __float128
  scalbn(__float128 __x, int __n) noexcept
  { return scalbnq(__x, __n); }

  inline int
  signbit(__float128 __x) noexcept
  { return signbitq(__x); }

  inline void
  sincos(__float128 __x, __float128 * __sin, __float128 * __cos)
  noexcept
  { return sincosq(__x, __sin, __cos); }

  inline __float128
  sin(__float128 __x) noexcept
  { return sinq(__x); }

  inline __float128
  sinh(__float128 __x) noexcept
  { return sinhq(__x); }

  inline __float128
  sqrt(__float128 __x) noexcept
  { return sqrtq(__x); }

  inline __float128
  tan(__float128 __x) noexcept
  { return tanq(__x); }

  inline __float128
  tanh(__float128 __x) noexcept
  { return tanhq(__x); }

  inline __float128
  tgamma(__float128 __x) noexcept
  { return tgammaq(__x); }

  inline __float128
  trunc(__float128 __x) noexcept
  { return truncq(__x); }

  inline __float128
  y0(__float128 __x) noexcept
  { return y0q(__x); }

  inline __float128
  y1(__float128 __x) noexcept
  { return y1q(__x); }

  inline __float128
  yn(int __n, __float128 __x) noexcept
  { return ynq(__n, __x); }

} // namespace std
} // extern "C++"

#endif // EMSR_HAVE_FLOAT128

#endif // FLOAT128_MATH_H
