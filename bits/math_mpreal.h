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

/** @file bits/math_mpreal.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_MATH_MPREAL_H
#define _GLIBCXX_BITS_MATH_MPREAL_H 1

#pragma GCC system_header


namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  inline mpfr::mpreal
  acos(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::acos(__x); }

  inline mpfr::mpreal
  asin(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::asin(__x); }

  inline mpfr::mpreal
  atan(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::atan(__x); }

  inline mpfr::mpreal
  atan2(const mpfr::mpreal& __y, const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::atan2(__y, __x); }

  inline mpfr::mpreal
  cbrt(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::cbrt(__x); }

  inline mpfr::mpreal
  ceil(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::ceil(__x); }

  inline mpfr::mpreal
  copysign(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::copysign(__x, __y); }

  inline mpfr::mpreal
  cos(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::cos(__x); }

  inline mpfr::mpreal
  cosh(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::cosh(__x); }

  inline mpfr::mpreal
  exp(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::exp(__x); }

  inline mpfr::mpreal
  erf(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::erf(__x); }

  inline mpfr::mpreal
  erfc(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::erfc(__x); }

  inline mpfr::mpreal
  expm1(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::expm1(__x); }

  inline mpfr::mpreal
  fabs(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::fabs(__x); }

  inline mpfr::mpreal
  abs(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::fabs(__x); }

  inline mpfr::mpreal
  fdim(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::dim(__x, __y); }

  // This is what mpreal has
  inline mpfr::mpreal
  dim(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::dim(__x, __y); }

  inline mpfr::mpreal
  floor(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::floor(__x); }

  inline mpfr::mpreal
  fma(const mpfr::mpreal& __m, const mpfr::mpreal& __x, const mpfr::mpreal& __b) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::fma(__m, __x, __b); }

  inline mpfr::mpreal
  fmax(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::fmax(__x, __y); }

  inline mpfr::mpreal
  fmin(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::fmin(__x, __y); }

  inline mpfr::mpreal
  fmod(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::fmod(__x, __y); }

  inline mpfr::mpreal
  frexp(const mpfr::mpreal& __x, int* __exp) _GLIBCXX_USE_NOEXCEPT
  {
    mp_exp_t __mpexp;
    auto __ret = mpfr::frexp(__x, &__mpexp);
    *__exp = static_cast<int>(__mpexp);
    return __ret;
  }

  inline mpfr::mpreal
  hypot(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::hypot(__x, __y); }

  inline mpfr::mpreal
  hypot(const mpfr::mpreal& __x, const mpfr::mpreal& __y, const mpfr::mpreal& __z) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::hypot(__x, __y, __z); }

  inline int
  isinf(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::isinf(__x); }

  //inline int
  //ilogb(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  //{ return mpfr::ilogb(__x); }

  inline int
  isnan(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::isnan(__x); }

  inline mpfr::mpreal
  j0(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::besselj0(__x); }

  inline mpfr::mpreal
  j1(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::besselj1(__x); }

  inline mpfr::mpreal
  jn(int __n, const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::besseljn(__n, __x); }

  inline mpfr::mpreal
  ldexp(const mpfr::mpreal& __x, int __exp) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::ldexp(__x, __exp); }

  inline mpfr::mpreal
  lgamma(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::lgamma(__x); }

  inline long long int
  llrint(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::rint(__x).toLLong(); }

  inline long long int
  llround(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::round(__x).toLLong(); }

#ifndef NO_LOGBQ
  inline mpfr::mpreal
  logb(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::logb(__x); }
#endif

  inline mpfr::mpreal
  log(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::log(__x); }

  inline mpfr::mpreal
  log10(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::log10(__x); }

  inline mpfr::mpreal
  log2(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::log2(__x); }

  inline mpfr::mpreal
  log1p(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::log1p(__x); }

  inline long int
  lrint(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::rint(__x).toLong(); }

  inline long int
  lround(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::round(__x).toLong(); }

  inline mpfr::mpreal
  modf(const mpfr::mpreal& __x, mpfr::mpreal* __iptr) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::modf(__x, *__iptr); }

  //inline mpfr::mpreal
  //nan(const char*) _GLIBCXX_USE_NOEXCEPT
  //{ return std::numeric_limits<mpfr::mpreal>::quiet_NaN(); }

  inline mpfr::mpreal
  nearbyint(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  {
    mp_rnd_t __rnd = mpfr::mpreal::get_default_rnd();

    switch (__rnd)
    {
    case GMP_RNDN:
      return mpfr::rint_round(__x);
    case GMP_RNDZ:
      return mpfr::rint_trunc(__x);
    case GMP_RNDU:
      return mpfr::rint_ceil(__x);
    case GMP_RNDD:
      return mpfr::rint_floor(__x);
    default:
      return mpfr::rint_round(__x);
    }
  }

  inline mpfr::mpreal
  nextafter(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::nexttoward(__x, __y); }

  inline mpfr::mpreal
  pow(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::pow(__x, __y); }

  inline mpfr::mpreal
  remainder(const mpfr::mpreal& __x, const mpfr::mpreal& __y) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::remainder(__x, __y); }

  inline mpfr::mpreal
  remquo(const mpfr::mpreal& __x, const mpfr::mpreal& __y, int* __n) _GLIBCXX_USE_NOEXCEPT
  {
    long __ln = *__n;
    auto __ret = mpfr::remquo(&__ln, __x, __y);
    *__n = static_cast<int>(__ln);
    return __ret;
  }

  inline mpfr::mpreal
  rint(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return static_cast<int>(mpfr::rint(__x).toLong()); }

  inline mpfr::mpreal
  round(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return static_cast<int>(mpfr::round(__x).toLong()); }

  inline mpfr::mpreal
  scalbln(const mpfr::mpreal& __x, long int __n) _GLIBCXX_USE_NOEXCEPT
  {
    auto __exp = mp_exp_t(__n);
    return mpfr::scalbn(__x, __exp);
  }

  inline mpfr::mpreal
  scalbn(const mpfr::mpreal& __x, int __n) _GLIBCXX_USE_NOEXCEPT
  {
    auto __exp = mp_exp_t(__n);
    return mpfr::scalbn(__x, __exp);
  }

  inline int
  signbit(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::signbit(__x); }

  inline void
  sincos(const mpfr::mpreal& __x, mpfr::mpreal* __sin, mpfr::mpreal* __cos)
  _GLIBCXX_USE_NOEXCEPT
  {
    *__sin = mpfr::sin(__x);
    *__cos = mpfr::cos(__x);
  }

  inline mpfr::mpreal
  sin(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::sin(__x); }

  inline mpfr::mpreal
  sinh(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::sinh(__x); }

  inline mpfr::mpreal
  sqrt(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::sqrt(__x); }

  inline mpfr::mpreal
  tan(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::tan(__x); }

  inline mpfr::mpreal
  tanh(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::tanh(__x); }

  inline mpfr::mpreal
  tgamma(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::tgamma(__x); }

  inline mpfr::mpreal
  trunc(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::trunc(__x); }

  inline mpfr::mpreal
  y0(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::bessely0(__x); }

  inline mpfr::mpreal
  y1(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::bessely1(__x); }

  inline mpfr::mpreal
  yn(int __n, const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return mpfr::besselyn(__n, __x); }

  inline mpfr::mpreal
  mod(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return __x * __x; }

  inline mpfr::mpreal
  arg(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  {
    auto __prec = __x.getPrecision();
    return __x < 0 ? mpfr::const_pi(__prec) : mpfr::mpreal{0, __prec};
  }

  inline mpfr::mpreal
  imag(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  {
    auto __prec = __x.getPrecision();
    return mpfr::mpreal{0, __prec};
  }

  inline mpfr::mpreal
  real(const mpfr::mpreal& __x) _GLIBCXX_USE_NOEXCEPT
  { return __x; }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_MATH_MPREAL_H
