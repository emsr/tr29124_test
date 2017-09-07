// -*- C++ -*- header.

// Copyright (C) 2017 Free Software Foundation, Inc.
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

/** @file bits/complex128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_COMPLEX128_MATH_H
#define _GLIBCXX_BITS_COMPLEX128_MATH_H 1

#pragma GCC system_header

#if _GLIBCXX_HAVE_FLOAT128_MATH

typedef _Complex float __attribute__((mode(TC))) __complex128;

// FIXME: This is a fake quadmath.h
#ifdef __cplusplus
extern "C" {
#endif
  __float128 strtoflt128 (const char*, char**) _GLIBCXX_NOEXCEPT;

  __float128 cabsq (__complex128) _GLIBCXX_NOEXCEPT;
  __float128 cargq (__complex128) _GLIBCXX_NOEXCEPT;
  __float128 cimagq (__complex128) _GLIBCXX_NOEXCEPT;
  __float128 crealq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 cacosq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 cacoshq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 casinq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 casinhq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 catanq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 catanhq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 ccosq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 ccoshq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 cexpq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 cexpiq (__float128) _GLIBCXX_NOEXCEPT;
  __complex128 clogq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 clog10q (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 conjq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 cpowq (__complex128, __complex128) _GLIBCXX_NOEXCEPT;
  __complex128 cprojq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 csinq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 csinhq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 csqrtq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 ctanq (__complex128) _GLIBCXX_NOEXCEPT;
  __complex128 ctanhq (__complex128) _GLIBCXX_NOEXCEPT;
#ifdef __cplusplus
}
#endif

namespace std
{

#if _GLIBCXX_USE_C99_COMPLEX

  inline __float128
  __complex_abs(const __complex128& __z)
  { return cabsq(__z); }

  inline __float128
  __complex_arg(const __complex128& __z)
  { return cargq(__z); }

  inline __complex128
  __complex_cos(const __complex128& __z)
  { return ccosq(__z); }

  inline __complex128
  __complex_cosh(const __complex128& __z)
  { return ccoshq(__z); }

  inline __complex128
  __complex_exp(const __complex128& __z)
  { return cexpq(__z); }

  inline __complex128
  __complex_log(const __complex128& __z)
  { return clogq(__z); }

  inline __complex128
  __complex_sin(const __complex128& __z)
  { return csinq(__z); }

  inline __complex128
  __complex_sinh(const __complex128& __z)
  { return csinhq(__z); }      

  inline __complex128
  __complex_sqrt(const __complex128& __z)
  { return csqrtq(__z); }

  inline __complex128
  __complex_tan(const __complex128& __z)
  { return ctanq(__z); }

  inline __complex128
  __complex_tanh(const __complex128& __z)
  { return ctanhq(__z); }

  inline __complex128
  __complex_pow(const __complex128& __x, const __complex128& __y)
  { return cpowq(__x, __y); }

#if __cplusplus >= 201103L

  inline __complex128
  __complex_acos(const __complex128& __z)
  { return cacosq(__z); }

  inline __complex128
  __complex_asin(const __complex128& __z)
  { return casinq(__z); }

  inline __complex128
  __complex_atan(const __complex128& __z)
  { return catanq(__z); }

  inline __complex128
  __complex_acosh(const __complex128& __z)
  { return cacoshq(__z); }

  inline __complex128
  __complex_asinh(const __complex128& __z)
  { return casinhq(__z); }

  inline __complex128
  __complex_atanh(const __complex128& __z)
  { return catanhq(__z); }

  inline __complex128
  __complex_proj(const __complex128& __z)
  { return cprojq(__z); }

#endif // C++11

#endif // _GLIBCXX_USE_C99_COMPLEX

} // namespace std

#endif // _GLIBCXX_HAVE_FLOAT128_MATH

#endif // _GLIBCXX_BITS_COMPLEX128_MATH_H
