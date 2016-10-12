// Math extensions -*- C++ -*-

// Copyright (C) 2013-2016 Free Software Foundation, Inc.
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

/** @file ext/math_const_mpreal.h
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_MATH_CONST_MPREAL_H
#define _EXT_MATH_CONST_MPREAL_H 1

#pragma GCC system_header

#if __cplusplus < 201103L
# include <bits/c++0x_warning.h>
#else

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION


  // The following functions mirror the constants abobe but also
  // admit generic programming with non-findamental types.
  // For fundamental types, these constexpr functions return
  // the appropriate constant above.
  // Developers of multi-precision types are encouraged to overload
  // these functions calling multiprecision "constant" functions as
  // available.

  template<>
    mpfr::mpreal
    __const_4_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 4 * mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_4_pi_div_3(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 4 * mpfr::const_pi(__proto.getPrecision()) / 3; }

  template<>
    mpfr::mpreal
    __const_2_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 2 * mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_pi_half(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_pi(__proto.getPrecision()) / 2; }

  template<>
    mpfr::mpreal
    __const_pi_third(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_pi(__proto.getPrecision()) / 3; }

  template<>
    mpfr::mpreal
    __const_pi_quarter(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_pi(__proto.getPrecision()) / 4; }

  template<>
    mpfr::mpreal
    __const_root_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::sqrt(mpfr::const_pi(__proto.getPrecision())); }

  template<>
    mpfr::mpreal
    __const_cbrt_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::cbrt(mpfr::const_pi(__proto.getPrecision())); }

  template<>
    mpfr::mpreal
    __const_root_pi_div_2(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::sqrt(mpfr::const_pi(__proto.getPrecision()) / 2); }

  template<>
    mpfr::mpreal
    __const_one_div_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 1 / mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_two_div_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 2 / mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_two_div_root_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 2 / mpfr::sqrt(mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_pi_sqr_div_6(mpfr::mpreal __proto = mpfr::mpreal{})
    {
      auto __pi = mpfr::const_pi(__proto.getPrecision());
      return __pi * __pi / 6;
    }

  template<>
    mpfr::mpreal
    __const_deg(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 180 / mpfr::const_pi(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_rad(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_pi(__proto.getPrecision()) / 180; }

  template<>
    mpfr::mpreal
    __const_e(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::exp(1, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_one_div_e(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 1 / mpfr::exp(1, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_log2_e(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::log2(mpfr::exp(1, __proto.getPrecision())); }

  template<>
    mpfr::mpreal
    __const_log10_e(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::log10(mpfr::exp(1, __proto.getPrecision())); }

  template<>
    mpfr::mpreal
    __const_ln_2(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_log2(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_ln_3(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::log(3, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_ln_10(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::log(10, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_ln_pi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::log(mpfr::const_pi(__proto.getPrecision())); }

  template<>
    mpfr::mpreal
    __const_gamma_e(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_euler(__proto.getPrecision());

  template<>
    mpfr::mpreal
    __const_phi(mpfr::mpreal __proto = mpfr::mpreal{})
    { return (1 + mpfr::sqrt(5, __proto.getPrecision())) / 2; }

  template<>
    mpfr::mpreal
    __const_catalan(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::const_catalan(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_root_2(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::sqrt(2, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_root_3(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::sqrt(3, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_root_5(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::sqrt(5, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_root_7(mpfr::mpreal __proto = mpfr::mpreal{})
    { return mpfr::sqrt(7, __proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __const_one_div_root_2(mpfr::mpreal __proto = mpfr::mpreal{})
    { return 1 / mpfr::sqrt(2, __proto.getPrecision()); }


_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // C++11

#endif // _EXT_MATH_CONST_MPREAL_H
