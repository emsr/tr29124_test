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

/** @file ext/math_const.h
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_MATH_CONST_H
#define _EXT_MATH_CONST_H 1

#pragma GCC system_header

#if __cplusplus < 201103L
# include <bits/c++0x_warning.h>
#else

#include <type_traits>
#include <limits>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // A class for math constants.
  template<typename _RealType>
    struct __math_constants
    {
      static_assert(std::is_floating_point<_RealType>::value,
		    "template argument not a floating point type");

      // Constant @f$ 4\pi @f$.
      static _GLIBCXX_CONSTEXPR _RealType __4_pi
      	= 12.56637061435917295385057353311801153678L;
      // Constant @f$ 4\pi/3 @f$.
      static _GLIBCXX_CONSTEXPR _RealType __4_pi_div_3
      	= 4.188790204786390984616857844372670512253L;
      // Constant @f$ 2\pi @f$.
      static _GLIBCXX_CONSTEXPR _RealType __2_pi
      	= 6.283185307179586476925286766559005768391L;
      // Constant @f$ \pi @f$.
      static _GLIBCXX_CONSTEXPR _RealType __pi
      	= 3.141592653589793238462643383279502884195L;
      // Constant @f$ \pi / 2 @f$.
      static _GLIBCXX_CONSTEXPR _RealType __pi_half
      	= 1.570796326794896619231321691639751442098L;
      // Constant @f$ \pi / 3 @f$.
      static _GLIBCXX_CONSTEXPR _RealType __pi_third
      	= 1.047197551196597746154214461093167628063L;
      // Constant @f$ \pi / 4 @f$.
      static _GLIBCXX_CONSTEXPR _RealType __pi_quarter
      	= 0.785398163397448309615660845819875721049L;
      // Constant @f$ \sqrt(\pi / 2) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_pi
      	= 1.772453850905516027298167483341145182798L;
      static _GLIBCXX_CONSTEXPR _RealType __cbrt_pi
      	= 1.464591887561523263020142527263790391736L;
      // Constant @f$ \sqrt(\pi / 2) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_pi_div_2
      	= 1.253314137315500251207882642405522626505L;
      // Constant @f$ 1 / \pi @f$.
      static _GLIBCXX_CONSTEXPR _RealType __one_div_pi
      	= 0.318309886183790671537767526745028724069L;
      // Constant @f$ 2 / \pi @f$.
      static _GLIBCXX_CONSTEXPR _RealType __two_div_pi
      	= 0.636619772367581343075535053490057448138L;
      // Constant @f$ 2 / \sqrt(\pi) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __two_div_root_pi
      	= 1.128379167095512573896158903121545171688L;
      // Constant: @f$ \pi^2/6 @f$.
      static _GLIBCXX_CONSTEXPR _RealType __pi_sqr_div_6
      	= 1.644934066848226436472415166646025189221L;
      // Constant: radians to degree conversion.
      static _GLIBCXX_CONSTEXPR _RealType __deg
      	= _RealType{180} / __pi;
      // Constant: degree to radians conversion.
      static _GLIBCXX_CONSTEXPR _RealType __rad
      	= __pi / _RealType{180};

      // Constant Euler's number @f$ e @f$.
      static _GLIBCXX_CONSTEXPR _RealType __e
      	= 2.718281828459045235360287471352662497759L;
      // Constant @f$ 1 / e @f$.
      static _GLIBCXX_CONSTEXPR _RealType __one_div_e
      	= 0.367879441171442321595523770161460867446L;
      // Constant @f$ \log_2(e) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __log2_e
      	= 1.442695040888963407359924681001892137427L;
      // Constant @f$ \log_10(e) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __log10_e
      	= 0.434294481903251827651128918916605082294L;
      // Constant @f$ \ln(2) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __ln_2
      	= 0.693147180559945309417232121458176568075L;
      // Constant @f$ \ln(3) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __ln_3
      	= 1.098612288668109691395245236922525704648L;
      // Constant @f$ \ln(10) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __ln_10
      	= 2.302585092994045684017991454684364207602L;
      ///  Constant @f$ \log(\pi) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __ln_pi
      	= 1.144729885849400174143427351353058711646L;

      // Constant Euler-Mascheroni @f$ \gamma_E @f$.
      static _GLIBCXX_CONSTEXPR _RealType __gamma_e
      	= 0.577215664901532860606512090082402431043L;
      // Constant Golden Ratio @f$ \phi @f$.
      static _GLIBCXX_CONSTEXPR _RealType __phi
      	= 1.618033988749894848204586834365638117720L;
      // Catalan's constant.
      static _GLIBCXX_CONSTEXPR _RealType __catalan
      	= 0.915965594177219015054603514932384110773L;
      // Constant @f$ \sqrt(2) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_2
      	= 1.414213562373095048801688724209698078569L;
      // Constant @f$ \sqrt(3) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_3
      	= 1.732050807568877293527446341505872366945L;
      // Constant @f$ \sqrt(3) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_3_div_2
      	= 0.866025403784438646763723170752936183473L;
      // Constant @f$ \sqrt(5) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_5
      	= 2.236067977499789696409173668731276235440L;
      // Constant @f$ \sqrt(7) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __root_7
      	= 2.645751311064590590501615753639260425706L;
      // Constant @f$ 1 / \sqrt(2) @f$.
      static _GLIBCXX_CONSTEXPR _RealType __one_div_root_2
      	= 0.707106781186547524400844362104849039285L;
    };

  // And the template definitions for the constants.
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__4_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__4_pi_div_3;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__2_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__pi_half;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__pi_third;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__pi_quarter;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__root_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__cbrt_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__root_pi_div_2;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__one_div_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__two_div_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__two_div_root_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__pi_sqr_div_6;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__deg;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__rad;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__e;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__one_div_e;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__log2_e;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__log10_e;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__ln_2;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__ln_3;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__ln_10;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__ln_pi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__gamma_e;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__phi;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__catalan;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__root_2;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__root_3;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__root_5;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__root_7;
  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType __math_constants<_RealType>::__one_div_root_2;

  // The following functions mirror the constants abobe but also
  // admit generic programming with non-findamental types.
  // For fundamental types, these constexpr functions return
  // the appropriate constant above.
  // Developers of multi-precision types are encouraged to overload
  // these functions calling multiprecision "constant" functions as
  // available.

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_4_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__4_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_4_pi_div_3(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__4_pi_div_3; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_2_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__const_2_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_pi_half(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__pi_half; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_pi_third(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__pi_third; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_pi_quarter(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__pi_quarter; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_root_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__root_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_cbrt_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__cbrt_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_root_pi_div_2(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__root_pi_div_2; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_one_div_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__one_div_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_two_div_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__two_div_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_two_div_root_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__two_div_root_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_pi_sqr_div_6(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__pi_sqr_div_6; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_deg(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__deg; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_rad(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__rad; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_e(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__e; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_one_div_e(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__one_div_e; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_log2_e(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__log2_e; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_log10_e(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__log10_e; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_ln_2(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__ln_2; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_ln_3(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__ln_3; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_ln_10(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__ln_10; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_ln_pi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__ln_pi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_gamma_e(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__gamma_e;

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_phi(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__phi; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_catalan(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__catalan; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_root_2(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__root_2; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_root_3(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__root_3; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_root_5(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__root_5; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_root_7(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__root_7; }

  template<typename _RealType>
    _GLIBCXX_CONSTEXPR _RealType
    __const_one_div_root_2(_RealType = _RealType{})
    _GLIBCXX_NOEXCEPT_IF(std::is_fundamental<_RealType>::value)
    { return __math_constants<_RealType>::__one_div_root_2; }


_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // C++11

#endif // _EXT_MATH_CONST_H
