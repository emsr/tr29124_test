// TR29124 math special functions -*- C++ -*-

// Copyright (C) 2015 Free Software Foundation, Inc.
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

/** @file bits/complex_util.h
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_COMPLEX_UTIL_H
#define _GLIBCXX_BITS_COMPLEX_UTIL_H 1

#pragma GCC system_header

#include <ratio>

  /**
   *  Return true if one component of a complex number is NaN.
   */
  template<typename _Tp>
    inline bool
    __isnan(const std::complex<_Tp> & z)
    { return __isnan(std::real(z)) || __is_nan(std::imag(z)); }


  /**
   *  Return a fraction as a real number.
   */
  template<intmax_t _Num, intmax_t _Den = 1, typename _Tp = double>
    inline constexpr _Tp
    __frac(const std::complex<_Tp> & z)
    {
      using __rat_t = std::ratio<_Num, _Den>;
      return _Tp(__rat_t::num) / _Tp(__rat_t::den);
    }


  /**
   *  Return the L1 norm modulus or the Manhattan metric distance of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __norm_L1(const std::complex<_Tp> & z)
    { return std::abs(std::real(z)) + std::abs(std::imag(z)); }

  /**
   *  Return the Linf norm modulus of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __norm_Linf(const std::complex<_Tp> & z)
    { return std::max(std::abs(std::real(z)), std::abs(std::imag(z))); }


  /**
   *  Carefully compute @c z1/z2 avoiding overflow and destructive underflow.
   *  If the quotient is successfully computed, then the logical value @c true
   *  is returned and the quotient is returned in @c z1dz2.
   *  Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Tp>
    bool
    __safe_div(std::complex<_Tp> __z1, std::complex<_Tp> __z2,
	     std::complex<_Tp> & __z1dz2);

  /**
   *  Carefully compute @c s/z2 avoiding overflow and destructive underflow.
   *  If the quotient is successfully computed, then the logical value @c true
   *  is returned and the quotient is returned in @c z1dz2.
   *  Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename _Tp>
    inline bool
    __safe_div(_Sp __s, std::complex<_Tp> __z,
	       std::complex<_Tp> & __sdz)
    { return __safe_div(std::complex<_Tp>(__s), __z, __sdz); }

  /**
   *  Carefully compute @c z1/s avoiding overflow and destructive underflow.
   *  If the quotient is successfully computed, then the logical value @c true
   *  is returned and the quotient is returned in @c z1dz2.
   *  Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename _Tp>
    inline bool
    __safe_div(std::complex<_Tp> __z, _Sp __s,
	       std::complex<_Tp> & __zds)
    { return __safe_div(__z, std::complex<_Tp>(__s), __zds); }

#include "complex_util.tcc"

#endif // _GLIBCXX_BITS_COMPLEX_UTIL_H
