// Special functions -*- C++ -*-

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

/** @file ext/complex_safe_math.h
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_EXT_COMPLEX_SAFE_MATH_H
#define _GLIBCXX_EXT_COMPLEX_SAFE_MATH_H 1

#pragma GCC system_header

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Carefully compute @c z1/z2 avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_div(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2);

  /**
   * Carefully compute @c s/z2 avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_div(_Sp __s, const std::complex<_Tp>& __z)
    { return __safe_div(std::complex<_Tp>(__s), __z); }

  /**
   * Carefully compute @c z1/s avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_div(const std::complex<_Tp>& __z, _Sp __s)
    { return __safe_div(__z, std::complex<_Tp>(__s)); }

  /**
   * @brief Carefully compute and return @c s1*s2 avoiding overflow.
   *	   If the product can be successfully computed it is returned.
   *	   Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    _Tp
    __safe_mul(_Tp __s1, _Tp __s2);

  /**
   * Carefully compute @c z1*z2 avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_mul(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2);

  /**
   * Carefully compute @c s*z avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_mul(_Sp __s, const std::complex<_Tp>& __z)
    { return __safe_mul(std::complex<_Tp>(__s), __z); }

  /**
   * Carefully compute @c z*s avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_mul(const std::complex<_Tp>& __z, _Sp __s)
    { return __safe_mul(__z, std::complex<_Tp>(__s)); }

  /**
   * Carefully compute @c z*z avoiding overflow.
   * If the square is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_sqr(const std::complex<_Tp>& __z);

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include <ext/complex_safe_math.tcc>

#endif // _GLIBCXX_EXT_COMPLEX_SAFE_MATH_H