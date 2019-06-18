
  /**
   * Return the L1 norm modulus or the Manhattan metric distance of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l1_norm(const std::complex<_Tp>& __z)
    { return std::abs(std::real(__z)) + std::abs(std::imag(__z)); }
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

#ifndef _GLIBCXX_EXT_COMPLEX_NORMS_H
#define _GLIBCXX_EXT_COMPLEX_NORMS_H 1

#pragma GCC system_header

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the L2 norm modulus or the Euclidean metric distance of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l2_norm(const std::complex<_Tp>& __z)
    { return std::norm(__z); }

  /**
   * Return the Linf norm modulus of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __linf_norm(const std::complex<_Tp>& __z)
    { return std::max(std::abs(std::real(__z)), std::abs(std::imag(__z))); }


  /**
   * Return the L1 norm modulus or the Manhattan metric distance of a real number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l1_norm(_Tp __x)
    { return std::abs(__x); }

  /**
   * Return the L2 norm modulus or the Euclidean metric distance of a real number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l2_norm(_Tp __x)
    { return std::abs(__x); }

  /**
   * Return the Linf norm modulus of a real number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __linf_norm(_Tp __x)
    { return std::abs(__x); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_EXT_COMPLEX_NORMS_H
