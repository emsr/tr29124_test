// Special functions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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
 */

#ifndef COMPLEX_SAFE_MATH_H
#define COMPLEX_SAFE_MATH_H 1

#include <complex>

namespace emsr
{

  /**
   * Carefully compute @c z1/z2 avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename Tp>
    std::complex<Tp>
    safe_div(const std::complex<Tp>& z1, const std::complex<Tp>& z2);

  /**
   * Carefully compute @c s/z2 avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename Tp>
    inline std::complex<Tp>
    safe_div(_Sp s, const std::complex<Tp>& z)
    { return safe_div(std::complex<Tp>(s), z); }

  /**
   * Carefully compute @c z1/s avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename Tp>
    inline std::complex<Tp>
    safe_div(const std::complex<Tp>& z, _Sp s)
    { return safe_div(z, std::complex<Tp>(s)); }

  /**
   * @brief Carefully compute and return @c s1*s2 avoiding overflow.
   *	   If the product can be successfully computed it is returned.
   *	   Otherwise, std::runtime_error is thrown.
   */
  template<typename Tp>
    Tp
    safe_mul(Tp s1, Tp s2);

  /**
   * Carefully compute @c z1*z2 avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename Tp>
    std::complex<Tp>
    safe_mul(const std::complex<Tp>& z1, const std::complex<Tp>& z2);

  /**
   * Carefully compute @c s*z avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Sp, typename Tp>
    inline std::complex<Tp>
    safe_mul(_Sp s, const std::complex<Tp>& z)
    { return safe_mul(std::complex<Tp>(s), z); }

  /**
   * Carefully compute @c z*s avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Sp, typename Tp>
    inline std::complex<Tp>
    safe_mul(const std::complex<Tp>& z, _Sp s)
    { return safe_mul(z, std::complex<Tp>(s)); }

  /**
   * Carefully compute @c z*z avoiding overflow.
   * If the square is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename Tp>
    std::complex<Tp>
    safe_sqr(const std::complex<Tp>& z);

} // namespace emsr

#include <emsr/complex_safe_math.tcc>

#endif // COMPLEX_SAFE_MATH_H
