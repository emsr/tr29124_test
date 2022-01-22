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

#ifndef COMPLEX_NORMS_H
#define COMPLEX_NORMS_H 1

#include <complex>

namespace emsr
{

  /**
   * Return the L1 norm modulus or the Manhattan metric distance of a complex number.
   */
  template<typename Tp>
    inline constexpr Tp
    l1_norm(const std::complex<Tp>& z)
    { return std::abs(std::real(z)) + std::abs(std::imag(z)); }

  /**
   * Return the L2 norm modulus or the Euclidean metric distance of a complex number.
   */
  template<typename Tp>
    inline constexpr Tp
    l2_norm(const std::complex<Tp>& z)
    { return std::norm(z); }

  /**
   * Return the Linf norm modulus of a complex number.
   */
  template<typename Tp>
    inline constexpr Tp
    linf_norm(const std::complex<Tp>& z)
    { return std::max(std::abs(std::real(z)), std::abs(std::imag(z))); }


  /**
   * Return the L1 norm modulus or the Manhattan metric distance of a real number.
   */
  template<typename Tp>
    inline constexpr Tp
    l1_norm(Tp x)
    { return std::abs(x); }

  /**
   * Return the L2 norm modulus or the Euclidean metric distance of a real number.
   */
  template<typename Tp>
    inline constexpr Tp
    l2_norm(Tp x)
    { return std::abs(x); }

  /**
   * Return the Linf norm modulus of a real number.
   */
  template<typename Tp>
    inline constexpr Tp
    linf_norm(Tp x)
    { return std::abs(x); }

} // namespace emsr

#endif // COMPLEX_NORMS_H
