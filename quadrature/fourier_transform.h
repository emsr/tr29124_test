// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2017 Free Software Foundation, Inc.
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
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.


#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H 1

namespace __gnu_cxx
{

/**
 * Fast Fourier Transform
 *
 * Discrete Fourier Transform can be regarded as evaluating a
 * polynomial of degree N-1 on the powers
 * omega^0, omega, omega^2, ..., omega^(N-1) where omega is a the
 * Nth root of unity.
 *
 * Given a polynomial of even degree
 *
 * p(t) = a_0 + a_1 t + a_2 t^2 + ...
 *
 * we find:
 *
 * p(t) = a_0 + a_2 t^2 + a_4 t^4 + ... + t(a_1 + a_3 t^2 + ...)
 *      = q_e(t^2) + t q_o (t^2)
 *
 * where q_e and q_o are polynomials formed with the even and odd
 * coefficients of p(t). Thus, we get
 *
 * p(1) = q_e(1) + omega^0 q_o(1)
 * p(omega) = q_e(omega^2) + omega^1 q_o(omega^2)
 * p(omega^2) = q_e(omega^4) + omega^2 q_o(omega^4)
 * ...
 *
 * Note how on the RHS the Fourier transforms of q_e and q_o appear. Thus:
 */

  /**
   * Discreet Fourier transform on complex data.
   */
  template<typename _Tp>
    void
    __discrete_fourier_transform(bool __do_forward,
				 std::vector<std::complex<_Tp>>& __z);

  /**
   * Fast Fourier Transform on complex data.
   */
  template<typename _Tp>
    void
    fast_fourier_transform(std::vector<std::complex<_Tp>>& __z);

  /**
   * Inverse Fast Fourier Transform on complex data.
   */
  template<typename _Tp>
    void
    inv_fast_fourier_transform(std::vector<std::complex<_Tp>>& __z);

  /**
   * Fast Fourier Transform on real data.
   */
  template<typename _Tp>
    void
    fast_fourier_transform(std::vector<_Tp>& __x);

  /**
   * Inverse Fast Fourier Transform on real data.
   */
  template<typename _Tp>
    void
    inv_fast_fourier_transform(std::vector<_Tp>& __x);

  /**
   * Fast Fourier Transform on input range.
   */
  template <typename _CmplxIter>
    void
    fast_fourier_transform(const _CmplxIter& __from, const _CmplxIter& __to);

  /**
   * Inverse Fast Fourier Transform on input range.
   */
  template <typename _CmplxIter>
    void
    inv_fast_fourier_transform(const _CmplxIter& __from, const _CmplxIter& __to);

} // namespace __gnu_cxx

#include "fourier_transform.tcc"

#endif // FOURIER_TRANSFORM_H
