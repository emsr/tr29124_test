// -*- C++ -*- header.

// Copyright (C) 2016 Free Software Foundation, Inc.
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

/** @file bits/specfun128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SPECFUN128_H
#define SPECFUN128_H 1

#pragma GCC system_header

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

namespace std
{

  inline __float128
  assoc_laguerreq(unsigned int __n, unsigned int __m, __float128 __x)
  { return __detail::__assoc_laguerre<__float128>(__n, __m, __x); }

  inline __float128
  assoc_legendreq(unsigned int __l, unsigned int __m, __float128 __x)
  { return __detail::__assoc_legendre_p<__float128>(__l, __m, __x); }

  inline __float128
  betaq(__float128 __x, __float128 __y)
  { return __detail::__beta<__float128>(__x, __y); }

  inline __float128
  comp_ellint_1q(__float128 __k)
  { return __detail::__comp_ellint_1<__float128>(__k); }

  inline __float128
  comp_ellint_2q(__float128 __k)
  { return __detail::__comp_ellint_2<__float128>(__k); }

  inline __float128
  comp_ellint_3q(__float128 __k, __float128 __nu)
  { return __detail::__comp_ellint_3<__float128>(__k, __nu); }

  inline __float128
  cyl_bessel_iq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_bessel_i<__float128>(__nu, __x); }

  inline __float128
  cyl_bessel_jq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_bessel_j<__float128>(__nu, __x); }

  inline __float128
  cyl_bessel_kq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_bessel_k<__float128>(__nu, __x); }

  inline __float128
  cyl_neumannq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_neumann_n<__float128>(__nu, __x); }

  inline __float128
  ellint_1q(__float128 __k, __float128 __phi)
  { return __detail::__ellint_1<__float128>(__k, __phi); }

  inline __float128
  ellint_2q(__float128 __k, __float128 __phi)
  { return __detail::__ellint_2<__float128>(__k, __phi); }

  inline __float128
  ellint_3q(__float128 __k, __float128 __nu, __float128 __phi)
  { return __detail::__ellint_3<__float128>(__k, __nu, __phi); }

  inline __float128
  expintq(__float128 __x)
  { return __detail::__expint<__float128>(__x); }

  inline __float128
  hermiteq(unsigned int __n, __float128 __x)
  { return __detail::__poly_hermite<__float128>(__n, __x); }

  inline __float128
  laguerreq(unsigned int __n, __float128 __x)
  { return __detail::__laguerre<__float128>(__n, __x); }

  inline __float128
  legendreq(unsigned int __n, __float128 __x)
  { return __detail::__poly_legendre_p<__float128>(__n, __x); }

  inline __float128
  riemann_zetaq(__float128 __s)
  { return __detail::__riemann_zeta<__float128>(__s); }

  inline __float128
  sph_besselq(unsigned int __n, __float128 __x)
  { return __detail::__sph_bessel<__float128>(__n, __x); }

  inline __float128
  sph_legendreq(unsigned int __l, unsigned int __m, __float128 __theta)
  { return __detail::__sph_legendre<__float128>(__l, __m, __theta); }

  inline __float128
  sph_neumannq(unsigned int __n, __float128 __x)
  { return __detail::__sph_neumann<__float128>(__n, __x); }

} // namespace std

namespace __gnu_cxx
{
} // namespace __gnu_cxx

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // SPECFUN128_H
