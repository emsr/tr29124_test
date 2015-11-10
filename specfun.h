// math special functions -*- C++ -*-

// Copyright (C) 2006-2015 Free Software Foundation, Inc.
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

/** @file bits/specfun.h
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SPECIAL_FUNCTION_H
#define _GLIBCXX_BITS_SPECIAL_FUNCTION_H 1

#pragma GCC system_header

#include <bits/stl_algobase.h>
#include <limits>
#if __cplusplus >= 201103L
#  include <type_traits>
#  define TR1NS std::
#else
#  include <tr1/type_traits>
#  include <tr1/cmath>
#  define TR1NS std::tr1::
#endif

#include <complex>

#include <bits/sf_gamma.tcc>
#include <bits/sf_bessel.tcc>
#include <bits/sf_beta.tcc>
#include <bits/sf_ellint.tcc>
#include <bits/sf_expint.tcc>
#include <bits/sf_hyperg.tcc>
#include <bits/sf_legendre.tcc>
#include <bits/sf_mod_bessel.tcc>
#include <bits/sf_hermite.tcc>
#include <bits/sf_laguerre.tcc>
#include <bits/sf_zeta.tcc>

namespace std _GLIBCXX_VISIBILITY(default)
{

  /**
   * @defgroup tr29124_math_spec_func Mathematical Special Functions
   * @ingroup numerics
   *
   * A collection of advanced mathematical special functions.
   * @{
   */

  inline float
  assoc_laguerref(unsigned int __n, unsigned int __m, float __x)
  { return __detail::__assoc_laguerre<float>(__n, __m, __x); }

  inline long double
  assoc_laguerrel(unsigned int __n, unsigned int __m, long double __x)
  { return __detail::__assoc_laguerre<long double>(__n, __m, __x); }

  ///  Associated Laguerre polynomials.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    assoc_laguerre(unsigned int __n, unsigned int __m, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__assoc_laguerre<__type>(__n, __m, __x);
    }

  inline float
  assoc_legendref(unsigned int __l, unsigned int __m, float __x)
  { return __detail::__assoc_legendre_p<float>(__l, __m, __x); }

  inline long double
  assoc_legendrel(unsigned int __l, unsigned int __m, long double __x)
  { return __detail::__assoc_legendre_p<long double>(__l, __m, __x); }

  ///  Associated Legendre functions.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    assoc_legendre(unsigned int __l, unsigned int __m, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__assoc_legendre_p<__type>(__l, __m, __x);
    }

  inline float
  betaf(float __x, float __y)
  { return __detail::__beta<float>(__x, __y); }

  inline long double
  betal(long double __x, long double __y)
  { return __detail::__beta<long double>(__x, __y); }

  ///  Beta functions.
  template<typename _Tpx, typename _Tpy>
    inline typename __gnu_cxx::__promote_2<_Tpx, _Tpy>::__type
    beta(_Tpx __x, _Tpy __y)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpx, _Tpy>::__type __type;
      return __detail::__beta<__type>(__x, __y);
    }

  inline float
  comp_ellint_1f(float __k)
  { return __detail::__comp_ellint_1<float>(__k); }

  inline long double
  comp_ellint_1l(long double __k)
  { return __detail::__comp_ellint_1<long double>(__k); }

  ///  Complete elliptic integrals of the first kind.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    comp_ellint_1(_Tp __k)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__comp_ellint_1<__type>(__k);
    }

  inline float
  comp_ellint_2f(float __k)
  { return __detail::__comp_ellint_2<float>(__k); }

  inline long double
  comp_ellint_2l(long double __k)
  { return __detail::__comp_ellint_2<long double>(__k); }

  ///  Complete elliptic integrals of the second kind.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    comp_ellint_2(_Tp __k)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__comp_ellint_2<__type>(__k);
    }

  inline float
  comp_ellint_3f(float __k, float __nu)
  { return __detail::__comp_ellint_3<float>(__k, __nu); }

  inline long double
  comp_ellint_3l(long double __k, long double __nu)
  { return __detail::__comp_ellint_3<long double>(__k, __nu); }

  ///  Complete elliptic integrals of the third kind.
  template<typename _Tp, typename _Tpn>
    inline typename __gnu_cxx::__promote_2<_Tp, _Tpn>::__type
    comp_ellint_3(_Tp __k, _Tpn __nu)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Tpn>::__type __type;
      return __detail::__comp_ellint_3<__type>(__k, __nu);
    }

  inline float
  cyl_bessel_if(float __nu, float __x)
  { return __detail::__cyl_bessel_i<float>(__nu, __x); }

  inline long double
  cyl_bessel_il(long double __nu, long double __x)
  { return __detail::__cyl_bessel_i<long double>(__nu, __x); }

  ///  Regular modified cylindrical Bessel functions.
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_bessel_i(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_bessel_i<__type>(__nu, __x);
    }

  inline float
  cyl_bessel_jf(float __nu, float __x)
  { return __detail::__cyl_bessel_j<float>(__nu, __x); }

  inline long double
  cyl_bessel_jl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_j<long double>(__nu, __x); }

  ///  Cylindrical Bessel functions (of the first kind).
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_bessel_j(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_bessel_j<__type>(__nu, __x);
    }

  inline float
  cyl_bessel_kf(float __nu, float __x)
  { return __detail::__cyl_bessel_k<float>(__nu, __x); }

  inline long double
  cyl_bessel_kl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_k<long double>(__nu, __x); }

  ///  Irregular modified cylindrical Bessel functions.
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_bessel_k(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_bessel_k<__type>(__nu, __x);
    }

  inline float
  cyl_neumannf(float __nu, float __x)
  { return __detail::__cyl_neumann_n<float>(__nu, __x); }

  inline long double
  cyl_neumannl(long double __nu, long double __x)
  { return __detail::__cyl_neumann_n<long double>(__nu, __x); }

  ///  Cylindrical Neumann functions.
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_neumann(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_neumann_n<__type>(__nu, __x);
    }

  inline float
  ellint_1f(float __k, float __phi)
  { return __detail::__ellint_1<float>(__k, __phi); }

  inline long double
  ellint_1l(long double __k, long double __phi)
  { return __detail::__ellint_1<long double>(__k, __phi); }

  ///  Incomplete elliptic integrals of the first kind.
  template<typename _Tp, typename _Tpp>
    inline typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type
    ellint_1(_Tp __k, _Tpp __phi)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type __type;
      return __detail::__ellint_1<__type>(__k, __phi);
    }

  inline float
  ellint_2f(float __k, float __phi)
  { return __detail::__ellint_2<float>(__k, __phi); }

  inline long double
  ellint_2l(long double __k, long double __phi)
  { return __detail::__ellint_2<long double>(__k, __phi); }

  ///  Incomplete elliptic integrals of the second kind.
  template<typename _Tp, typename _Tpp>
    inline typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type
    ellint_2(_Tp __k, _Tpp __phi)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type __type;
      return __detail::__ellint_2<__type>(__k, __phi);
    }

  inline float
  ellint_3f(float __k, float __nu, float __phi)
  { return __detail::__ellint_3<float>(__k, __nu, __phi); }

  inline long double
  ellint_3l(long double __k, long double __nu, long double __phi)
  { return __detail::__ellint_3<long double>(__k, __nu, __phi); }

  ///  Incomplete elliptic integrals of the third kind.
  template<typename _Tp, typename _Tpn, typename _Tpp>
    inline typename __gnu_cxx::__promote_3<_Tp, _Tpn, _Tpp>::__type
    ellint_3(_Tp __k, _Tpn __nu, _Tpp __phi)
    {
      typedef typename __gnu_cxx::__promote_3<_Tp, _Tpn, _Tpp>::__type __type;
      return __detail::__ellint_3<__type>(__k, __nu, __phi);
    }

  inline float
  expintf(float __x)
  { return __detail::__expint<float>(__x); }

  inline long double
  expintl(long double __x)
  { return __detail::__expint<long double>(__x); }

  ///  Exponential integrals.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    expint(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__expint<__type>(__x);
    }

  inline float
  hermitef(unsigned int __n, float __x)
  { return __detail::__poly_hermite<float>(__n, __x); }

  inline long double
  hermitel(unsigned int __n, long double __x)
  { return __detail::__poly_hermite<long double>(__n, __x); }

  ///  Hermite polynomials.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    hermite(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__poly_hermite<__type>(__n, __x);
    }

  inline float
  laguerref(unsigned int __n, float __x)
  { return __detail::__laguerre<float>(__n, __x); }

  inline long double
  laguerrel(unsigned int __n, long double __x)
  { return __detail::__laguerre<long double>(__n, __x); }

  ///  Laguerre polynomials.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    laguerre(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__laguerre<__type>(__n, __x);
    }

  inline float
  legendref(unsigned int __n, float __x)
  { return __detail::__poly_legendre_p<float>(__n, __x); }

  inline long double
  legendrel(unsigned int __n, long double __x)
  { return __detail::__poly_legendre_p<long double>(__n, __x); }

  ///  Legendre polynomials.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    legendre(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__poly_legendre_p<__type>(__n, __x);
    }

  inline float
  riemann_zetaf(float __x)
  { return __detail::__riemann_zeta<float>(__x); }

  inline long double
  riemann_zetal(long double __x)
  { return __detail::__riemann_zeta<long double>(__x); }

  ///  Riemann zeta function.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    riemann_zeta(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__riemann_zeta<__type>(__x);
    }

  inline float
  sph_besself(unsigned int __n, float __x)
  { return __detail::__sph_bessel<float>(__n, __x); }

  inline long double
  sph_bessell(unsigned int __n, long double __x)
  { return __detail::__sph_bessel<long double>(__n, __x); }

  ///  Spherical Bessel functions.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_bessel(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_bessel<__type>(__n, __x);
    }

  inline float
  sph_legendref(unsigned int __l, unsigned int __m, float __theta)
  { return __detail::__sph_legendre<float>(__l, __m, __theta); }

  inline long double
  sph_legendrel(unsigned int __l, unsigned int __m, long double __theta)
  { return __detail::__sph_legendre<long double>(__l, __m, __theta); }

  ///  Spherical associated Legendre functions.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_legendre(unsigned int __l, unsigned int __m, _Tp __theta)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_legendre<__type>(__l, __m, __theta);
    }

  inline float
  sph_neumannf(unsigned int __n, float __x)
  { return __detail::__sph_neumann<float>(__n, __x); }

  inline long double
  sph_neumannl(unsigned int __n, long double __x)
  { return __detail::__sph_neumann<long double>(__n, __x); }

  ///  Spherical Neumann functions.
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_neumann(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_neumann<__type>(__n, __x);
    }

  /* @} */ // tr29124_math_spec_func

} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{

  inline float
  conf_hypergf(float __a, float __c, float __x)
  { return std::__detail::__conf_hyperg<float>(__a, __c, __x); }

  inline long double
  conf_hypergl(long double __a, long double __c, long double __x)
  { return std::__detail::__conf_hyperg<long double>(__a, __c, __x); }

  ///  Confluent hypergeometric functions.
  template<typename _Tpa, typename _Tpc, typename _Tp>
    inline typename __gnu_cxx::__promote_3<_Tpa, _Tpc, _Tp>::__type
    conf_hyperg(_Tpa __a, _Tpc __c, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_3<_Tpa, _Tpc, _Tp>::__type;
      return std::__detail::__conf_hyperg<__type>(__a, __c, __x);
    }

  inline float
  hypergf(float __a, float __b, float __c, float __x)
  { return std::__detail::__hyperg<float>(__a, __b, __c, __x); }

  inline long double
  hypergl(long double __a, long double __b, long double __c, long double __x)
  { return std::__detail::__hyperg<long double>(__a, __b, __c, __x); }

  ///  Hypergeometric functions.
  template<typename _Tpa, typename _Tpb, typename _Tpc, typename _Tp>
    inline typename __gnu_cxx::__promote_4<_Tpa, _Tpb, _Tpc, _Tp>::__type
    hyperg(_Tpa __a, _Tpb __b, _Tpc __c, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_4<_Tpa, _Tpb, _Tpc, _Tp>::__type;
      return std::__detail::__hyperg<__type>(__a, __b, __c, __x);
    }


  inline std::complex<float>
  sph_hankel_h1f(int __n, float __z)
  { return std::__detail::__sph_hankel_h1<float>(__n, __z); }

  inline std::complex<long double>
  sph_hankel_h1l(int __n, long double __z)
  { return std::__detail::__sph_hankel_h1<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote<_Tp>::__type>
    sph_hankel_h1(int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_hankel_h1<__type>(__n, __z);
    }


  inline std::complex<float>
  sph_hankel_h2f(int __n, float __z)
  { return std::__detail::__sph_hankel_h2<float>(__n, __z); }

  inline std::complex<long double>
  sph_hankel_h2l(int __n, std::complex<long double> __z)
  { return std::__detail::__sph_hankel_h2<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote<_Tp>::__type>
    sph_hankel_h2(int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_hankel_h2<__type>(__n, __z);
    }


  inline float
  sph_besself(int __n, float __z)
  { return std::__detail::__sph_bessel<float>(__n, __z); }

  inline long double
  sph_bessell(int __n, long double __z)
  { return std::__detail::__sph_bessel<long double>(__n, __z); }

  template<typename _Tp>
    inline __gnu_cxx::__promote<_Tp>::__type
    sph_bessel(int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_bessel<__type>(__n, __z);
    }


  inline float
  sph_neumannf(int __n, float __z)
  { return std::__detail::__sph_neumann<float>(__n, __z); }

  inline long double
  sph_neumannl(int __n, long double __z)
  { return std::__detail::__sph_neumann<long double>(__n, __z); }

  template<typename _Tp>
    inline __gnu_cxx::__promote<_Tp>::__type
    sph_neumann(int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_neumann<__type>(__n, __z);
    }


  inline float
  airy_aif(float __z)
  { return std::__detail::__airy_ai<float>(__z); }

  inline long double
  airy_ail(long double __z)
  { return std::__detail::__airy_ai<long double>(__z); }

  template<typename _Tp>
    inline __gnu_cxx::__promote<_Tp>::__type
    airy_ai(_Tp __z)
    {
      using __type = __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__airy_ai<__type>(__z);
    }


  inline float
  airy_bif(float __z)
  { return std::__detail::__airy_bi<float>(__z); }

  inline long double
  airy_bil(long double __z)
  { return std::__detail::__airy_bi<long double>(__z); }

  template<typename _Tp>
    inline __gnu_cxx::__promote<_Tp>::__type
    airy_bi(_Tp __z)
    {
      using __type = __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail:__airy_bi<__type>(__z);
    }

} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SPECIAL_FUNCTION_H
