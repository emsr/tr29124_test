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
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SPECFUN_H
#define _GLIBCXX_BITS_SPECFUN_H 1

#define __STDCPP_MATH_SPEC_FUNCS__ 201003L

#pragma GCC visibility push(default)

#include <bits/c++config.h>

#if __STDCPP_WANT_MATH_SPEC_FUNCS__ == 0
# error include <cmath> and define __STDCPP_WANT_MATH_SPEC_FUNCS__
#endif

#pragma GCC system_header

#include <limits>
#include <bits/stl_algobase.h>
#include <bits/specfun_util.h>

#if __cplusplus >= 201103L
#  include <type_traits>
#  include <complex>
#  include <bits/complex_util.h>
#  include <bits/sf_gamma.tcc>
#  include <bits/sf_bessel.tcc>
#  include <bits/sf_beta.tcc>
#  include <bits/sf_chebyshev.tcc>
#  include <bits/sf_dawson.tcc>
#  include <bits/sf_ellint.tcc>
#  include <bits/sf_expint.tcc>
#  include <bits/sf_fresnel.tcc>
#  include <bits/sf_gegenbauer.tcc>
#  include <bits/sf_hermite.tcc>
#  include <bits/sf_hyperg.tcc>
#  include <bits/sf_hypint.tcc>
#  include <bits/sf_jacobi.tcc>
#  include <bits/sf_laguerre.tcc>
#  include <bits/sf_legendre.tcc>
#  include <bits/sf_mod_bessel.tcc>
#  include <bits/sf_theta.tcc>
#  include <bits/sf_trigint.tcc>
#  include <bits/sf_zeta.tcc>
#else
#  include <tr1/type_traits>
#  include <tr1/cmath>
#  define _GLIBCXX_MATH_NS ::std::tr1::
#  include <tr1/gamma.tcc>
#  include <tr1/bessel_function.tcc>
#  include <tr1/beta_function.tcc>
#  include <tr1/ell_integral.tcc>
#  include <tr1/exp_integral.tcc>
#  include <tr1/hypergeometric.tcc>
#  include <tr1/legendre_function.tcc>
#  include <tr1/modified_bessel_func.tcc>
#  include <tr1/poly_hermite.tcc>
#  include <tr1/poly_laguerre.tcc>
#  include <tr1/riemann_zeta.tcc>
#endif

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @defgroup tr29124_math_spec_func Mathematical Special Functions
   * @ingroup numerics
   *
   * A collection of advanced mathematical special functions.
   * @{
   */

  //  Associated Laguerre polynomials

  inline float
  assoc_laguerref(unsigned int __n, unsigned int __m, float __x)
  { return __detail::__assoc_laguerre<float>(__n, __m, __x); }

  inline long double
  assoc_laguerrel(unsigned int __n, unsigned int __m, long double __x)
  { return __detail::__assoc_laguerre<long double>(__n, __m, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    assoc_laguerre(unsigned int __n, unsigned int __m, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__assoc_laguerre<__type>(__n, __m, __x);
    }

  //  Associated Legendre functions

  inline float
  assoc_legendref(unsigned int __l, unsigned int __m, float __x)
  { return __detail::__assoc_legendre_p<float>(__l, __m, __x); }

  inline long double
  assoc_legendrel(unsigned int __l, unsigned int __m, long double __x)
  { return __detail::__assoc_legendre_p<long double>(__l, __m, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    assoc_legendre(unsigned int __l, unsigned int __m, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__assoc_legendre_p<__type>(__l, __m, __x);
    }

  //  Beta functions

  inline float
  betaf(float __x, float __y)
  { return __detail::__beta<float>(__x, __y); }

  inline long double
  betal(long double __x, long double __y)
  { return __detail::__beta<long double>(__x, __y); }

  template<typename _Tpx, typename _Tpy>
    inline __gnu_cxx::__promote_num_t<_Tpx, _Tpy>
    beta(_Tpx __x, _Tpy __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpx, _Tpy>;
      return __detail::__beta<__type>(__x, __y);
    }

  //  Complete elliptic integrals of the first kind

  inline float
  comp_ellint_1f(float __k)
  { return __detail::__comp_ellint_1<float>(__k); }

  inline long double
  comp_ellint_1l(long double __k)
  { return __detail::__comp_ellint_1<long double>(__k); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    comp_ellint_1(_Tp __k)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__comp_ellint_1<__type>(__k);
    }

  //  Complete elliptic integrals of the second kind

  inline float
  comp_ellint_2f(float __k)
  { return __detail::__comp_ellint_2<float>(__k); }

  inline long double
  comp_ellint_2l(long double __k)
  { return __detail::__comp_ellint_2<long double>(__k); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    comp_ellint_2(_Tp __k)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__comp_ellint_2<__type>(__k);
    }

  //  Complete elliptic integrals of the third kind

  inline float
  comp_ellint_3f(float __k, float __nu)
  { return __detail::__comp_ellint_3<float>(__k, __nu); }

  inline long double
  comp_ellint_3l(long double __k, long double __nu)
  { return __detail::__comp_ellint_3<long double>(__k, __nu); }

  template<typename _Tp, typename _Tpn>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tpn>
    comp_ellint_3(_Tp __k, _Tpn __nu)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tpn>;
      return __detail::__comp_ellint_3<__type>(__k, __nu);
    }

  //  Regular modified cylindrical Bessel functions

  inline float
  cyl_bessel_if(float __nu, float __x)
  { return __detail::__cyl_bessel_i<float>(__nu, __x); }

  inline long double
  cyl_bessel_il(long double __nu, long double __x)
  { return __detail::__cyl_bessel_i<long double>(__nu, __x); }

  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    cyl_bessel_i(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return __detail::__cyl_bessel_i<__type>(__nu, __x);
    }

  //  Cylindrical Bessel functions (of the first kind)

  inline float
  cyl_bessel_jf(float __nu, float __x)
  { return __detail::__cyl_bessel_j<float>(__nu, __x); }

  inline long double
  cyl_bessel_jl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_j<long double>(__nu, __x); }

  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    cyl_bessel_j(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return __detail::__cyl_bessel_j<__type>(__nu, __x);
    }

  //  Irregular modified cylindrical Bessel functions

  inline float
  cyl_bessel_kf(float __nu, float __x)
  { return __detail::__cyl_bessel_k<float>(__nu, __x); }

  inline long double
  cyl_bessel_kl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_k<long double>(__nu, __x); }

  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    cyl_bessel_k(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return __detail::__cyl_bessel_k<__type>(__nu, __x);
    }

  //  Cylindrical Neumann functions

  inline float
  cyl_neumannf(float __nu, float __x)
  { return __detail::__cyl_neumann_n<float>(__nu, __x); }

  inline long double
  cyl_neumannl(long double __nu, long double __x)
  { return __detail::__cyl_neumann_n<long double>(__nu, __x); }

  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    cyl_neumann(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return __detail::__cyl_neumann_n<__type>(__nu, __x);
    }

  //  Incomplete elliptic integrals of the first kind

  inline float
  ellint_1f(float __k, float __phi)
  { return __detail::__ellint_1<float>(__k, __phi); }

  inline long double
  ellint_1l(long double __k, long double __phi)
  { return __detail::__ellint_1<long double>(__k, __phi); }

  template<typename _Tp, typename _Tpp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tpp>
    ellint_1(_Tp __k, _Tpp __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tpp>;
      return __detail::__ellint_1<__type>(__k, __phi);
    }

  //  Incomplete elliptic integrals of the second kind

  inline float
  ellint_2f(float __k, float __phi)
  { return __detail::__ellint_2<float>(__k, __phi); }

  inline long double
  ellint_2l(long double __k, long double __phi)
  { return __detail::__ellint_2<long double>(__k, __phi); }

  template<typename _Tp, typename _Tpp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tpp>
    ellint_2(_Tp __k, _Tpp __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tpp>;
      return __detail::__ellint_2<__type>(__k, __phi);
    }

  //  Incomplete elliptic integrals of the third kind

  inline float
  ellint_3f(float __k, float __nu, float __phi)
  { return __detail::__ellint_3<float>(__k, __nu, __phi); }

  inline long double
  ellint_3l(long double __k, long double __nu, long double __phi)
  { return __detail::__ellint_3<long double>(__k, __nu, __phi); }

  template<typename _Tp, typename _Tpn, typename _Tpp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tpn, _Tpp>
    ellint_3(_Tp __k, _Tpn __nu, _Tpp __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tpn, _Tpp>;
      return __detail::__ellint_3<__type>(__k, __nu, __phi);
    }

  //  Exponential integrals

  inline float
  expintf(float __x)
  { return __detail::__expint<float>(__x); }

  inline long double
  expintl(long double __x)
  { return __detail::__expint<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    expint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__expint<__type>(__x);
    }

  //  Hermite polynomials

  inline float
  hermitef(unsigned int __n, float __x)
  { return __detail::__poly_hermite<float>(__n, __x); }

  inline long double
  hermitel(unsigned int __n, long double __x)
  { return __detail::__poly_hermite<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    hermite(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__poly_hermite<__type>(__n, __x);
    }

  //  Laguerre polynomials

  inline float
  laguerref(unsigned int __n, float __x)
  { return __detail::__laguerre<float>(__n, __x); }

  inline long double
  laguerrel(unsigned int __n, long double __x)
  { return __detail::__laguerre<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    laguerre(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__laguerre<__type>(__n, __x);
    }

  //  Legendre polynomials

  inline float
  legendref(unsigned int __n, float __x)
  { return __detail::__poly_legendre_p<float>(__n, __x); }

  inline long double
  legendrel(unsigned int __n, long double __x)
  { return __detail::__poly_legendre_p<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    legendre(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__poly_legendre_p<__type>(__n, __x);
    }

  //  Riemann zeta functions

  inline float
  riemann_zetaf(float __s)
  { return __detail::__riemann_zeta<float>(__s); }

  inline long double
  riemann_zetal(long double __s)
  { return __detail::__riemann_zeta<long double>(__s); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    riemann_zeta(_Tp __s)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__riemann_zeta<__type>(__s);
    }

  //  Spherical Bessel functions

  inline float
  sph_besself(unsigned int __n, float __x)
  { return __detail::__sph_bessel<float>(__n, __x); }

  inline long double
  sph_bessell(unsigned int __n, long double __x)
  { return __detail::__sph_bessel<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sph_bessel(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__sph_bessel<__type>(__n, __x);
    }

  //  Spherical associated Legendre functions

  inline float
  sph_legendref(unsigned int __l, unsigned int __m, float __theta)
  { return __detail::__sph_legendre<float>(__l, __m, __theta); }

  inline long double
  sph_legendrel(unsigned int __l, unsigned int __m, long double __theta)
  { return __detail::__sph_legendre<long double>(__l, __m, __theta); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sph_legendre(unsigned int __l, unsigned int __m, _Tp __theta)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__sph_legendre<__type>(__l, __m, __theta);
    }

  //  Spherical Neumann functions

  inline float
  sph_neumannf(unsigned int __n, float __x)
  { return __detail::__sph_neumann<float>(__n, __x); }

  inline long double
  sph_neumannl(unsigned int __n, long double __x)
  { return __detail::__sph_neumann<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sph_neumann(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return __detail::__sph_neumann<__type>(__n, __x);
    }

  /* @} */ // tr29124_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{

  //  Confluent hypergeometric functions

  inline float
  conf_hypergf(float __a, float __c, float __x)
  { return std::__detail::__conf_hyperg<float>(__a, __c, __x); }

  inline long double
  conf_hypergl(long double __a, long double __c, long double __x)
  { return std::__detail::__conf_hyperg<long double>(__a, __c, __x); }

  template<typename _Tpa, typename _Tpc, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpa, _Tpc, _Tp>
    conf_hyperg(_Tpa __a, _Tpc __c, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpa, _Tpc, _Tp>;
      return std::__detail::__conf_hyperg<__type>(__a, __c, __x);
    }

  //  Hypergeometric functions

  inline float
  hypergf(float __a, float __b, float __c, float __x)
  { return std::__detail::__hyperg<float>(__a, __b, __c, __x); }

  inline long double
  hypergl(long double __a, long double __b, long double __c, long double __x)
  { return std::__detail::__hyperg<long double>(__a, __b, __c, __x); }

  template<typename _Tpa, typename _Tpb, typename _Tpc, typename _Tp>
    inline typename __gnu_cxx::__promote_num_t<_Tpa, _Tpb, _Tpc, _Tp>
    hyperg(_Tpa __a, _Tpb __b, _Tpc __c, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpa, _Tpb, _Tpc, _Tp>;
      return std::__detail::__hyperg<__type>(__a, __b, __c, __x);
    }

#if __cplusplus >= 201103L

  // Sinus cardinal functions

  inline float
  sincf(float __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<float>::quiet_NaN();
    else if (std::abs(__x) == std::numeric_limits<float>::infinity())
      return 0.0F;
    else
      return __x == 0.0F
	   ? 1.0F
	   : std::sin/*f*/(__x) / __x;
  }

  inline long double
  sincl(long double __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<long double>::quiet_NaN();
    else if (std::abs(__x) == std::numeric_limits<long double>::infinity())
      return 0.0L;
    else
      return __x == 0.0L
	   ? 1.0L
	   : std::sin/*l*/(__x) / __x;
  }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinc(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      if (__isnan(__x))
        return std::numeric_limits<__type>::quiet_NaN();
      else if (std::abs(__x) == std::numeric_limits<__type>::infinity())
	return __type(0);
      else
        return __type(__x) == __type(0)
             ? __type(1)
             : std::sin(__type(__x)) / __type(__x);
    }

  //  Logarithmic integrals

  inline float
  logintf(float __x)
  { return std::__detail::__logint<float>(__x); }

  inline long double
  logintl(long double __x)
  { return std::__detail::__logint<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    logint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__logint<__type>(__x);
    }

  //  Sine integrals

  inline float
  sinintf(float __x)
  { return std::__detail::__sincosint<float>(__x).first; }

  inline long double
  sinintl(long double __x)
  { return std::__detail::__sincosint<long double>(__x).first; }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sincosint<__type>(__x).first;
    }

  //  Cosine integrals

  inline float
  cosintf(float __x)
  { return std::__detail::__sincosint<float>(__x).second; }

  inline long double
  cosintl(long double __x)
  { return std::__detail::__sincosint<long double>(__x).second; }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    cosint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sincosint<__type>(__x).second;
    }

  //  Hyperbolic sine integrals

  inline float
  sinhintf(float __x)
  { return std::__detail::__sinhint<float>(__x); }

  inline long double
  sinhintl(long double __x)
  { return std::__detail::__sinhint<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinhint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sinhint<__type>(__x);
    }

  //  Hyperbolic cosine integrals

  inline float
  coshintf(float __x)
  { return std::__detail::__coshint<float>(__x); }

  inline long double
  coshintl(long double __x)
  { return std::__detail::__coshint<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    coshint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__coshint<__type>(__x);
    }

  //  Chebyshev polynomials of the first kind

  inline float
  chebyshev_tf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_t<float>(__n, __x); }

  inline long double
  chebyshev_tl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_t<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_t(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_t<__type>(__n, __x);
    }

  //  Chebyshev polynomials of the second kind

  inline float
  chebyshev_uf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_u<float>(__n, __x); }

  inline long double
  chebyshev_ul(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_u<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_u(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_u<__type>(__n, __x);
    }

  //  Chebyshev polynomials of the third kind

  inline float
  chebyshev_vf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_v<float>(__n, __x); }

  inline long double
  chebyshev_vl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_v<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_v(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_v<__type>(__n, __x);
    }

  //  Chebyshev polynomials of the fourth kind

  inline float
  chebyshev_wf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_w<float>(__n, __x); }

  inline long double
  chebyshev_wl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_w<long double>(__n, __x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_w(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_w<__type>(__n, __x);
    }

  //  Jacobi polynomials

  inline float
  jacobif(unsigned __n, float __alpha, float __beta, float __x)
  { return std::__detail::__poly_jacobi<float>(__n, __alpha, __beta, __x); }

  inline long double
  jacobil(unsigned __n, long double __alpha, long double __beta, long double __x)
  { return std::__detail::__poly_jacobi<long double>(__n, __alpha, __beta, __x); }

  template<typename _Talpha, typename _Tbeta, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Talpha, _Tbeta, _Tp>
    jacobi(unsigned __n, _Talpha __alpha, _Tbeta __beta, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Talpha, _Tbeta, _Tp>;
      return std::__detail::__poly_jacobi<__type>(__n, __alpha, __beta, __x);
    }

  //  Gegenbauer polynomials

  inline float
  gegenbauerf(unsigned int __n, float __alpha, float __x)
  { return std::__detail::__gegenbauer_poly<float>(__n, __alpha, __x); }

  inline long double
  gegenbauerl(unsigned int __n, long double __alpha, long double __x)
  { return std::__detail::__gegenbauer_poly<long double>(__n, __alpha, __x); }

  template<typename _Talpha, typename _Tp>
    inline typename __gnu_cxx::__promote_num_t<_Talpha, _Tp>
    gegenbauer(unsigned int __n, _Talpha __alpha, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Talpha, _Tp>;
      return std::__detail::__gegenbauer_poly<__type>(__n, __alpha, __x);
    }

  //  Zernike polynomials

  inline float
  zernikef(unsigned int __n, int __m, float __rho, float __phi)
  { return std::__detail::__poly_radial_jacobi(__n, std::abs(__m), __rho)
         * (__m >= 0 ? std::cos(__m * __phi) : std::sin(__m * __phi)); }

  inline long double
  zernikel(unsigned int __n, int __m, long double __rho, long double __phi)
  { return std::__detail::__poly_radial_jacobi(__n, std::abs(__m), __rho)
         * (__m >= 0 ? std::cos(__m * __phi) : std::sin(__m * __phi)); }

  template<typename _Trho, typename _Tphi>
    inline __gnu_cxx::__promote_num_t<_Trho, _Tphi>
    zernike(unsigned int __n, int __m, _Trho __rho, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Trho, _Tphi>;
      return std::__detail::__poly_radial_jacobi<__type>(__n, std::abs(__m), __rho)
           * (__m >= 0 ? std::cos(__m * __phi) : std::sin(__m * __phi));
    }

  //  Radial polynomials

  inline float
  radpolyf(unsigned int __n, unsigned int __m, float __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  inline long double
  radpolyl(unsigned int __n, unsigned int __m, long double __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    radpoly(unsigned int __n, unsigned int __m, double __rho)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__poly_radial_jacobi<__type>(__n, __m, __rho);
    }

  //  Hyperbolic sinus cardinal functions

  inline float
  sinhcf(float __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<float>::quiet_NaN();
    else
      return __x == 0.0F ? 1.0F : std::sinh/*f*/(__x) / __x;
  }

  inline long double
  sinchl(long double __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<long double>::quiet_NaN();
    else
      return __x == 0.0L ? 1.0L : std::sinh/*l*/(__x) / __x;
  }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinch(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      if (__isnan(__x))
        return std::numeric_limits<__type>::quiet_NaN();
      else
        return __type(__x) == __type(0)
             ? __type(1)
             : std::sinh(__type(__x)) / __type(__x);
    }

  //  Dawson's integrals

  inline float
  dawsonintf(float __x)
  { return std::__detail::__dawson<float>(__x); }

  inline long double
  dawsonintl(long double __x)
  { return std::__detail::__dawson<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    dawsonint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__dawson<__type>(__x);
    }

  //  Cylindrical Hankel functions of the first kind

  inline std::complex<float>
  cyl_hankel_h1f(float __nu, float __z)
  { return std::__detail::__cyl_hankel_h1<float>(__nu, __z); }

  inline std::complex<long double>
  cyl_hankel_h1l(long double __nu, long double __z)
  { return std::__detail::__cyl_hankel_h1<long double>(__nu, __z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tpnu, _Tp>>
    cyl_hankel_h1(_Tpnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_h1<__type>(__nu, __z);
    }

  //  Cylindrical Hankel functions of the second kind

  inline std::complex<float>
  cyl_hankel_h2f(float __nu, float __z)
  { return std::__detail::__cyl_hankel_h2<float>(__nu, __z); }

  inline std::complex<long double>
  cyl_hankel_h2l(long double __nu, long double __z)
  { return std::__detail::__cyl_hankel_h2<long double>(__nu, __z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tpnu, _Tp>>
    cyl_hankel_h2(_Tpnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_h2<__type>(__nu, __z);
    }

  //  Spherical Hankel functions of the first kind

  inline std::complex<float>
  sph_hankel_h1f(unsigned int __n, float __z)
  { return std::__detail::__sph_hankel_h1<float>(__n, __z); }

  inline std::complex<long double>
  sph_hankel_h1l(unsigned int __n, long double __z)
  { return std::__detail::__sph_hankel_h1<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    sph_hankel_h1(unsigned int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sph_hankel_h1<__type>(__n, __z);
    }

  //  Spherical Hankel functions of the second kind

  inline std::complex<float>
  sph_hankel_h2f(unsigned int __n, float __z)
  { return std::__detail::__sph_hankel_h2<float>(__n, __z); }

  inline std::complex<long double>
  sph_hankel_h2l(unsigned int __n, long double __z)
  { return std::__detail::__sph_hankel_h2<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    sph_hankel_h2(unsigned int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sph_hankel_h2<__type>(__n, __z);
    }

  //  Modified spherical Bessel functions of the first kind

  inline float
  sph_bessel_if(unsigned int __n, float __x)
  {
    float __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<float>(__n, __x,
        				  __i_n, __k_n, __ip_n, __kp_n);
    return __i_n;
  }

  inline long double
  sph_bessel_il(unsigned int __n, long double __x)
  {
    long double __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<long double>(__n, __x,
        					__i_n, __k_n, __ip_n, __kp_n);
    return __i_n;
  }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sph_bessel_i(unsigned int __n, double __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __i_n, __k_n, __ip_n, __kp_n;
      std::__detail::__sph_bessel_ik<__type>(__n, __x,
        				     __i_n, __k_n, __ip_n, __kp_n);
      return __i_n;
    }

  //  Modified spherical Bessel functions of the second kind

  inline float
  sph_bessel_kf(unsigned int __n, float __x)
  {
    float __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<float>(__n, __x,
        				  __i_n, __k_n, __ip_n, __kp_n);
    return __k_n;
  }

  inline long double
  sph_bessel_kl(unsigned int __n, long double __x)
  {
    long double __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<long double>(__n, __x,
        					__i_n, __k_n, __ip_n, __kp_n);
    return __k_n;
  }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sph_bessel_k(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __i_n, __k_n, __ip_n, __kp_n;
      std::__detail::__sph_bessel_ik<__type>(__n, __x,
        				     __i_n, __k_n, __ip_n, __kp_n);
      return __k_n;
    }

  //  Airy functions of the first kind

  inline float
  airy_aif(float __x)
  {
    float __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<float>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Ai;
  }

  inline long double
  airy_ail(long double __x)
  {
    long double __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<long double>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Ai;
  }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    airy_ai(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __Ai, __Bi, __Aip, __Bip;
      std::__detail::__airy<__type>(__x, __Ai, __Bi, __Aip, __Bip);
      return __Ai;
    }

  //  Airy functions of the second kind

  inline float
  airy_bif(float __x)
  {
    float __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<float>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Bi;
  }

  inline long double
  airy_bil(long double __x)
  {
    long double __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<long double>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Bi;
  }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    airy_bi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __Ai, __Bi, __Aip, __Bip;
      std::__detail::__airy<__type>(__x, __Ai, __Bi, __Aip, __Bip);
      return __Bi;
    }

  //  Upper incomplete gamma functions

  inline float
  gamma_uf(float __n, float __x)
  { return std::__detail::__gamma_u<float>(__n, __x); }

  inline long double
  gamma_ul(long double __n, long double __x)
  { return std::__detail::__gamma_u<long double>(__n, __x); }

  template<typename _Tn, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tn, _Tp>
    gamma_u(_Tn __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tn, _Tp>;
      return std::__detail::__gamma_u<__type>(__n, __x);
    }

  //  Lower incomplete gamma functions

  inline float
  gamma_lf(float __n, float __x)
  { return std::__detail::__gamma_l<float>(__n, __x); }

  inline long double
  gamma_ll(long double __n, long double __x)
  { return std::__detail::__gamma_l<long double>(__n, __x); }

  template<typename _Tn, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tn, _Tp>
    gamma_l(_Tn __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tn, _Tp>;
      return std::__detail::__gamma_l<__type>(__n, __x);
    }

  //  Digamma functions

  inline float
  digammaf(float __z)
  { return std::__detail::__psi<float>(__z); }

  inline long double
  digammal(long double __z)
  { return std::__detail::__psi<long double>(__z); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    digamma(_Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__psi<__type>(__z);
    }

  //  Dilogarithm functions

  inline float
  dilogf(float __x)
  { return std::__detail::__dilog<float>(__x); }

  inline long double
  dilogl(long double __x)
  { return std::__detail::__dilog<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    dilog(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__dilog<__type>(__x);
    }

  //  Complete Carlson elliptic R_F functions

  inline float
  comp_ellint_rf(float __x, float __y)
  { return std::__detail::__comp_ellint_rf<float>(__x, __y); }

  inline long double
  comp_ellint_rf(long double __x, long double __y)
  { return std::__detail::__comp_ellint_rf<long double>(__x, __y); }

  template<typename _Tx, typename _Ty>
    inline __gnu_cxx::__promote_num_t<_Tx, _Ty>
    comp_ellint_rf(_Tx __x, _Ty __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tx, _Ty>;
      return std::__detail::__comp_ellint_rf<__type>(__x, __y);
    }

  //  Carlson elliptic R_F functions

  inline float
  ellint_rff(float __x, float __y, float __z)
  { return std::__detail::__ellint_rf<float>(__x, __y, __z); }

  inline long double
  ellint_rfl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rf<long double>(__x, __y, __z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>
    ellint_rf(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rf<__type>(__x, __y, __z);
    }

  //  Carlson elliptic R_C functions

  inline float
  ellint_rcf(float __x, float __y)
  { return std::__detail::__ellint_rc<float>(__x, __y); }

  inline long double
  ellint_rcl(long double __x, long double __y)
  { return std::__detail::__ellint_rc<long double>(__x, __y); }

  template<typename _Tp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up>
    ellint_rc(_Tp __x, _Up __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up>;
      return std::__detail::__ellint_rc<__type>(__x, __y);
    }

  //  Carlson elliptic R_J functions

  inline float
  ellint_rjf(float __x, float __y, float __z, float __p)
  { return std::__detail::__ellint_rj<float>(__x, __y, __z, __p); }

  inline long double
  ellint_rjl(long double __x, long double __y, long double __z, long double __p)
  { return std::__detail::__ellint_rj<long double>(__x, __y, __z, __p); }

  template<typename _Tp, typename _Up, typename _Vp, typename _Wp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp, _Wp>
    ellint_rj(_Tp __x, _Up __y, _Vp __z, _Wp __p)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp, _Wp>;
      return std::__detail::__ellint_rj<__type>(__x, __y, __z, __p);
    }

  //  Carlson elliptic R_D functions

  inline float
  ellint_rdf(float __x, float __y, float __z)
  { return std::__detail::__ellint_rd<float>(__x, __y, __z); }

  inline long double
  ellint_rdl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rd<long double>(__x, __y, __z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>
    ellint_rd(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rd<__type>(__x, __y, __z);
    }

  //  Complete Carlson elliptic R_G functions

  inline float
  comp_ellint_rg(float __x, float __y)
  { return std::__detail::__comp_ellint_rg<float>(__x, __y); }

  inline long double
  comp_ellint_rg(long double __x, long double __y)
  { return std::__detail::__comp_ellint_rg<long double>(__x, __y); }

  template<typename _Tx, typename _Ty>
    inline __gnu_cxx::__promote_num_t<_Tx, _Ty>
    comp_ellint_rg(_Tx __x, _Ty __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tx, _Ty>;
      return std::__detail::__comp_ellint_rg<__type>(__x, __y);
    }

  //  Carlson elliptic R_G functions

  inline float
  ellint_rgf(float __x, float __y, float __z)
  { return std::__detail::__ellint_rg<float>(__x, __y, __z); }

  inline long double
  ellint_rgl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rg<long double>(__x, __y, __z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>
    ellint_rg(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rg<__type>(__x, __y, __z);
    }

  //  Hurwitz zeta functions

  inline float
  hurwitz_zetaf(float __s, float __a)
  { return std::__detail::__hurwitz_zeta<float>(__s, __a); }

  inline long double
  hurwitz_zetal(long double __s, long double __a)
  { return std::__detail::__hurwitz_zeta<long double>(__s, __a); }

  template<typename _Tp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up>
    hurwitz_zeta(_Tp __s, _Up __a)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up>;
      return std::__detail::__hurwitz_zeta<__type>(__s, __a);
    }

#endif // __cplusplus >= 201103L

} // namespace __gnu_cxx

#pragma GCC visibility pop

#endif // _GLIBCXX_BITS_SPECFUN_H
