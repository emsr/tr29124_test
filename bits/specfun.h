// math special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
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
 * You should not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SPECFUN_H
#define _GLIBCXX_BITS_SPECFUN_H 1

#pragma GCC visibility push(default)

#include <bits/c++config.h>

#define __STDCPP_MATH_SPEC_FUNCS__ 201003L

#define __cpp_lib_math_special_functions 201603L

#if __cplusplus <= 201402L && __STDCPP_WANT_MATH_SPEC_FUNCS__ == 0
# error include <cmath> and define __STDCPP_WANT_MATH_SPEC_FUNCS__
#endif

#pragma GCC system_header

#include <limits>
#include <bits/stl_algobase.h>
#include <bits/specfun_util.h>

#if __cplusplus >= 201103L
#  include <type_traits>
#  include <complex>
#  include <bits/numeric_limits.h>
#  include <bits/complex_util.h>
#  include <bits/sf_gamma.tcc>
#  include <bits/sf_bessel.tcc>
#  include <bits/sf_beta.tcc>
#  include <bits/sf_cardinal.tcc>
#  include <bits/sf_chebyshev.tcc>
#  include <bits/sf_dawson.tcc>
#  include <bits/sf_ellint.tcc>
#  include <bits/sf_expint.tcc>
#  include <bits/sf_fresnel.tcc>
#  include <bits/sf_gegenbauer.tcc>
#  include <bits/sf_hyperg.tcc>
#  include <bits/sf_hypint.tcc>
#  include <bits/sf_jacobi.tcc>
#  include <bits/sf_laguerre.tcc>
#  include <bits/sf_legendre.tcc>
#  include <bits/sf_hydrogen.tcc> // Needs __sph_legendre.
#  include <bits/sf_mod_bessel.tcc>
#  include <bits/sf_hermite.tcc> // Needs __airy.
#  include <bits/sf_theta.tcc>
#  include <bits/sf_trigint.tcc>
#  include <bits/sf_zeta.tcc>
#  include <bits/sf_owens_t.tcc>
#  include <bits/sf_polylog.tcc>
#  include <bits/sf_airy.tcc>
#  include <bits/sf_hankel.tcc>
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

  // Associated Laguerre polynomials

  /**
   * Return the associated Laguerre polynomial of order @c n,
   * degree @c m: @f$ L_n^m(x) @f$ for @c float argument.
   *
   * @see assoc_laguerre for more details.
   */
  inline float
  assoc_laguerref(unsigned int __n, unsigned int __m, float __x)
  { return __detail::__assoc_laguerre<float>(__n, __m, __x); }

  /**
   * Return the associated Laguerre polynomial of order @c n,
   * degree @c m: @f$ L_n^m(x) @f$.
   *
   * @see assoc_laguerre for more details.
   */
  inline long double
  assoc_laguerrel(unsigned int __n, unsigned int __m, long double __x)
  { return __detail::__assoc_laguerre<long double>(__n, __m, __x); }

  /**
   * Return the associated Laguerre polynomial of order @c n,
   * degree @c m: @f$ L_n^m(x) @f$.
   *
   * The associated Laguerre function of real degree @f$ \alpha @f$,
   * @f$ L_n^\alpha(x) @f$, is defined by
   * @f[
   * 	 L_n^\alpha(x) = \frac{(\alpha + 1)_n}{n!}
   * 			 {}_1F_1(-n; \alpha + 1; x)
   * @f]
   * where @f$ (\alpha)_n @f$ is the Pochhammer symbol and
   * @f$ {}_1F_1(a; c; x) @f$ is the confluent hypergeometric function.
   *
   * The associated Laguerre polynomial is defined for integral
   * degree @f$ \alpha = m @f$ by:
   * @f[
   * 	 L_n^m(x) = (-1)^m \frac{d^m}{dx^m} L_{n + m}(x)
   * @f]
   * where the Laguerre polynomial is defined by:
   * @f[
   * 	 L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @param __n The order of the Laguerre function.
   * @param __m The degree of the Laguerre function.
   * @param __x The argument of the Laguerre function.
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    assoc_laguerre(unsigned int __n, unsigned int __m, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__assoc_laguerre<__type>(__n, __m, __x);
    }

  // Associated Legendre functions

  /**
   * Return the associated Legendre function of degree @c l and order @c m
   * for @c float argument.
   *
   * @see assoc_legendre for more details.
   */
  inline float
  assoc_legendref(unsigned int __l, unsigned int __m, float __x)
  { return __detail::__assoc_legendre_p<float>(__l, __m, __x); }

  /**
   * Return the associated Legendre function of degree @c l and order @c m.
   *
   * @see assoc_legendre for more details.
   */
  inline long double
  assoc_legendrel(unsigned int __l, unsigned int __m, long double __x)
  { return __detail::__assoc_legendre_p<long double>(__l, __m, __x); }

  /**
   * Return the associated Legendre function of degree @c l and order @c m.
   *
   * The associated Legendre function is derived from the Legendre function
   * @f$ P_l(x) @f$ by the Rodrigues formula:
   * @f[
   *   P_l^m(x) = (1 - x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)
   * @f]
   *
   * @param  __l  The degree of the associated Legendre function.
   * 		@f$ l >= 0 @f$.
   * @param  __m  The order of the associated Legendre function.
   * 		@f$ m <= l @f$.
   * @param  __x  The argument of the associated Legendre function.
   * 		@f$ |x| <= 1 @f$.
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    assoc_legendre(unsigned int __l, unsigned int __m, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__assoc_legendre_p<__type>(__l, __m, __x);
    }

  // Beta functions

  /**
   * Return the beta function, \f$B(a,b)\f$, for @c float parameters @c a, @c b.
   *
   * @see beta for more details.
   */
  inline float
  betaf(float __a, float __b)
  { return __detail::__beta<float>(__a, __b); }

  /**
   * Return the beta function, \f$B(a,b)\f$, for long double
   * parameters @c a, @c b.
   *
   * @see beta for more details.
   */
  inline long double
  betal(long double __a, long double __b)
  { return __detail::__beta<long double>(__a, __b); }

  /**
   * Return the beta function, \f$B(a,b)\f$, for real parameters @c a, @c b.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   */
  template<typename _Tpa, typename _Tpb>
    inline typename __gnu_cxx::__promote_2<_Tpa, _Tpb>::__type
    beta(_Tpa __a, _Tpb __b)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpa, _Tpb>::__type __type;
      return __detail::__beta<__type>(__a, __b);
    }

  // Complete elliptic integrals of the first kind

  /**
   * Return the complete elliptic integral of the first kind @f$ E(k) @f$
   * for @c float modulus @c k.
   *
   * @see comp_ellint_1 for details.
   */
  inline float
  comp_ellint_1f(float __k)
  { return __detail::__comp_ellint_1<float>(__k); }

  /**
   * Return the complete elliptic integral of the first kind @f$ E(k) @f$
   * for long double modulus @c k.
   *
   * @see comp_ellint_1 for details.
   */
  inline long double
  comp_ellint_1l(long double __k)
  { return __detail::__comp_ellint_1<long double>(__k); }

  /**
   * Return the complete elliptic integral of the first kind
   * @f$ K(k) @f$ for real modulus @c k.
   *
   * The complete elliptic integral of the first kind is defined as
   * @f[
   *   K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   * 					     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   * first kind.
   *
   * @param  __k  The modulus
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    comp_ellint_1(_Tp __k)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__comp_ellint_1<__type>(__k);
    }

  // Complete elliptic integrals of the second kind

  /**
   * Return the complete elliptic integral of the second kind @f$ E(k) @f$
   * for @c float modulus @c k.
   *
   * @see comp_ellint_2 for details.
   */
  inline float
  comp_ellint_2f(float __k)
  { return __detail::__comp_ellint_2<float>(__k); }

  /**
   * Return the complete elliptic integral of the second kind @f$ E(k) @f$
   * for long double modulus @c k.
   *
   * @see comp_ellint_2 for details.
   */
  inline long double
  comp_ellint_2l(long double __k)
  { return __detail::__comp_ellint_2<long double>(__k); }

  /**
   * Return the complete elliptic integral of the second kind @f$ E(k) @f$
   * for real modulus @c k.
   *
   * The complete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   * @f]
   *
   * @param  __k  The modulus
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    comp_ellint_2(_Tp __k)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__comp_ellint_2<__type>(__k);
    }

  // Complete elliptic integrals of the third kind

  /**
   * @brief Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) @f$ for @c float modulus @c k.
   *
   * @see comp_ellint_3 for details.
   */
  inline float
  comp_ellint_3f(float __k, float __nu)
  { return __detail::__comp_ellint_3<float>(__k, __nu); }

  /**
   * @brief Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) @f$ for @c long double modulus @c k.
   *
   * @see comp_ellint_3 for details.
   */
  inline long double
  comp_ellint_3l(long double __k, long double __nu)
  { return __detail::__comp_ellint_3<long double>(__k, __nu); }

  /**
   * Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ for real modulus @c k.
   *
   * The complete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu) = \int_0^{\pi/2}
   * 		     \frac{d\theta}
   * 		   {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   *
   * @param  __k  The modulus of the elliptic function.
   * @param  __nu  The argument of the elliptic function.
   */
  template<typename _Tp, typename _Tpn>
    inline typename __gnu_cxx::__promote_2<_Tp, _Tpn>::__type
    comp_ellint_3(_Tp __k, _Tpn __nu)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Tpn>::__type __type;
      return __detail::__comp_ellint_3<__type>(__k, __nu);
    }

  // Regular modified cylindrical Bessel functions

  /**
   * Return the regular modified Bessel function \f$ I_{\nu}(x) \f$
   * of @c float order \f$ \nu \f$ and argument f$ x \f$.
   *
   * @see cyl_bessel_i for setails.
   */
  inline float
  cyl_bessel_if(float __nu, float __x)
  { return __detail::__cyl_bessel_i<float>(__nu, __x); }

  /**
   * Return the regular modified Bessel function \f$ I_{\nu}(x) \f$
   * of @c long double order \f$ \nu \f$ and argument f$ x \f$.
   *
   * @see cyl_bessel_i for setails.
   */
  inline long double
  cyl_bessel_il(long double __nu, long double __x)
  { return __detail::__cyl_bessel_i<long double>(__nu, __x); }

  /**
   * Return the regular modified Bessel function \f$ I_{\nu}(x) \f$
   * of real order \f$ \nu \f$ and argument f$ x \f$.
   *
   * The regular modified cylindrical Bessel function is:
   * @f[
   *  I_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{(x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @param  __nu  The order of the regular modified Bessel function.
   * @param  __x   The argument of the regular modified Bessel function.
   */
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_bessel_i(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_bessel_i<__type>(__nu, __x);
    }

  // Cylindrical Bessel functions (of the first kind)

  /**
   * Return the Bessel function of the first kind \f$ J_{\nu}(x) \f$
   * of @c float order \f$ \nu \f$ and argument \f$ x \f$.
   *
   * @see cyl_bessel_j for setails.
   */
  inline float
  cyl_bessel_jf(float __nu, float __x)
  { return __detail::__cyl_bessel_j<float>(__nu, __x); }

  /**
   * Return the Bessel function of the first kind \f$ J_{\nu}(x) \f$
   * of @c long double order \f$ \nu \f$ and argument \f$ x \f$.
   *
   * @see cyl_bessel_j for setails.
   */
  inline long double
  cyl_bessel_jl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_j<long double>(__nu, __x); }

  /**
   * Return the Bessel function @f$ J_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x @f$.
   *
   * The cylindrical Bessel function is:
   * @f[
   *    J_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @param  __nu  The order of the Bessel function.
   * @param  __x   The argument of the Bessel function.
   */
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_bessel_j(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_bessel_j<__type>(__nu, __x);
    }

  // Irregular modified cylindrical Bessel functions

  /**
   * Return the irregular modified Bessel function \f$ K_{\nu}(x) \f$
   * of @c float order \f$ \nu \f$ for @c and argument \f$ x \f$.
   *
   * @see cyl_bessel_k for setails.
   */
  inline float
  cyl_bessel_kf(float __nu, float __x)
  { return __detail::__cyl_bessel_k<float>(__nu, __x); }

  /**
   * Return the irregular modified Bessel function \f$ K_{\nu}(x) \f$
   * of @c long double order \f$ \nu \f$ for @c and argument \f$ x \f$.
   *
   * @see cyl_bessel_k for setails.
   */
  inline long double
  cyl_bessel_kl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_k<long double>(__nu, __x); }

  /**
   * Return the irregular modified Bessel function \f$ K_{\nu}(x) \f$
   * of real order \f$ \nu \f$ and argument \f$ x \f$.
   *
   * The irregular modified Bessel function is defined by:
   * @f[
   * 	K_{\nu}(x) = \frac{\pi}{2}
   * 		     \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin \nu\pi}
   * @f]
   * where for integral \f$ \nu = n \f$ a limit is taken:
   * \f$ lim_{\nu \to n} \f$.
   * For negative argument we have simply:
   * @f[
   * 	K_{-\nu}(x) = K_{\nu}(x)
   * @f]
   *
   * @param  __nu  The order of the irregular modified Bessel function.
   * @param  __x   The argument of the irregular modified Bessel function.
   */
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_bessel_k(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_bessel_k<__type>(__nu, __x);
    }

  // Cylindrical Neumann functions

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of @c float order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * @see cyl_neumann for setails.
   */
  inline float
  cyl_neumannf(float __nu, float __x)
  { return __detail::__cyl_neumann_n<float>(__nu, __x); }

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of @c long double order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * @see cyl_neumann for setails.
   */
  inline long double
  cyl_neumannl(long double __nu, long double __x)
  { return __detail::__cyl_neumann_n<long double>(__nu, __x); }

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * The Neumann function is defined by:
   * @f[
   *    N_{\nu}(x) = \frac{J_{\nu}(x) \cos \nu\pi - J_{-\nu}(x)}
   *                      {\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   *
   * @param  __nu  The order of the Neumann function.
   * @param  __x   The argument of the Neumann function.
   */
  template<typename _Tpnu, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type
    cyl_neumann(_Tpnu __nu, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type __type;
      return __detail::__cyl_neumann_n<__type>(__nu, __x);
    }

  // Incomplete elliptic integrals of the first kind

  /**
   * Return the incomplete elliptic integral of the first kind @f$ E(k,\phi) @f$
   * for @c float modulus @f$ k @f$ and angle @f$ \phi @f$.
   *
   * @see ellint_1 for details.
   */
  inline float
  ellint_1f(float __k, float __phi)
  { return __detail::__ellint_1<float>(__k, __phi); }

  /**
   * Return the incomplete elliptic integral of the first kind @f$ E(k,\phi) @f$
   * for @c long double modulus @f$ k @f$ and angle @f$ \phi @f$.
   *
   * @see ellint_1 for details.
   */
  inline long double
  ellint_1l(long double __k, long double __phi)
  { return __detail::__ellint_1<long double>(__k, __phi); }

  /**
   * Return the incomplete elliptic integral of the first kind @f$ F(k,\phi) @f$
   * for @c real modulus @f$ k @f$ and angle @f$ \phi @f$.
   *
   * The incomplete elliptic integral of the first kind is defined as
   * @f[
   *   F(k,\phi) = \int_0^{\phi}\frac{d\theta}
   * 				     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * For  @f$ \phi= \pi/2 @f$ this becomes the complete elliptic integral of
   * the first kind, @f$ K(k) @f$.  @see comp_ellint_1.
   *
   * @param  __k  The modulus of the elliptic function.
   * @param  __phi  The integral limit argument of the elliptic function.
   */
  template<typename _Tp, typename _Tpp>
    inline typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type
    ellint_1(_Tp __k, _Tpp __phi)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type __type;
      return __detail::__ellint_1<__type>(__k, __phi);
    }

  // Incomplete elliptic integrals of the second kind

  /**
   * @brief Return the incomplete elliptic integral of the second kind
   * @f$ E(k,\phi) @f$ for @c float argument.
   *
   * @see ellint_2 for details.
   */
  inline float
  ellint_2f(float __k, float __phi)
  { return __detail::__ellint_2<float>(__k, __phi); }

  /**
   * @brief Return the incomplete elliptic integral of the second kind
   * @f$ E(k,\phi) @f$.
   *
   * @see ellint_2 for details.
   */
  inline long double
  ellint_2l(long double __k, long double __phi)
  { return __detail::__ellint_2<long double>(__k, __phi); }

  /**
   * Return the incomplete elliptic integral of the second kind
   * @f$ E(k,\phi) @f$.
   *
   * The incomplete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\phi) = \int_0^{\phi} \sqrt{1 - k^2 sin^2\theta}
   * @f]
   * For  @f$ \phi= \pi/2 @f$ this becomes the complete elliptic integral of
   * the second kind, @f$ E(k) @f$.  @see comp_ellint_2.
   *
   * @param  __k  The argument of the elliptic function.
   * @param  __phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the second kind.
   */
  template<typename _Tp, typename _Tpp>
    inline typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type
    ellint_2(_Tp __k, _Tpp __phi)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Tpp>::__type __type;
      return __detail::__ellint_2<__type>(__k, __phi);
    }

  // Incomplete elliptic integrals of the third kind

  /**
   * @brief Return the incomplete elliptic integral of the third kind
   * @f$ \Pi(k,\nu,\phi) @f$ for @c float argument.
   *
   * @see ellint_3 for details.
   */
  inline float
  ellint_3f(float __k, float __nu, float __phi)
  { return __detail::__ellint_3<float>(__k, __nu, __phi); }

  /**
   * @brief Return the incomplete elliptic integral of the third kind
   * @f$ \Pi(k,\nu,\phi) @f$.
   *
   * @see ellint_3 for details.
   */
  inline long double
  ellint_3l(long double __k, long double __nu, long double __phi)
  { return __detail::__ellint_3<long double>(__k, __nu, __phi); }

  /**
   * @brief Return the incomplete elliptic integral of the third kind
   * @f$ \Pi(k,\nu,\phi) @f$.
   *
   * The incomplete elliptic integral of the third kind is defined by:
   * @f[
   *   \Pi(k,\nu,\phi) = \int_0^{\phi}
   * 			 \frac{d\theta}
   * 			 {(1 - \nu \sin^2\theta)
   * 			  \sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   * For  @f$ \phi= \pi/2 @f$ this becomes the complete elliptic integral of
   * the third kind, @f$ \Pi(k,\nu) @f$.  @see comp_ellint_3.
   *
   * @param  __k  The modulus of the elliptic function.
   * @param  __nu  The second argument of the elliptic function.
   * @param  __phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the third kind.
   */
  template<typename _Tp, typename _Tpn, typename _Tpp>
    inline typename __gnu_cxx::__promote_3<_Tp, _Tpn, _Tpp>::__type
    ellint_3(_Tp __k, _Tpn __nu, _Tpp __phi)
    {
      typedef typename __gnu_cxx::__promote_3<_Tp, _Tpn, _Tpp>::__type __type;
      return __detail::__ellint_3<__type>(__k, __nu, __phi);
    }

  // Exponential integrals

  /**
   * Return the exponential integral @f$ Ei(x) @f$ for @c float argument @c x.
   *
   * @see expint for details.
   */
  inline float
  expintf(float __x)
  { return __detail::__expint<float>(__x); }

  /**
   * Return the exponential integral @f$ Ei(x) @f$
   * for @c long double argument @c x.
   *
   * @see expint for details.
   */
  inline long double
  expintl(long double __x)
  { return __detail::__expint<long double>(__x); }

  /**
   * Return the exponential integral @f$ Ei(x) @f$ for @c lreal argument @c x.
   *
   * The exponential integral is given by
   * \f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * \f]
   *
   * @param  __x  The argument of the exponential integral function.
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    expint(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__expint<__type>(__x);
    }

  // Hermite polynomials

  /**
   * Return the Hermite polynomial of order n, @f$ H_n(x) @f$,
   * for float argument @c x.
   *
   * @see hermite for details.
   */
  inline float
  hermitef(unsigned int __n, float __x)
  { return __detail::__poly_hermite<float>(__n, __x); }

  /**
   * Return the Hermite polynomial of order n, @f$ H_n(x) @f$,
   * for @c long double argument @c x.
   *
   * @see hermite for details.
   */
  inline long double
  hermitel(unsigned int __n, long double __x)
  { return __detail::__poly_hermite<long double>(__n, __x); }

  /**
   * Return the Hermite polynomial of order n, @f$ H_n(x) @f$,
   * for @c real argument @c x.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * The Hermite polynomial obeys a reflection formula:
   * @f[
   *   H_n(-x) = (-1)^n H_n(x)
   * @f]
   *
   * @param __n The order
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    hermite(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__poly_hermite<__type>(__n, __x);
    }

  // Laguerre polynomials

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$
   * of degree @c n and @c float argument @c x.
   *
   * @see laguerre for more details.
   */
  inline float
  laguerref(unsigned int __n, float __x)
  { return __detail::__laguerre<float>(__n, __x); }

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$
   * of degree @c n and @c long double argument @c x.
   *
   * @see laguerre for more details.
   */
  inline long double
  laguerrel(unsigned int __n, long double __x)
  { return __detail::__laguerre<long double>(__n, __x); }

  /**
   * Returns the Laguerre polynomial
   * of degree @c n, and argument @c x: @f$ L_n(x) @f$.
   *
   * The Laguerre polynomial is defined by:
   * @f[
   * 	 L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @param __n The order of the Laguerre function.
   * @param __x The argument of the Laguerre function.
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    laguerre(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__laguerre<__type>(__n, __x);
    }

  // Legendre polynomials

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of degree @f$ l @f$
   * for @c float argument.
   *
   * @see legendre for more details.
   */
  inline float
  legendref(unsigned int __l, float __x)
  { return __detail::__poly_legendre_p<float>(__l, __x); }

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of degree @f$ l @f$
   * for @c long double argument.
   *
   * @see legendre for more details.
   */
  inline long double
  legendrel(unsigned int __l, long double __x)
  { return __detail::__poly_legendre_p<long double>(__l, __x); }

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of degree @f$ l @f$
   * for @c real argument.
   *
   * The Legendre function of order @f$ l @f$ and argument @f$ x @f$,
   * @f$ P_l(x) @f$, is defined by:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   *
   * @param  __l  The degree @f$l >= 0@f$
   * @param  __x  The argument @f$|x| <= 1@f$
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    legendre(unsigned int __l, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__poly_legendre_p<__type>(__l, __x);
    }

  // Riemann zeta functions

  /**
   * Return the Riemann zeta function @f$ \zeta(s) @f$
   * for @c float argument @f$ s @f$.
   *
   * @see riemann_zeta for more details.
   */
  inline float
  riemann_zetaf(float __s)
  { return __detail::__riemann_zeta<float>(__s); }

  /**
   * Return the Riemann zeta function @f$ \zeta(s) @f$
   * for @c long double argument @f$ s @f$.
   *
   * @see riemann_zeta for more details.
   */
  inline long double
  riemann_zetal(long double __s)
  { return __detail::__riemann_zeta<long double>(__s); }

  /**
   * Return the Riemann zeta function @f$ \zeta(s) @f$
   * for real argument @f$ s @f$.
   *
   * The Riemann zeta function is defined by:
   * @f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} k^{-s} for s > 1
   * 		   \frac{(2\pi)^s}{pi} sin(\frac{\pi s}{2})
   * 		   \Gamma (1 - s) \zeta (1 - s) for s < 1
   * @f]
   * For s < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   * @f]
   *
   * @param __s The argument @c s != 1
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    riemann_zeta(_Tp __s)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__riemann_zeta<__type>(__s);
    }

  // Spherical Bessel functions

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of order n
   * for @c float argument.
   *
   * @see sph_bessel for more details.
   */
  inline float
  sph_besself(unsigned int __n, float __x)
  { return __detail::__sph_bessel<float>(__n, __x); }

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of order n
   * for @c long double argument.
   *
   * @see sph_bessel for more details.
   */
  inline long double
  sph_bessell(unsigned int __n, long double __x)
  { return __detail::__sph_bessel<long double>(__n, __x); }

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of order n
   * for real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  j_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} J_{n+1/2}(x)
   * @f]
   *
   * @param  __n  The non-negative integral order
   * @param  __x  The non-negative real argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_bessel(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_bessel<__type>(__n, __x);
    }

  // Spherical associated Legendre functions

  /**
   * Return the spherical Legendre function of non-negative integral
   * degree @c l and order @c m and float angle @f$ \theta @f$ in radians.
   *
   * @see sph_legendre for details.
   */
  inline float
  sph_legendref(unsigned int __l, unsigned int __m, float __theta)
  { return __detail::__sph_legendre<float>(__l, __m, __theta); }

  /**
   * Return the spherical Legendre function of non-negative integral
   * degree @c l and order @c m and @c long double angle @f$ \theta @f$
   * in radians.
   *
   * @see sph_legendre for details.
   */
  inline long double
  sph_legendrel(unsigned int __l, unsigned int __m, long double __theta)
  { return __detail::__sph_legendre<long double>(__l, __m, __theta); }

  /**
   * Return the spherical Legendre function of non-negative integral
   * degree @c l and order @c m and real angle @f$ \theta @f$ in radians.
   *
   * The spherical Legendre function is defined by
   * @f[
   *  Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   *                              \frac{(l-m)!}{(l+m)!}]
   *                   P_l^m(\cos\theta) \exp^{im\phi}
   * @f]
   *
   * @param __l The non-negative order @f$ l >= 0 @f$.
   * @param __m The non-negative degree @f$ m >= 0 @f$ and @f$ m <= l @f$.
   * @param __theta The radian polar angle argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_legendre(unsigned int __l, unsigned int __m, _Tp __theta)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_legendre<__type>(__l, __m, __theta);
    }

  // Spherical Neumann functions

  /**
   * Return the spherical Neumann function of non-negative integral order @c n
   * and non-negative @c float argument @c x.
   *
   * @see sph_neumann for details.
   */
  inline float
  sph_neumannf(unsigned int __n, float __x)
  { return __detail::__sph_neumann<float>(__n, __x); }

  /**
   * Return the spherical Neumann function of non-negative integral order @c n
   * and non-negative @c long double argument @c x.
   *
   * @see sph_neumann for details.
   */
  inline long double
  sph_neumannl(unsigned int __n, long double __x)
  { return __detail::__sph_neumann<long double>(__n, __x); }

  /**
   * Return the spherical Neumann function of non-negative integral order @c n
   * and non-negative real argument @c x.
   *
   * The spherical Neumann function is defined by
   * @f[
   *    n_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} N_{n+1/2}(x)
   * @f]
   *
   * @param  __n  The non-negative integral order
   * @param  __x  The non-negative real argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_neumann(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_neumann<__type>(__n, __x);
    }

  /* @} */ // tr29124_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
  /**
   * @defgroup gnu_math_spec_func Extended Mathematical Special Functions
   * @ingroup numerics
   *
   * A collection of advanced mathematical special functions.
   * @{
   */

  // Confluent hypergeometric functions

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of @c float numeratorial parameter @c a, denominatorial parameter @c c,
   * and argument @c x.
   *
   * @see conf_hyperg for details.
   */
  inline float
  conf_hypergf(float __a, float __c, float __x)
  { return std::__detail::__conf_hyperg<float>(__a, __c, __x); }

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of @c long double numeratorial parameter @c a,
   * denominatorial parameter @c c, and argument @c x.
   *
   * @see conf_hyperg for details.
   */
  inline long double
  conf_hypergl(long double __a, long double __c, long double __x)
  { return std::__detail::__conf_hyperg<long double>(__a, __c, __x); }

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of real numeratorial parameter @c a, denominatorial parameter @c c,
   * and argument @c x.
   *
   * The confluent hypergeometric function is defined by
   * @f[
   *    {}_1F_1(a;c;x) = \sum_{n=0}^{\infty} \frac{(a)_n x^n}{(c)_n n!}
   * @f]
   * where the Pochhammer symbol is @f$ (x)_k = (x)(x+1)...(x+k-1) @f$,
   * @f$ (x)_0 = 1 @f$
   *
   * @param __a The numeratorial parameter
   * @param __c The denominatorial parameter
   * @param __x The argument
   */
  template<typename _Tpa, typename _Tpc, typename _Tp>
    inline typename __gnu_cxx::__promote_3<_Tpa, _Tpc, _Tp>::__type
    conf_hyperg(_Tpa __a, _Tpc __c, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_3<_Tpa, _Tpc, _Tp>::__type __type;
      return std::__detail::__conf_hyperg<__type>(__a, __c, __x);
    }

  // Hypergeometric functions

  /**
   * Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * of @ float numeratorial parameters @c a and @c b,
   * denominatorial parameter @c c, and argument @c x.
   *
   * @see hyperg for details.
   */
  inline float
  hypergf(float __a, float __b, float __c, float __x)
  { return std::__detail::__hyperg<float>(__a, __b, __c, __x); }

  /**
   * Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * of @ long double numeratorial parameters @c a and @c b,
   * denominatorial parameter @c c, and argument @c x.
   *
   * @see hyperg for details.
   */
  inline long double
  hypergl(long double __a, long double __b, long double __c, long double __x)
  { return std::__detail::__hyperg<long double>(__a, __b, __c, __x); }

  /**
   * Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * of real numeratorial parameters @c a and @c b,
   * denominatorial parameter @c c, and argument @c x.
   *
   * The hypergeometric function is defined by
   * @f[
   *    {}_2F_1(a;c;x) = \sum_{n=0}^{\infty} \frac{(a)_n (b)_n x^n}{(c)_n n!}
   * @f]
   * where the Pochhammer symbol is @f$ (x)_k = (x)(x+1)...(x+k-1) @f$,
   * @f$ (x)_0 = 1 @f$
   *
   * @param __a The first numeratorial parameter
   * @param __b The second numeratorial parameter
   * @param __c The denominatorial parameter
   * @param __x The argument
   */
  template<typename _Tpa, typename _Tpb, typename _Tpc, typename _Tp>
    inline typename __gnu_cxx::__promote_4<_Tpa, _Tpb, _Tpc, _Tp>::__type
    hyperg(_Tpa __a, _Tpb __b, _Tpc __c, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_4<_Tpa, _Tpb, _Tpc, _Tp>
		::__type __type;
      return std::__detail::__hyperg<__type>(__a, __b, __c, __x);
    }

#if __cplusplus >= 201103L

  // Confluent hypergeometric limit functions

  /**
   * Return the confluent hypergeometric limit function @f$ {}_0F_1(;c;x) @f$
   * of @c float numeratorial parameter @c c and argument @c x.
   *
   * @see conf_hyperg_lim for details.
   */
  inline float
  conf_hyperg_limf(float __c, float __x)
  { return std::__detail::__conf_hyperg_lim<float>(__c, __x); }

  /**
   * Return the confluent hypergeometric limit function @f$ {}_0F_1(;c;x) @f$
   * of @c long double numeratorial parameter @c c and argument @c x.
   *
   * @see conf_hyperg_lim for details.
   */
  inline long double
  conf_hyperg_liml(long double __c, long double __x)
  { return std::__detail::__conf_hyperg_lim<long double>(__c, __x); }

  /**
   * Return the confluent hypergeometric limit function @f$ {}_0F_1(;c;x) @f$
   * of real numeratorial parameter @c c and argument @c x.
   *
   * The confluent hypergeometric limit function is defined by
   * @f[
   *    {}_0F_1(;c;x) = \sum_{n=0}^{\infty} \frac{x^n}{(c)_n n!}
   * @f]
   * where the Pochhammer symbol is @f$ (x)_k = (x)(x+1)...(x+k-1) @f$,
   * @f$ (x)_0 = 1 @f$
   *
   * @param __c The denominatorial parameter
   * @param __x The argument
   */
  template<typename _Tpc, typename _Tp>
    inline typename __gnu_cxx::__promote_2<_Tpc, _Tp>::__type
    conf_hyperg_lim(_Tpc __c, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote_2<_Tpc, _Tp>::__type __type;
      return std::__detail::__conf_hyperg_lim<__type>(__c, __x);
    }

  // Unnormalized sinus cardinal functions

  inline float
  sinc_pif(float __x)
  { return std::__detail::__sinc_pi<float>(__x); }

  inline long double
  sinc_pil(long double __x)
  { return std::__detail::__sinc_pi<long double>(__x); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinc_pi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sinc_pi<__type>(__x);
    }

  // Normalized sinus cardinal functions

  inline float
  sincf(float __x)
  { return std::__detail::__sinc<float>(__x); }

  inline long double
  sincl(long double __x)
  { return std::__detail::__sinc<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinc(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sinc<__type>(__x);
    }

  // Logarithmic integrals

  /**
   * Return the logarithmic integral of argument @c x.
   *
   * @see logint for details.
   */
  inline float
  logintf(float __x)
  { return std::__detail::__logint<float>(__x); }

  /**
   * Return the logarithmic integral of argument @c x.
   *
   * @see logint for details.
   */
  inline long double
  logintl(long double __x)
  { return std::__detail::__logint<long double>(__x); }

  /**
   * Return the logarithmic integral of argument @c x.
   *
   * The logarithmic integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    logint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__logint<__type>(__x);
    }

  // Sine integrals

  /**
   * Return the sine integral of argument @c x.
   *
   * @see sinint for details.
   */
  inline float
  sinintf(float __x)
  { return std::__detail::__sincosint<float>(__x).first; }

  /**
   * Return the sine integral of argument @c x.
   *
   * @see sinint for details.
   */
  inline long double
  sinintl(long double __x)
  { return std::__detail::__sincosint<long double>(__x).first; }

  /**
   * Return the sine integral of argument @c x.
   *
   * The sine integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sincosint<__type>(__x).first;
    }

  // Cosine integrals

  /**
   * Return the cosine integral of argument @c x.
   *
   * @see cosint for details.
   */
  inline float
  cosintf(float __x)
  { return std::__detail::__sincosint<float>(__x).second; }

  /**
   * Return the cosine integral of argument @c x.
   *
   * @see cosint for details.
   */
  inline long double
  cosintl(long double __x)
  { return std::__detail::__sincosint<long double>(__x).second; }

  /**
   * Return the cosine integral of argument @c x.
   *
   * The cosine integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    cosint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sincosint<__type>(__x).second;
    }

  // Hyperbolic sine integrals

  /**
   * Return the hyperbolic sine integral of argument @c x.
   *
   * @see sinhint for details.
   */
  inline float
  sinhintf(float __x)
  { return std::__detail::__sinhint<float>(__x); }

  /**
   * Return the hyperbolic sine integral of argument @c x.
   *
   * @see sinhint for details.
   */
  inline long double
  sinhintl(long double __x)
  { return std::__detail::__sinhint<long double>(__x); }

  /**
   * Return the hyperbolic sine integral of argument @c x.
   *
   * The sine hyperbolic integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinhint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sinhint<__type>(__x);
    }

  // Hyperbolic cosine integrals

  inline float
  coshintf(float __x)
  { return std::__detail::__coshint<float>(__x); }

  /**
   * Return the hyperbolic cosine integral of argument @c x.
   *
   * @see coshint for details.
   */
  inline long double
  coshintl(long double __x)
  { return std::__detail::__coshint<long double>(__x); }

  /**
   * Return the hyperbolic cosine integral of argument @c x.
   *
   * The hyperbolic cosine integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    coshint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__coshint<__type>(__x);
    }

  // Slots for Jacobi elliptic function tuple.
  enum
  {
    _GLIBCXX_JACOBI_SN,
    _GLIBCXX_JACOBI_CN,
    _GLIBCXX_JACOBI_DN
  };

  // Jacobi elliptic sn functions.

  /**
   * Return the Jacobi elliptic @c sn integral of modulus @c k and argument @c u.
   *
   * @see jacobi_sn for details.
   */
  inline float
  jacobi_snf(float __k, float __u)
  {
    return std::get<_GLIBCXX_JACOBI_SN>
		(std::__detail::__jacobi_sncndn<float>(__k, __u));
  }

  /**
   * Return the Jacobi elliptic @c sn integral of modulus @c k and argument @c u.
   *
   * @see jacobi_sn for details.
   */
  inline long double
  jacobi_snl(long double __k, long double __u)
  {
    return std::get<_GLIBCXX_JACOBI_SN>
		(std::__detail::__jacobi_sncndn<long double>(__k, __u));
  }

  /**
   * Return the Jacobi elliptic @c sn integral of modulus @c k and argument @c u.
   *
   * The Jacobi elliptic @c sn integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus
   * @param __u The argument
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Kp, _Up>
    jacobi_sn(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_num_t<_Kp, _Up>;
      return std::get<_GLIBCXX_JACOBI_SN>
		(std::__detail::__jacobi_sncndn<__type>(__k, __u));
    }

  // Jacobi elliptic cn functions.

  /**
   * Return the Jacobi elliptic @c cn integral of modulus @c k and argument @c u.
   *
   * @see jacobi_cn for details.
   */
  inline float
  jacobi_cnf(float __k, float __u)
  {
    return std::get<_GLIBCXX_JACOBI_CN>
		(std::__detail::__jacobi_sncndn<float>(__k, __u));
  }

  /**
   * Return the Jacobi elliptic @c cn integral of modulus @c k and argument @c u.
   *
   * @see jacobi_cn for details.
   */
  inline long double
  jacobi_cnl(long double __k, long double __u)
  {
    return std::get<_GLIBCXX_JACOBI_CN>
		(std::__detail::__jacobi_sncndn<long double>(__k, __u));
  }

  /**
   * Return the Jacobi elliptic @c cn integral of modulus @c k and argument @c u.
   *
   * The Jacobi elliptic @c cn integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus
   * @param __u The argument
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Kp, _Up>
    jacobi_cn(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_num_t<_Kp, _Up>;
      return std::get<_GLIBCXX_JACOBI_CN>
		(std::__detail::__jacobi_sncndn<__type>(__k, __u));
    }

  // Jacobi elliptic dn functions.

  /**
   * Return the Jacobi elliptic @c dn integral of modulus @c k and argument @c u.
   *
   * @see jacobi_dn for details.
   */
  inline float
  jacobi_dnf(float __k, float __u)
  {
    return std::get<_GLIBCXX_JACOBI_DN>
		(std::__detail::__jacobi_sncndn<float>(__k, __u));
  }

  /**
   * Return the Jacobi elliptic @c dn integral of modulus @c k and argument @c u.
   *
   * @see jacobi_dn for details.
   */
  inline long double
  jacobi_dnl(long double __k, long double __u)
  {
    return std::get<_GLIBCXX_JACOBI_DN>
		(std::__detail::__jacobi_sncndn<long double>(__k, __u));
  }

  /**
   * Return the Jacobi elliptic @c dn integral of modulus @c k and argument @c u.
   *
   * The Jacobi elliptic @c dn integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus
   * @param __u The argument
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Kp, _Up>
    jacobi_dn(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_num_t<_Kp, _Up>;
      return std::get<_GLIBCXX_JACOBI_DN>
		(std::__detail::__jacobi_sncndn<__type>(__k, __u));
    }

  // Chebyshev polynomials of the first kind

  /**
   * Return the Chebyshev polynomials of the first kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_t for details.
   */
  inline float
  chebyshev_tf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_t<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the first kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_t for details.
   */
  inline long double
  chebyshev_tl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_t<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the first kind of order @c n
   * and argument @c x.
   *
   * The Chebyshev polynomials of the first kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __n
   * @param __x
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_t(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_t<__type>(__n, __x);
    }

  // Chebyshev polynomials of the second kind

  /**
   * Return the Chebyshev polynomials of the second kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_u for details.
   */
  inline float
  chebyshev_uf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_u<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the second kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_u for details.
   */
  inline long double
  chebyshev_ul(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_u<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the second kind of order @c n
   * and argument @c x.
   *
   * The Chebyshev polynomials of the second kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __n
   * @param __x
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_u(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_u<__type>(__n, __x);
    }

  // Chebyshev polynomials of the third kind

  /**
   * Return the Chebyshev polynomials of the third kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_v for details.
   */
  inline float
  chebyshev_vf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_v<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the third kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_v for details.
   */
  inline long double
  chebyshev_vl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_v<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the third kind of order @c n
   * and argument @c x.
   *
   * The Chebyshev polynomials of the third kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __n
   * @param __x
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_v(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_v<__type>(__n, __x);
    }

  // Chebyshev polynomials of the fourth kind

  /**
   * Return the Chebyshev polynomials of the fourth kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_w for details.
   */
  inline float
  chebyshev_wf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_w<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the fourth kind of order @c n
   * and argument @c x.
   *
   * @see chebyshev_w for details.
   */
  inline long double
  chebyshev_wl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_w<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the fourth kind of order @c n
   * and argument @c x.
   *
   * The Chebyshev polynomials of the fourth kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __n
   * @param __x
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    chebyshev_w(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__chebyshev_w<__type>(__n, __x);
    }

  // Jacobi polynomials

  inline float
  jacobif(unsigned __n, float __alpha, float __beta, float __x)
  { return std::__detail::__poly_jacobi<float>(__n, __alpha, __beta, __x); }

  inline long double
  jacobil(unsigned __n, long double __alpha, long double __beta, long double __x)
  { return std::__detail::__poly_jacobi<long double>(__n, __alpha, __beta, __x); }

  /**
   * 
   */
  template<typename _Talpha, typename _Tbeta, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Talpha, _Tbeta, _Tp>
    jacobi(unsigned __n, _Talpha __alpha, _Tbeta __beta, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Talpha, _Tbeta, _Tp>;
      return std::__detail::__poly_jacobi<__type>(__n, __alpha, __beta, __x);
    }

  // Gegenbauer polynomials

  inline float
  gegenbauerf(unsigned int __n, float __alpha, float __x)
  { return std::__detail::__gegenbauer_poly<float>(__n, __alpha, __x); }

  inline long double
  gegenbauerl(unsigned int __n, long double __alpha, long double __x)
  { return std::__detail::__gegenbauer_poly<long double>(__n, __alpha, __x); }

  /**
   * 
   */
  template<typename _Talpha, typename _Tp>
    inline typename __gnu_cxx::__promote_num_t<_Talpha, _Tp>
    gegenbauer(unsigned int __n, _Talpha __alpha, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Talpha, _Tp>;
      return std::__detail::__gegenbauer_poly<__type>(__n, __alpha, __x);
    }

  // Zernike polynomials

  inline float
  zernikef(unsigned int __n, int __m, float __rho, float __phi)
  { return std::__detail::__zernike<float>(__n, __m, __rho, __phi); }

  inline long double
  zernikel(unsigned int __n, int __m, long double __rho, long double __phi)
  { return std::__detail::__zernike<long double>(__n, __m, __rho, __phi); }

  /**
   * 
   */
  template<typename _Trho, typename _Tphi>
    inline __gnu_cxx::__promote_num_t<_Trho, _Tphi>
    zernike(unsigned int __n, int __m, _Trho __rho, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Trho, _Tphi>;
      return std::__detail::__zernike<__type>(__n, __m, __rho, __phi);
    }

  // Radial polynomials

  inline float
  radpolyf(unsigned int __n, unsigned int __m, float __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  inline long double
  radpolyl(unsigned int __n, unsigned int __m, long double __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    radpoly(unsigned int __n, unsigned int __m, _Tp __rho)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__poly_radial_jacobi<__type>(__n, __m, __rho);
    }

  // Unnormalized hyperbolic sinus cardinal functions

  inline float
  sinhc_pif(float __x)
  { return std::__detail::__sinhc_pi<float>(__x); }

  inline long double
  sinhc_pil(long double __x)
  { return std::__detail::__sinhc_pi<long double>(__x); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinhc_pi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sinhc_pi<__type>(__x);
    }

  // Normalized hyperbolic sinus cardinal functions

  inline float
  sinhcf(float __x)
  { return std::__detail::__sinhc<float>(__x); }

  inline long double
  sinhcl(long double __x)
  { return std::__detail::__sinhc<long double>(__x); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sinhc(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sinhc<__type>(__x);
    }

  // Cylindrical Hankel functions of the first kind

  inline std::complex<float>
  cyl_hankel_1f(float __nu, float __z)
  { return std::__detail::__cyl_hankel_1<float>(__nu, __z); }

  inline std::complex<long double>
  cyl_hankel_1l(long double __nu, long double __z)
  { return std::__detail::__cyl_hankel_1<long double>(__nu, __z); }

  /**
   * 
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tpnu, _Tp>>
    cyl_hankel_1(_Tpnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_1<__type>(__nu, __z);
    }

  // Cylindrical Hankel functions of the second kind

  inline std::complex<float>
  cyl_hankel_2f(float __nu, float __z)
  { return std::__detail::__cyl_hankel_2<float>(__nu, __z); }

  inline std::complex<long double>
  cyl_hankel_2l(long double __nu, long double __z)
  { return std::__detail::__cyl_hankel_2<long double>(__nu, __z); }

  /**
   * 
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tpnu, _Tp>>
    cyl_hankel_2(_Tpnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_2<__type>(__nu, __z);
    }

  // Spherical Hankel functions of the first kind

  inline std::complex<float>
  sph_hankel_1f(unsigned int __n, float __z)
  { return std::__detail::__sph_hankel_1<float>(__n, __z); }

  inline std::complex<long double>
  sph_hankel_1l(unsigned int __n, long double __z)
  { return std::__detail::__sph_hankel_1<long double>(__n, __z); }

  /**
   * 
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    sph_hankel_1(unsigned int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sph_hankel_1<__type>(__n, __z);
    }

  // Spherical Hankel functions of the second kind

  inline std::complex<float>
  sph_hankel_2f(unsigned int __n, float __z)
  { return std::__detail::__sph_hankel_2<float>(__n, __z); }

  inline std::complex<long double>
  sph_hankel_2l(unsigned int __n, long double __z)
  { return std::__detail::__sph_hankel_2<long double>(__n, __z); }

  /**
   * 
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    sph_hankel_2(unsigned int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sph_hankel_2<__type>(__n, __z);
    }

  // Modified spherical Bessel functions of the first kind

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

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    sph_bessel_i(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __i_n, __k_n, __ip_n, __kp_n;
      std::__detail::__sph_bessel_ik<__type>(__n, __x,
        				     __i_n, __k_n, __ip_n, __kp_n);
      return __i_n;
    }

  // Modified spherical Bessel functions of the second kind

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

  /**
   * 
   */
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

  // Airy functions of the first kind

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

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    airy_ai(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __Ai, __Bi, __Aip, __Bip;
      std::__detail::__airy<__type>(__x, __Ai, __Bi, __Aip, __Bip);
      return __Ai;
    }

  // Airy functions of the second kind

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

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    airy_bi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      __type __Ai, __Bi, __Aip, __Bip;
      std::__detail::__airy<__type>(__x, __Ai, __Bi, __Aip, __Bip);
      return __Bi;
    }

  // Upper incomplete gamma functions

  inline float
  gamma_uf(float __n, float __x)
  { return std::__detail::__gamma_u<float>(__n, __x); }

  inline long double
  gamma_ul(long double __n, long double __x)
  { return std::__detail::__gamma_u<long double>(__n, __x); }

  /**
   * 
   */
  template<typename _Tn, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tn, _Tp>
    gamma_u(_Tn __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tn, _Tp>;
      return std::__detail::__gamma_u<__type>(__n, __x);
    }

  // Lower incomplete gamma functions

  inline float
  gamma_lf(float __n, float __x)
  { return std::__detail::__gamma_l<float>(__n, __x); }

  inline long double
  gamma_ll(long double __n, long double __x)
  { return std::__detail::__gamma_l<long double>(__n, __x); }

  /**
   * 
   */
  template<typename _Tn, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tn, _Tp>
    gamma_l(_Tn __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tn, _Tp>;
      return std::__detail::__gamma_l<__type>(__n, __x);
    }

  // Digamma functions

  inline float
  digammaf(float __z)
  { return std::__detail::__psi<float>(__z); }

  inline long double
  digammal(long double __z)
  { return std::__detail::__psi<long double>(__z); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    digamma(_Tp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__psi<__type>(__z);
    }

  // Dilogarithm functions

  inline float
  dilogf(float __x)
  { return std::__detail::__dilog<float>(__x); }

  inline long double
  dilogl(long double __x)
  { return std::__detail::__dilog<long double>(__x); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    dilog(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__dilog<__type>(__x);
    }

  // Complete Carlson elliptic R_F functions

  inline float
  comp_ellint_rf(float __x, float __y)
  { return std::__detail::__comp_ellint_rf<float>(__x, __y); }

  inline long double
  comp_ellint_rf(long double __x, long double __y)
  { return std::__detail::__comp_ellint_rf<long double>(__x, __y); }

  /**
   * 
   */
  template<typename _Tx, typename _Ty>
    inline __gnu_cxx::__promote_num_t<_Tx, _Ty>
    comp_ellint_rf(_Tx __x, _Ty __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tx, _Ty>;
      return std::__detail::__comp_ellint_rf<__type>(__x, __y);
    }

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$.
   *
   * @see ellint_rf for details.
   */
  // Carlson elliptic R_F functions

  inline float
  ellint_rff(float __x, float __y, float __z)
  { return std::__detail::__ellint_rf<float>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$.
   *
   * @see ellint_rf for details.
   */
  inline long double
  ellint_rfl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rf<long double>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * of the first kind.
   *
   * The Carlson elliptic function of the first kind is defined by:
   * @f[
   *    R_F(x,y,z) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}}
   * @f]
   *
   * @param  __x  The first of three symmetric arguments.
   * @param  __y  The second of three symmetric arguments.
   * @param  __z  The third of three symmetric arguments.
   */
  template<typename _Tp, typename _Up, typename _Vp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>
    ellint_rf(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rf<__type>(__x, __y, __z);
    }

  // Carlson elliptic R_C functions

  /**
   * Return the Carlson elliptic function @f$ R_C(x,y) @f$.
   *
   * @see ellint_rc for details.
   */
  inline float
  ellint_rcf(float __x, float __y)
  { return std::__detail::__ellint_rc<float>(__x, __y); }

  /**
   * Return the Carlson elliptic function @f$ R_C(x,y) @f$.
   *
   * @see ellint_rc for details.
   */
  inline long double
  ellint_rcl(long double __x, long double __y)
  { return std::__detail::__ellint_rc<long double>(__x, __y); }

  /**
   * Return the Carlson elliptic function @f$ R_C(x,y) = R_F(x,y,y) @f$
   * where @f$ R_F(x,y,z) @f$ is the Carlson elliptic function
   * of the first kind.
   *
   * The Carlson elliptic function is defined by:
   * @f[
   *    R_C(x,y) = \frac{1}{2} \int_0^\infty
   *               \frac{dt}{(t + x)^{1/2}(t + y)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *    by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first argument.
   * @param  __y  The second argument.
   */
  template<typename _Tp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up>
    ellint_rc(_Tp __x, _Up __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up>;
      return std::__detail::__ellint_rc<__type>(__x, __y);
    }

  // Carlson elliptic R_J functions

  /**
   * Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$.
   *
   * @see ellint_rj for details.
   */
  inline float
  ellint_rjf(float __x, float __y, float __z, float __p)
  { return std::__detail::__ellint_rj<float>(__x, __y, __z, __p); }

  /**
   * Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$.
   *
   * @see ellint_rj for details.
   */
  inline long double
  ellint_rjl(long double __x, long double __y, long double __z, long double __p)
  { return std::__detail::__ellint_rj<long double>(__x, __y, __z, __p); }

  /**
   * Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$
   * 	     of the third kind.
   *
   * The Carlson elliptic function of the third kind is defined by:
   * @f[
   *    R_J(x,y,z,p) = \frac{3}{2} \int_0^\infty
   *                   \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}(t + p)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first of three symmetric arguments.
   * @param  __y  The second of three symmetric arguments.
   * @param  __z  The third of three symmetric arguments.
   * @param  __p  The fourth argument.
   */
  template<typename _Tp, typename _Up, typename _Vp, typename _Wp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp, _Wp>
    ellint_rj(_Tp __x, _Up __y, _Vp __z, _Wp __p)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp, _Wp>;
      return std::__detail::__ellint_rj<__type>(__x, __y, __z, __p);
    }

  // Carlson elliptic R_D functions

  /**
   * Return the Carlson elliptic function @f$ R_D(x,y,z) @f$.
   *
   * @see ellint_rd for details.
   */
  inline float
  ellint_rdf(float __x, float __y, float __z)
  { return std::__detail::__ellint_rd<float>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function @f$ R_D(x,y,z) @f$.
   *
   * @see ellint_rd for details.
   */
  inline long double
  ellint_rdl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rd<long double>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function of the second kind
   * @f$ R_D(x,y,z) = R_J(x,y,z,z) @f$ where
   * @f$ R_J(x,y,z,p) @f$ is the Carlson elliptic function
   * of the third kind.
   *
   * The Carlson elliptic function of the second kind is defined by:
   * @f[
   *    R_D(x,y,z) = \frac{3}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{3/2}}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *    by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first of two symmetric arguments.
   * @param  __y  The second of two symmetric arguments.
   * @param  __z  The third argument.
   */
  template<typename _Tp, typename _Up, typename _Vp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>
    ellint_rd(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rd<__type>(__x, __y, __z);
    }

  // Complete Carlson elliptic R_G functions

  /**
   * Return the Carlson complementary elliptic function @f$ R_G(x,y) @f$.
   *
   * @see comp_ellint_rg for details.
   */
  inline float
  comp_ellint_rg(float __x, float __y)
  { return std::__detail::__comp_ellint_rg<float>(__x, __y); }

  /**
   * Return the Carlson complementary elliptic function @f$ R_G(x,y) @f$.
   *
   * @see comp_ellint_rg for details.
   */
  inline long double
  comp_ellint_rg(long double __x, long double __y)
  { return std::__detail::__comp_ellint_rg<long double>(__x, __y); }

  /**
   * 
   */
  template<typename _Tx, typename _Ty>
    inline __gnu_cxx::__promote_num_t<_Tx, _Ty>
    comp_ellint_rg(_Tx __x, _Ty __y)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tx, _Ty>;
      return std::__detail::__comp_ellint_rg<__type>(__x, __y);
    }

  // Carlson elliptic R_G functions

  /**
   * Return the Carlson elliptic function @f$ R_G(x,y) @f$.
   *
   * @see ellint_rg for details.
   */
  inline float
  ellint_rgf(float __x, float __y, float __z)
  { return std::__detail::__ellint_rg<float>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function @f$ R_G(x,y) @f$.
   *
   * @see ellint_rg for details.
   */
  inline long double
  ellint_rgl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rg<long double>(__x, __y, __z); }

  /**
   * Return the symmetric Carlson elliptic function of the second kind
   * @f$ R_G(x,y,z) @f$.
   *
   * The Carlson symmetric elliptic function of the second kind is defined by:
   * @f[
   *    R_G(x,y,z) = \frac{1}{4} \int_0^\infty
   *            dt t [(t + x)(t + y)(t + z)]^{-1/2}
   *         (\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z})
   *  @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *    by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first of three symmetric arguments.
   * @param  __y  The second of three symmetric arguments.
   * @param  __z  The third of three symmetric arguments.
   */
  template<typename _Tp, typename _Up, typename _Vp>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>
    ellint_rg(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rg<__type>(__x, __y, __z);
    }

  // Hurwitz zeta functions

  inline float
  hurwitz_zetaf(float __s, float __a)
  { return std::__detail::__hurwitz_zeta<float>(__s, __a); }

  inline long double
  hurwitz_zetal(long double __s, long double __a)
  { return std::__detail::__hurwitz_zeta<long double>(__s, __a); }

  /**
   * Return the Hurwitz zeta function of argument @c s, and parameter @c a.
   *
   * The the Hurwitz zeta function is defined by
   * @f[
   *    \zeta(s, a) = 
   * @f]
   *
   * @param __s The argument
   * @param __a The parameter
   */
  template<typename _Tp, typename _Up>
    inline __gnu_cxx::__promote_num_t<_Tp, _Up>
    hurwitz_zeta(_Tp __s, _Up __a)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Up>;
      return std::__detail::__hurwitz_zeta<__type>(__s, __a);
    }

  // Digamma or psi functions

  inline float
  psif(float __x)
  { return std::__detail::__psi<float>(__x); }

  inline long double
  psil(long double __x)
  { return std::__detail::__psi<long double>(__x); }

  /**
   * Return the psi or digamma function of argument @c x.
   *
   * The the psi or digamma function is defined by
   * @f[
   *    \psi(x) = 
   * @f]
   *
   * @param __x The parameter
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    psi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__psi<__type>(__x);
    }

  // Incomplete beta functions

  /**
   * Return the regularized incomplete beta function of parameters @c a, @c b,
   * and argument @c x.
   *
   * See ibeta for details.
   */
  inline float
  ibetaf(float __a, float __b, float __x)
  { return std::__detail::__beta_inc<float>(__a, __b, __x); }

  /**
   * Return the regularized incomplete beta function of parameters @c a, @c b,
   * and argument @c x.
   *
   * See ibeta for details.
   */
  inline long double
  ibetal(long double __a, long double __b, long double __x)
  { return std::__detail::__beta_inc<long double>(__a, __b, __x); }

  /**
   * Return the regularized incomplete beta function of parameters @c a, @c b,
   * and argument @c x.
   *
   * The regularized incomplete beta function is defined by
   * @f[
   *    I_x(a, b) = \frac{B_x(a,b)}{B(a,b)}
   * @f]
   * where
   * @f[
   *   B_x(a,b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt
   * @f]
   * is the non-regularized beta function and @f$ B(a,b) @f$
   * is the usual beta function.
   * @f]
   *
   * @param __a The first parameter
   * @param __b The second parameter
   * @param __x The argument
   */
  template<typename _Ta, typename _Tb, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Ta, _Tb, _Tp>
    ibeta(_Ta __a, _Tb __b, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Ta, _Tb, _Tp>;
      return std::__detail::__beta_inc<__type>(__a, __b, __x);
    }

  // Complementary incomplete beta functions

  inline float
  ibetacf(float __a, float __b, float __x)
  { return 1.0F - ibetaf(__a, __b, __x); }

  inline long double
  ibetacl(long double __a, long double __b, long double __x)
  { return 1.0L - ibetal(__a, __b, __x); }

  /**
   * Return the regularized complementary incomplete beta function
   * of parameters @c a, @c b, and argument @c x.
   *
   * The regularized complementary incomplete beta function is defined by
   * @f[
   *    I_x(a, b) = I_x(a, b)
   * @f]
   *
   * @param __a The parameter
   * @param __b The parameter
   * @param __x The argument
   */
  template<typename _Ta, typename _Tb, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Ta, _Tb, _Tp>
    ibetac(_Ta __a, _Tb __b, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Ta, _Tb, _Tp>;
      return __type(1) - ibeta<__type>(__a, __b, __x);
    }

  // Fresnel sine integral

  inline float
  fresnel_sf(float __x)
  { return std::imag(std::__detail::__fresnel<float>(__x)); }

  inline long double
  fresnel_sl(long double __x)
  { return std::imag(std::__detail::__fresnel<long double>(__x)); }

  /**
   * Return the Fresnel sine integral of argument @c x.
   *
   * The Fresnel sine integral is defined by
   * @f[
   *    S(x) = \int_0^x \sin(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    fresnel_s(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::imag(std::__detail::__fresnel<__type>(__x));
    }

  // Fresnel cosine integral

  inline float
  fresnel_cf(float __x)
  { return std::real(std::__detail::__fresnel<float>(__x)); }

  inline long double
  fresnel_cl(long double __x)
  { return std::real(std::__detail::__fresnel<long double>(__x)); }

  /**
   * Return the Fresnel cosine integral of argument @c x.
   *
   * The Fresnel cosine integral is defined by
   * @f[
   *    C(x) = \int_0^x \cos(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    fresnel_c(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::real(std::__detail::__fresnel<__type>(__x));
    }

  // Dawson integral

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for @c float argument @c x.
   *
   * @see dawson for details.
   */
  inline float
  dawsonf(float __x)
  { return std::__detail::__dawson<float>(__x); }

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for @c long double argument @c x.
   *
   * @see dawson for details.
   */
  inline long double
  dawsonl(long double __x)
  { return std::__detail::__dawson<long double>(__x); }

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for real argument @c x.
   *
   * The Dawson integral is defined by:
   * @f[
   *    F(x) = e^{-x^2}\int_0^x e^{y^2}dy
   * @f]
   * and it's derivative is:
   * @f[
   *    F'(x) = 1 - 2xF(x)
   * @f]
   *
   * @param __x The argument @f$ -inf < x < inf @f$.
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    dawson(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__dawson<__type>(__x);
    }

  inline float
  expint_e1f(float __x)
  { return std::__detail::__expint_E1<float>(__x); }

  inline long double
  expint_e1l(long double __x)
  { return std::__detail::__expint_E1<long double>(__x); }

  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    expint_e1(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__expint_E1<__type>(__x);
    }

  inline float
  expint_enf(unsigned int __n, float __x)
  { return std::__detail::__expint<float>(__n, __x); }

  inline long double
  expint_enl(unsigned int __n, long double __x)
  { return std::__detail::__expint<long double>(__n, __x); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    expint_en(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__expint<__type>(__n, __x);
    }

  //  Log upper Pochhammer symbol

  inline float
  lpochhammer_uf(float __a, float __n)
  { return std::__detail::__log_pochhammer_u<float>(__a, __n); }

  inline long double
  lpochhammer_ul(long double __a, long double __n)
  { return std::__detail::__log_pochhammer_u<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tn>
    lpochhammer_u(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tn>;
      return std::__detail::__log_pochhammer_u<__type>(__a, __n);
    }

  //  Log lower Pochhammer symbol

  inline float
  lpochhammer_lf(float __a, float __n)
  { return std::__detail::__log_pochhammer_l<float>(__a, __n); }

  inline long double
  lpochhammer_ll(long double __a, long double __n)
  { return std::__detail::__log_pochhammer_l<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tn>
    lpochhammer_l(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tn>;
      return std::__detail::__log_pochhammer_l<__type>(__a, __n);
    }

  //  Upper Pochhammer symbols (see boost::rising_factorial)

  inline float
  pochhammer_uf(float __a, float __n)
  { return std::__detail::__pochhammer_u<float>(__a, __n); }

  inline long double
  pochhammer_ul(long double __a, long double __n)
  { return std::__detail::__pochhammer_u<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tn>
    pochhammer_u(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tn>;
      return std::__detail::__pochhammer_u<__type>(__a, __n);
    }

  //  Lower Pochhammer symbols (see boost::falling_factorial)

  inline float
  pochhammer_lf(float __a, float __n)
  { return std::__detail::__pochhammer_l<float>(__a, __n); }

  inline long double
  pochhammer_ll(long double __a, long double __n)
  { return std::__detail::__pochhammer_l<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tn>
    pochhammer_l(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tn>;
      return std::__detail::__pochhammer_l<__type>(__a, __n);
    }

  // Factorial

  inline float
  factorialf(unsigned int __n)
  { return std::__detail::__factorial<float>(__n); }

  inline long double
  factoriall(unsigned int __n)
  { return std::__detail::__factorial<long double>(__n); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    factorial(unsigned int __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__factorial<__type>(__n);
    }

  // Double factorial

  inline float
  double_factorialf(int __n)
  { return std::__detail::__double_factorial<float>(__n); }

  inline long double
  double_factoriall(int __n)
  { return std::__detail::__double_factorial<long double>(__n); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    double_factorial(int __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__double_factorial<__type>(__n);
    }

  // Log factorial

  inline float
  lfactorialf(unsigned int __n)
  { return std::__detail::__log_factorial<float>(__n); }

  inline long double
  lfactoriall(unsigned int __n)
  { return std::__detail::__log_factorial<long double>(__n); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    lfactorial(unsigned int __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__log_factorial<__type>(__n);
    }

  // Log double factorial

  inline float
  ldouble_factorialf(int __n)
  { return std::__detail::__log_double_factorial<float>(__n); }

  inline long double
  ldouble_factoriall(int __n)
  { return std::__detail::__log_double_factorial<long double>(__n); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    ldouble_factorial(int __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__log_double_factorial<__type>(__n);
    }

  // Binomial coefficient

  inline float
  bincoeff(unsigned int __n, unsigned int __k)
  { return std::__detail::__bincoef<float>(__n, __k); }

  inline long double
  bincoefl(unsigned int __n, unsigned int __k)
  { return std::__detail::__bincoef<long double>(__n, __k); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    bincoef(unsigned int __n, unsigned int __k)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__bincoef<__type>(__n, __k);
    }

  // Log binomial coefficient

  inline float
  lbincoeff(unsigned int __n, unsigned int __k)
  { return std::__detail::__log_bincoef<float>(__n, __k); }

  inline long double
  lbincoefl(unsigned int __n, unsigned int __k)
  { return std::__detail::__log_bincoef<long double>(__n, __k); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    lbincoef(unsigned int __n, unsigned int __k)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__log_bincoef<__type>(__n, __k);
    }

  // Bernoulli numbers

  inline float
  bernoullif(unsigned int __n)
  { return std::__detail::__bernoulli<float>(__n); }

  inline long double
  bernoullil(unsigned int __n)
  { return std::__detail::__bernoulli<long double>(__n); }

  /**
   * Return the Bernoulli number of integer order @c n.
   *
   * The Bernoulli numbers are defined by
   * @f[
   *    
   * @f]
   *
   * @param __n The order.
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    bernoulli(unsigned int __n)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__bernoulli<__type>(__n);
    }

  // Legendre functions of the second kind

  inline float
  legendre_qf(unsigned int __n, float __x)
  { return std::__detail::__poly_legendre_q<float>(__n, __x); }

  inline long double
  legendre_ql(unsigned int __n, long double __x)
  { return std::__detail::__poly_legendre_q<long double>(__n, __x); }

  /**
   * 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    legendre_q(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__poly_legendre_q<__type>(__n, __x);
    }

  // Scaled lower incomplete gamma

  inline float
  gamma_pf(float __a, float __x)
  { return std::__detail::__gamma_p<float>(__a, __x); }

  inline long double
  gamma_pl(long double __a, long double __x)
  { return std::__detail::__gamma_p<long double>(__a, __x); }

  /**
   * 
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Ta, _Tp>
    gamma_p(_Ta __a, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Ta, _Tp>;
      return std::__detail::__gamma_p<__type>(__a, __x);
    }

  // Scaled upper incomplete gamma

  inline float
  gamma_qf(float __a, float __x)
  { return std::__detail::__gamma_q<float>(__a, __x); }

  inline long double
  gamma_ql(long double __a, long double __x)
  { return std::__detail::__gamma_q<long double>(__a, __x); }

  /**
   * 
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Ta, _Tp>
    gamma_q(_Ta __a,_Tp  __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Ta, _Tp>;
      return std::__detail::__gamma_q<__type>(__a, __x);
    }

  // Jacobi zeta functions.

  inline float
  jacobi_zetaf(float __k, float __phi)
  { return std::__detail::__jacobi_zeta<float>(__k, __phi); }

  inline long double
  jacobi_zetal(long double __k, long double __phi)
  { return std::__detail::__jacobi_zeta<long double>(__k, __phi); }

  /**
   * Return the Jacobi zeta function of @c k and @f$ @c \phi @f$.
   *
   * The Jacobi zeta function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus
   * @param __phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline __gnu_cxx::__promote_num_t<_Tk, _Tphi>
    jacobi_zeta(_Tk __k, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tk, _Tphi>;
      return std::__detail::__jacobi_zeta<__type>(__k, __phi);
    }

  // Heuman lambda functions.

  inline float
  heuman_lambdaf(float __k, float __phi)
  { return std::__detail::__heuman_lambda<float>(__k, __phi); }

  inline long double
  heuman_lambdal(long double __k, long double __phi)
  { return std::__detail::__heuman_lambda<long double>(__k, __phi); }

  /**
   * Return the Heuman lambda function of @c k and @f$ @c \phi @f$.
   *
   * The complete Heuman lambda function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus
   * @param __phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline __gnu_cxx::__promote_num_t<_Tk, _Tphi>
    heuman_lambda(_Tk __k, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tk, _Tphi>;
      return std::__detail::__heuman_lambda<__type>(__k, __phi);
    }

  // Complete Legendre elliptic integral D.

  inline float
  comp_ellint_df(float __k)
  { return std::__detail::__comp_ellint_d<float>(__k); }

  inline long double
  comp_ellint_dl(long double __k)
  { return std::__detail::__comp_ellint_d<long double>(__k); }

  /**
   * Return the complete Legendre elliptic integral D of @c k and @c @f$ \phi @f$.
   *
   * The complete Legendre elliptic integral D is defined by
   * @f[
   *    D(k) = \int_0^{\pi/2} \frac{\sin^2\theta d\theta}{\sqrt{1-k^2sin2\theta}}
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   */
  template<typename _Tk>
    inline __gnu_cxx::__promote_num_t<_Tk>
    comp_ellint_d(_Tk __k)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tk>;
      return std::__detail::__comp_ellint_d<__type>(__k);
    }

  // Legendre elliptic integrals D.

  inline float
  ellint_df(float __k, float __phi)
  { return std::__detail::__ellint_d<float>(__k, __phi); }

  inline long double
  ellint_dl(long double __k, long double __phi)
  { return std::__detail::__ellint_d<long double>(__k, __phi); }

  /**
   * Return the incomplete Legendre elliptic integral D of @c k and @c @f$ \phi @f$.
   *
   * The Legendre elliptic integral D is defined by
   * @f[
   *    D(k,\phi) = \int_0^\phi \frac{\sin^2\theta d\theta}{\sqrt{1-k^2sin2\theta}}
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline __gnu_cxx::__promote_num_t<_Tk, _Tphi>
    ellint_d(_Tk __k, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tk, _Tphi>;
      return std::__detail::__ellint_d<__type>(__k, __phi);
    }

  // Bulirsch elliptic integrals of the first kind.

  inline float
  ellint_el1f(float __x, float __k_c)
  { return std::__detail::__ellint_el1<float>(__x, __k_c); }

  inline long double
  ellint_el1l(long double __x, long double __k_c)
  { return std::__detail::__ellint_el1<long double>(__x, __k_c); }

  /**
   * Return the Bulirsch elliptic integral of the first kind of ...
   *
   * The Bulirsch elliptic integral of the first kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The argument
   * @param __k_c The complementary modulus @f$ k_c = \sqrt(1 - k^2) @f$
   */
  template<typename _Tp, typename _Tk>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tk>
    ellint_el1(_Tp __x, _Tk __k_c)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tk>;
      return std::__detail::__ellint_el1<__type>(__x, __k_c);
    }

  // Bulirsch elliptic integrals of the second kind.

  inline float
  ellint_el2f(float __x, float __k_c, float __a, float __b)
  { return std::__detail::__ellint_el2<float>(__x, __k_c, __a, __b); }

  inline long double
  ellint_el2l(long double __x, long double __k_c,
	      long double __a, long double __b)
  { return std::__detail::__ellint_el2<long double>(__x, __k_c, __a, __b); }

  /**
   * Return the Bulirsch elliptic integral of the second kind of ...
   *
   * The Bulirsch elliptic integral of the second kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The argument
   * @param __k_c The complementary modulus @f$ k_c = \sqrt(1 - k^2) @f$
   * @param __a The 
   * @param __b The 
   */
  template<typename _Tp, typename _Tk, typename _Ta, typename _Tb>
    inline __gnu_cxx::__promote_num_t<_Tp, _Tk, _Ta, _Tb>
    ellint_el2(_Tp __x, _Tk __k_c, _Ta __a, _Tb __b)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp, _Tk, _Ta, _Tb>;
      return std::__detail::__ellint_el2<__type>(__x, __k_c, __a, __b);
    }

  // Bulirsch elliptic integrals of the third kind.

  inline float
  ellint_el3f(float __x, float __k_c, float __p)
  { return std::__detail::__ellint_el3<float>(__x, __k_c, __p); }

  inline long double
  ellint_el3l(long double __x, long double __k_c, long double __p)
  { return std::__detail::__ellint_el3<long double>(__x, __k_c, __p); }

  /**
   * Return the Bulirsch elliptic integral of the third kind of ...
   *
   * The Bulirsch elliptic integral of the third kind is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x The 
   * @param __k_c The complementary modulus @f$ k_c = \sqrt(1 - k^2) @f$
   * @param __p The 
   */
  template<typename _Tx, typename _Tk, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tx, _Tk, _Tp>
    ellint_el3(_Tx __x, _Tk __k_c, _Tp __p)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tx, _Tk, _Tp>;
      return std::__detail::__ellint_el3<__type>(__x, __k_c, __p);
    }

  // Bulirsch complete elliptic integrals.

  inline float
  ellint_celf(float __k_c, float __p, float __a, float __b)
  { return std::__detail::__ellint_cel<float>(__k_c, __p, __a, __b); }

  inline long double
  ellint_cell(long double __k_c, long double __p,
	      long double __a, long double __b)
  { return std::__detail::__ellint_cel<long double>(__k_c, __p, __a, __b); }

  /**
   * Return the Bulirsch complete elliptic integral of ...
   *
   * The Bulirsch complete elliptic integral is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k_c The complementary modulus @f$ k_c = \sqrt(1 - k^2) @f$
   * @param __p The 
   * @param __a The 
   * @param __b The 
   */
  template<typename _Tk, typename _Tp, typename _Ta, typename _Tb>
    inline __gnu_cxx::__promote_num_t<_Tk, _Tp, _Ta, _Tb>
    ellint_cel(_Tk __k_c, _Tp __p, _Ta __a, _Tb __b)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tk, _Tp, _Ta, _Tb>;
      return std::__detail::__ellint_cel<__type>(__k_c, __p, __a, __b);
    }

  // Cylindrical Hankel functions of the first kind.

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * of complex order @c @f$ \nu @f$ and complex argument @c x.
   *
   * @see cyl_hankel_1 for more details.
   */
  inline std::complex<float>
  cyl_hankel_1f(std::complex<float> __nu, std::complex<float> __x)
  { return std::__detail::__cyl_hankel_1<float>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * of complex order @c @f$ \nu @f$ and complex argument @c x.
   *
   * @see cyl_hankel_1 for more details.
   */
  inline std::complex<long double>
  cyl_hankel_1l(std::complex<long double> __nu, std::complex<long double> __x)
  { return std::__detail::__cyl_hankel_1<long double>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * of complex order @c @f$ \nu @f$ and complex argument @c x.
   *
   * The cylindrical Hankel function of the first kind is defined by
   * @f[
   *    H^{(1)}_\nu(x) = J_\nu(x) + i N_\nu(x)
   * @f]
   *
   * @param __nu The complex order
   * @param __x The complex argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tpnu, _Tp>>
    cyl_hankel_1(std::complex<_Tpnu> __nu, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_1<__type>(__nu, __x);
    }

  // Cylindrical Hankel functions of the second kind.

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * of complex order @c @f$ \nu @f$ and complex argument @c x.
   *
   * @see cyl_hankel_2 for more details.
   */
  inline std::complex<float>
  cyl_hankel_2f(std::complex<float> __nu, std::complex<float> __x)
  { return std::__detail::__cyl_hankel_2<float>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * of complex order @c @f$ \nu @f$ and complex argument @c x.
   *
   * @see cyl_hankel_2 for more details.
   */
  inline std::complex<long double>
  cyl_hankel_2l(std::complex<long double> __nu, std::complex<long double> __x)
  { return std::__detail::__cyl_hankel_2<long double>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * of complex order @c @f$ \nu @f$ and complex argument @c x.
   *
   * The cylindrical Hankel function of the second kind is defined by
   * @f[
   *    H^{(2)}_\nu(x) = J_\nu(x) - i N_\nu(x)
   * @f]
   *
   * @param __nu The complex order
   * @param __x The complex argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tpnu, _Tp>>
    cyl_hankel_2(std::complex<_Tpnu> __nu, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_2<__type>(__nu, __x);
    }

  // Spherical Hankel functions of the first kind.

  inline std::complex<float>
  sph_hankel_1f(unsigned int __n, std::complex<float> __x)
  { return std::__detail::__sph_hankel_1<float>(__n, __x); }

  inline std::complex<long double>
  sph_hankel_1l(unsigned int __n, std::complex<long double> __x)
  { return std::__detail::__sph_hankel_1<long double>(__n, __x); }

  /**
   * Return the complex spherical Hankel function of the first kind
   * of non-negative order @c n and complex argument @c x.
   *
   * The spherical Hankel function of the first kind is defined by
   * @f[
   *    h^{(1)}_n(x) = j_n(x) + i n_n(x)
   * @f]
   *
   * @param __n The integral order >= 0
   * @param __x The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    sph_hankel_1(unsigned int __n, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sph_hankel_1<__type>(__n, __x);
    }

  // Spherical Hankel functions of the second kind.

  inline std::complex<float>
  sph_hankel_2f(unsigned int __n, std::complex<float> __x)
  { return std::__detail::__sph_hankel_2<float>(__n, __x); }

  inline std::complex<long double>
  sph_hankel_2l(unsigned int __n, std::complex<long double> __x)
  { return std::__detail::__sph_hankel_2<long double>(__n, __x); }

  /**
   * Return the complex spherical Hankel function of the second kind
   * of non-negative order @c n and complex argument @c x.
   *
   * The spherical Hankel function of the second kind is defined by
   * @f[
   *    h^{(2)}_n(x) = j_n(x) - i n_n(x)
   * @f]
   *
   * @param __n The integral order >= 0
   * @param __x The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    sph_hankel_2(unsigned int __n, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__sph_hankel_2<__type>(__n, __x);
    }

  // Spherical harmonic functions

  /**
   * Return the complex spherical harmonic function of degree @c l, order @c m,
   * and real zenith angle @f$ \theta @f$, and real azimuth angle @c @f$ \phi @f$.
   *
   * @see sph_harmonic for details.
   */
  inline std::complex<float>
  sph_harmonicf(unsigned int __l, int __m,
		float __theta, float __phi)
  { return std::__detail::__sph_harmonic<float>(__l, __m, __theta, __phi); }

  /**
   * Return the complex spherical harmonic function of degree @c l, order @c m,
   * and real zenith angle @f$ \theta @f$, and real azimuth angle @c @f$ \phi @f$.
   *
   * @see sph_harmonic for details.
   */
  inline std::complex<long double>
  sph_harmonicl(unsigned int __l, int __m,
		long double __theta, long double __phi)
  {
    return std::__detail::__sph_harmonic<long double>(__l, __m, __theta, __phi);
  }

  /**
   * Return the complex spherical harmonic function of degree @c l, order @c m,
   * and real zenith angle @f$ \theta @f$, and real azimuth angle @c @f$ \phi @f$.
   *
   * The spherical harmonic function is defined by:
   * @f[
   *    Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   *                                \frac{(l-m)!}{(l+m)!}]
   *                     P_l^{|m|}(\cos\theta) \exp^{im\phi}
   * @f]
   *
   * @param __l The order
   * @param __m The degree
   * @param __theta The zenith angle in radians
   * @param __phi The azimuth angle in radians
   */
  template<typename _Ttheta, typename _Tphi>
    inline std::complex<__gnu_cxx::__promote_num_t<_Ttheta, _Tphi>>
    sph_harmonic(unsigned int __l, int __m, _Ttheta __theta, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_num_t<_Ttheta, _Tphi>;
      return std::__detail::__sph_harmonic<__type>(__l, __m, __theta, __phi);
    }

  // Polylogarithm functions

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @c w.
   *
   * @see polylog for details.
   */
  inline std::complex<float>
  polylogf(float __s, std::complex<float> __w)
  { return std::__detail::__polylog<float>(__s, __w); }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @c w.
   *
   * @see polylog for details.
   */
  inline std::complex<long double>
  polylogl(long double __s, std::complex<long double> __w)
  { return std::__detail::__polylog<long double>(__s, __w); }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @c w.
   *
   * The polylogarithm function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __s 
   * @param __w 
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    polylog(_Tp __s, std::complex<_Tp> __w)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__polylog<__type>(__s, __w);
    }

  // Dirichlet eta function

  /**
   * Return the Dirichlet eta function of real argument @c x.
   *
   * @see dirichlet_eta for details.
   */
  inline float
  dirichlet_etaf(float __x)
  { return std::__detail::__dirichlet_eta<float>(__x); }

  /**
   * Return the Dirichlet eta function of real argument @c x.
   *
   * @see dirichlet_eta for details.
   */
  inline long double
  dirichlet_etal(long double __x)
  { return std::__detail::__dirichlet_eta<long double>(__x); }

  /**
   * Return the Dirichlet eta function of real argument @c x.
   *
   * The Dirichlet eta function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x 
   */
  template<typename _Tp>
    inline _Tp
    dirichlet_eta(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__dirichlet_eta<__type>(__x);
    }

  // Dirichlet beta function

  inline float
  dirichlet_betaf(float __x)
  { return std::__detail::__dirichlet_beta<float>(__x); }

  inline long double
  dirichlet_betal(long double __x)
  { return std::__detail::__dirichlet_beta<long double>(__x); }

  /**
   * Return the Dirichlet beta function of real argument @c x.
   *
   * The Dirichlet beta function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __x 
   */
  template<typename _Tp>
    inline _Tp
    dirichlet_beta(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__dirichlet_beta<__type>(__x);
    }

  // Clausen S functions

  /**
   * Return the Clausen sine function of order @c m and real argument @c x.
   *
   * @see clausen_s for details.
   */
  inline float
  clausen_sf(unsigned int __m, float __w)
  { return std::__detail::__clausen_s<float>(__m, __w); }

  /**
   * Return the Clausen sine function of order @c m and real argument @c x.
   *
   * @see clausen_s for details.
   */
  inline long double
  clausen_sl(unsigned int __m, long double __w)
  { return std::__detail::__clausen_s<long double>(__m, __w); }

  /**
   * Return the Clausen sine function of order @c m and real argument @c x.
   *
   * The Clausen sine function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __m 
   * @param __w 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    clausen_s(unsigned int __m, _Tp __w)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__clausen_s<__type>(__m, __w);
    }

  // Clausen C functions

  /**
   * Return the Clausen cosine function of order @c m and real argument @c x.
   *
   * @see clausen_c for details.
   */
  inline float
  clausen_cf(unsigned int __m, float __w)
  { return std::__detail::__clausen_c<float>(__m, __w); }

  /**
   * Return the Clausen cosine function of order @c m and real argument @c x.
   *
   * @see clausen_c for details.
   */
  inline long double
  clausen_cl(unsigned int __m, long double __w)
  { return std::__detail::__clausen_c<long double>(__m, __w); }

  /**
   * Return the Clausen cosine function of order @c m and real argument @c x.
   *
   * The Clausen cosine function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __m 
   * @param __w 
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    clausen_c(unsigned int __m, _Tp __w)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__clausen_c<__type>(__m, __w);
    }

  // Clausen functions - real argument

  /**
   * Return the Clausen function of integer order @c m and complex argument @c w.
   *
   * @see clausen for details.
   */
  inline float
  clausenf(unsigned int __m, float __w)
  { return std::__detail::__clausen<float>(__m, __w); }

  /**
   * Return the Clausen function of integer order @c m and complex argument @c w.
   *
   * @see clausen for details.
   */
  inline long double
  clausenl(unsigned int __m, long double __w)
  { return std::__detail::__clausen<long double>(__m, __w); }

  /**
   * Return the Clausen function of integer order @c m and complex argument @c w.
   *
   * The Clausen function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __m 
   * @param __w The complex argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tp>
    clausen(unsigned int __m, _Tp __w)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__clausen<__type>(__m, __w);
    }

  // Clausen functions - complex argument

  inline std::complex<float>
  clausenf(unsigned int __m, std::complex<float> __w)
  { return std::__detail::__clausen<float>(__m, __w); }

  inline std::complex<long double>
  clausenl(unsigned int __m, std::complex<long double> __w)
  { return std::__detail::__clausen<long double>(__m, __w); }

  /**
   * 
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_num_t<_Tp>>
    clausen(unsigned int __m, std::complex<_Tp> __w)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__clausen<__type>(__m, __w);
    }

  // Exponential theta_1 functions.

  /**
   * Return the exponential theta-1 function of period @c nu and argument @c x.
   *
   * @see theta_1 for details.
   */
  inline float
  theta_1f(float __nu, float __x)
  { return std::__detail::__theta_1<float>(__nu, __x); }

  /**
   * Return the exponential theta-1 function of period @c nu and argument @c x.
   *
   * @see theta_1 for details.
   */
  inline long double
  theta_1l(long double __nu, long double __x)
  { return std::__detail::__theta_1<long double>(__nu, __x); }

  /**
   * Return the exponential theta-1 function of period @c nu and argument @c x.
   *
   * The Neville theta-1 function is defined by
   * @f[
   *    \theta_1(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j - 1/2)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    theta_1(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__theta_1<__type>(__nu, __x);
    }

  // Exponential theta_2 functions.

  /**
   * Return the exponential theta-2 function of period @c nu and argument @c x.
   *
   * @see theta_2 for details.
   */
  inline float
  theta_2f(float __nu, float __x)
  { return std::__detail::__theta_2<float>(__nu, __x); }

  /**
   * Return the exponential theta-2 function of period @c nu and argument @c x.
   *
   * @see theta_2 for details.
   */
  inline long double
  theta_2l(long double __nu, long double __x)
  { return std::__detail::__theta_2<long double>(__nu, __x); }

  /**
   * Return the exponential theta-2 function of period @c nu and argument @c x.
   *
   * The exponential theta-2 function is defined by
   * @f[
   *    \theta_2(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 2) argument
   * @param __x The argument
   */
  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    theta_2(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__theta_2<__type>(__nu, __x);
    }

  // Exponential theta_3 functions.

  /**
   * Return the exponential theta-3 function of period @c nu and argument @c x.
   *
   * @see theta_3 for details.
   */
  inline float
  theta_3f(float __nu, float __x)
  { return std::__detail::__theta_3<float>(__nu, __x); }

  /**
   * Return the exponential theta-3 function of period @c nu and argument @c x.
   *
   * @see theta_3 for details.
   */
  inline long double
  theta_3l(long double __nu, long double __x)
  { return std::__detail::__theta_3<long double>(__nu, __x); }

  /**
   * Return the exponential theta-3 function of period @c nu and argument @c x.
   *
   * The exponential theta-3 function is defined by
   * @f[
   *    \theta_3(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    \exp\left( \frac{-(\nu+j)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 1) argument
   * @param __x The argument
   */
  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    theta_3(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__theta_3<__type>(__nu, __x);
    }

  // Exponential theta_4 functions.

  /**
   * Return the exponential theta-4 function of period @c nu and argument @c x.
   *
   * @see theta_4 for details.
   */
  inline float
  theta_4f(float __nu, float __x)
  { return std::__detail::__theta_4<float>(__nu, __x); }

  /**
   * Return the exponential theta-4 function of period @c nu and argument @c x.
   *
   * @see theta_4 for details.
   */
  inline long double
  theta_4l(long double __nu, long double __x)
  { return std::__detail::__theta_4<long double>(__nu, __x); }

  /**
   * Return the exponential theta-4 function of period @c nu and argument @c x.
   *
   * The exponential theta-4 function is defined by
   * @f[
   *    \theta_4(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *                 \exp\left( \frac{-(\nu + j + 1/2)^2}{x} \right)
   * @f]
   *
   * @param __nu The periodic (period = 1) argument
   * @param __x The argument
   */
  template<typename _Tpnu, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpnu, _Tp>
    theta_4(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpnu, _Tp>;
      return std::__detail::__theta_4<__type>(__nu, __x);
    }

  // Elliptic nome function

  inline float
  ellnomef(float __k)
  { return std::__detail::__ellnome<float>(__k); }

  inline long double
  ellnomel(long double __k)
  { return std::__detail::__ellnome<long double>(__k); }

  /**
   * Return the elliptic nome function of modulus @c k.
   *
   * The elliptic nome function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   */
  template<typename _Tp>
    inline _Tp
    ellnome(_Tp __k)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tp>;
      return std::__detail::__ellnome<__type>(__k);
    }

  // Neville theta_s functions.

  /**
   * Return the Neville theta-s function of modulus @c k and argument @c x.
   *
   * @see theta_s for details.
   */
  inline float
  theta_sf(float __k, float __x)
  { return std::__detail::__theta_s<float>(__k, __x); }

  /**
   * Return the Neville theta-s function of modulus @c k and argument @c x.
   *
   * @see theta_s for details.
   */
  inline long double
  theta_sl(long double __k, long double __x)
  { return std::__detail::__theta_s<long double>(__k, __x); }

  /**
   * Return the Neville theta-s function of modulus @c k and argument @c x.
   *
   * The Neville theta-s function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __x The argument
   */
  template<typename _Tpk, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpk, _Tp>
    theta_s(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpk, _Tp>;
      return std::__detail::__theta_s<__type>(__k, __x);
    }

  // Neville theta_c functions.

  /**
   * Return the Neville theta-c function of modulus @c k and argument @c x.
   *
   * @see theta_c for details.
   */
  inline float
  theta_cf(float __k, float __x)
  { return std::__detail::__theta_c<float>(__k, __x); }

  /**
   * Return the Neville theta-c function of modulus @c k and argument @c x.
   *
   * @see theta_c for details.
   */
  inline long double
  theta_cl(long double __k, long double __x)
  { return std::__detail::__theta_c<long double>(__k, __x); }

  /**
   * Return the Neville theta-c function of modulus @c k and argument @c x.
   *
   * The Neville theta-c function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __x The argument
   */
  template<typename _Tpk, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpk, _Tp>
    theta_c(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpk, _Tp>;
      return std::__detail::__theta_c<__type>(__k, __x);
    }

  // Neville theta_d functions.

  /**
   * Return the Neville theta-d function of modulus @c k and argument @c x.
   *
   * @see theta_d for details.
   */
  inline float
  theta_df(float __k, float __x)
  { return std::__detail::__theta_d<float>(__k, __x); }

  /**
   * Return the Neville theta-d function of modulus @c k and argument @c x.
   *
   * @see theta_d for details.
   */
  inline long double
  theta_dl(long double __k, long double __x)
  { return std::__detail::__theta_d<long double>(__k, __x); }

  /**
   * Return the Neville theta-d function of modulus @c k and argument @c x.
   *
   * The Neville theta-d function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __x The argument
   */
  template<typename _Tpk, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpk, _Tp>
    theta_d(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpk, _Tp>;
      return std::__detail::__theta_d<__type>(__k, __x);
    }

  // Neville theta_n functions.

  /**
   * Return the Neville theta-n function of modulus @c k and argument @c x.
   *
   * @see theta_n for details.
   */
  inline float
  theta_nf(float __k, float __x)
  { return std::__detail::__theta_n<float>(__k, __x); }

  /**
   * Return the Neville theta-n function of modulus @c k and argument @c x.
   *
   * @see theta_n for details.
   */
  inline long double
  theta_nl(long double __k, long double __x)
  { return std::__detail::__theta_n<long double>(__k, __x); }

  /**
   * Return the Neville theta-n function of modulus @c k and argument @c x.
   *
   * The Neville theta-n function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __x The argument
   */
  template<typename _Tpk, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tpk, _Tp>
    theta_n(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tpk, _Tp>;
      return std::__detail::__theta_n<__type>(__k, __x);
    }

  // Owens T functions.

  /**
   * Return the Owens T function function of thing @c h and argument @c a.
   *
   * @see owens_t for details.
   */
  inline float
  owens_tf(float __h, float __a)
  { return std::__detail::__owens_t<float>(__h, __a); }

  /**
   * Return the Owens T function function of thing @c h and argument @c a.
   *
   * @see owens_t for details.
   */
  inline long double
  owens_tl(long double __h, long double __a)
  { return std::__detail::__owens_t<long double>(__h, __a); }

  /**
   * Return the Owens T function of thing1 @c h and thing2 @c a.
   *
   * The Owens T function is defined by
   * @f[
   *    T(h,a) = \frac{1}{2\pi}\int_0^a
   *           \frac{\exp\left[-\frac{1}{2}h^2(1+x^2)\right]}{1+x^2} dx
   * @f]
   *
   * @param __h The shape factor
   * @param __a The integration lomit
   */
  template<typename _Tph, typename _Tpa>
    inline __gnu_cxx::__promote_num_t<_Tph, _Tpa>
    owens_t(_Tph __h, _Tpa __a)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tph, _Tpa>;
      return std::__detail::__owens_t<__type>(__h, __a);
    }
/*
  // Fermi-Dirac integrals.

  inline float
  fermi_diracf(float __s, float __x)
  { return std::__detail::__fermi_dirac<float>(__s, __x); }

  inline long double
  fermi_diracl(long double __s, long double __x)
  { return std::__detail::__fermi_dirac<long double>(__s, __x); }

  template<typename _Tps, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tps, _Tp>
    fermi_dirac(_Tps __s, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tps, _Tp>;
      return std::__detail::__fermi_dirac<__type>(__s, __x);
    }

  // Bose-Einstein integrals.

  inline float
  bose_einsteinf(float __s, float __x)
  { return std::__detail::__bose_einstein<float>(__s, __x); }

  inline long double
  bose_einsteinl(long double __s, long double __x)
  { return std::__detail::__bose_einstein<long double>(__s, __x); }

  template<typename _Tps, typename _Tp>
    inline __gnu_cxx::__promote_num_t<_Tps, _Tp>
    bose_einstein(_Tps __s, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_num_t<_Tps, _Tp>;
      return std::__detail::__bose_einstein<__type>(__s, __x);
    }
*/
#endif // __cplusplus >= 201103L

  /* @} */ // gnu_math_spec_func

} // namespace __gnu_cxx

#pragma GCC visibility pop

#endif // _GLIBCXX_BITS_SPECFUN_H
