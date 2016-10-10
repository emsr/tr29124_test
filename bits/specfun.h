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
 * Do not attempt to use it directly. @headername{cmath}
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
#  include <bits/numeric_limits.h>
#  include <bits/complex_util.h>
#  include <bits/sf_trig.tcc>
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
#  include <bits/sf_hermite.tcc>
#  include <bits/sf_theta.tcc>
#  include <bits/sf_trigint.tcc>
#  include <bits/sf_zeta.tcc>
#  include <bits/sf_owens_t.tcc>
#  include <bits/sf_polylog.tcc>
#  include <bits/sf_airy.tcc>
#  include <bits/sf_hankel.tcc>
#  include <bits/sf_distributions.tcc>
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
   * @defgroup math_spec_func C++ Mathematical Special Functions
   * @ingroup numerics
   *
   * A collection of advanced mathematical special functions.
   * @{
   */

  /**
   * @mainpage Mathematical Special Functions
   *
   * @section intro Introduction and History
   * The first significant library upgrade on the road to C++2011,
   * <a href="http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1836.pdf">
   * TR1</a>, included a set of 23 mathematical functions that significantly
   * extended the standard transcendental functions inherited from C and declared
   * in @<cmath@>.
   *
   * Although most components from TR1 were eventually adopted for C++11 these
   * math functions were left behind out of concern for implementability.
   * The math functions were published as a separate international standard
   * <a href="http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2010/n3060.pdf">
   * IS 29124 - Extensions to the C++ Library to Support Mathematical Special
   * Functions</a>.
   *
   * For C++17 these functions were incorporated into the main standard.
   *
   * @section contents Contents
   * The following functions are implemented in namespace @c std:
   * - @ref std::assoc_laguerre "assoc_laguerre - Associated Laguerre functions"
   * - @ref std::assoc_legendre "assoc_legendre - Associated Legendre functions"
   * - @ref std::beta "beta - Beta functions"
   * - @ref std::comp_ellint_1 "comp_ellint_1 - Complete elliptic functions of the first kind"
   * - @ref std::comp_ellint_2 "comp_ellint_2 - Complete elliptic functions of the second kind"
   * - @ref std::comp_ellint_3 "comp_ellint_3 - Complete elliptic functions of the third kind"
   * - @ref std::cyl_bessel_i "cyl_bessel_i - Regular modified cylindrical Bessel functions"
   * - @ref std::cyl_bessel_j "cyl_bessel_j - Cylindrical Bessel functions of the first kind"
   * - @ref std::cyl_bessel_k "cyl_bessel_k - Irregular modified cylindrical Bessel functions"
   * - @ref std::cyl_neumann "cyl_neumann - Cylindrical Neumann functions or Cylindrical Bessel functions of the second kind"
   * - @ref std::ellint_1 "ellint_1 - Incomplete elliptic functions of the first kind"
   * - @ref std::ellint_2 "ellint_2 - Incomplete elliptic functions of the second kind"
   * - @ref std::ellint_3 "ellint_3 - Incomplete elliptic functions of the third kind"
   * - @ref std::expint "expint - The exponential integral"
   * - @ref std::hermite "hermite - Hermite polynomials"
   * - @ref std::laguerre "laguerre - Laguerre functions"
   * - @ref std::legendre "legendre - Legendre polynomials"
   * - @ref std::riemann_zeta "riemann_zeta - The Riemann zeta function"
   * - @ref std::sph_bessel "sph_bessel - Spherical Bessel functions"
   * - @ref std::sph_legendre "sph_legendre - Spherical Legendre functions"
   * - @ref std::sph_neumann "sph_neumann - Spherical Neumann functions"
   *
   * The hypergeometric functions were stricken from the TR29124 and C++17
   * versions of this math library because of implementation concerns.
   * However, since they were in the TR1 version and since they are popular
   * we kept them as an extension in namespace @c __gnu_cxx:
   * - @ref __gnu_cxx::conf_hyperg "conf_hyperg - Confluent hypergeometric functions"
   * - @ref __gnu_cxx::hyperg "hyperg - Hypergeometric functions"
   *
   * In addition a large number of new functions are added as extensions:
   * - @ref __gnu_cxx::airy_ai "airy_ai - Airy functions of the first kind"
   * - @ref __gnu_cxx::airy_bi "airy_bi - Airy functions of the second kind"
   * - @ref __gnu_cxx::bincoef "bincoef - Binomial coefficients"
   * - @ref __gnu_cxx::bose_einstein "bose_einstein - Bose-Einstein integrals"
   * - @ref __gnu_cxx::chebyshev_t "chebyshev_t - Chebyshev polynomials of the first kind"
   * - @ref __gnu_cxx::chebyshev_u "chebyshev_u - Chebyshev polynomials of the second kind"
   * - @ref __gnu_cxx::chebyshev_v "chebyshev_v - Chebyshev polynomials of the third kind"
   * - @ref __gnu_cxx::chebyshev_w "chebyshev_w - Chebyshev polynomials of the fourth kind"
   * - @ref __gnu_cxx::clausen "clausen - Clausen integrals"
   * - @ref __gnu_cxx::clausen_c "clausen_c - Clausen cosine integrals"
   * - @ref __gnu_cxx::clausen_s "clausen_s - Clausen sine integrals"
   * - @ref __gnu_cxx::comp_ellint_d "comp_ellint_d - Incomplete Legendre D elliptic integral"
   * - @ref __gnu_cxx::conf_hyperg_lim "conf_hyperg_lim - Confluent hypergeometric limit functions"
   * - @ref __gnu_cxx::cos_pi "cos_pi - Reperiodized cosine function."
   * - @ref __gnu_cxx::cosh_pi "cosh_pi - Reperiodized hyperbolic cosine function."
   * - @ref __gnu_cxx::coshint "coshint - Hyperbolic cosine integral"
   * - @ref __gnu_cxx::cosint "cosint - Cosine integral"
   * - @ref __gnu_cxx::cyl_hankel_1 "cyl_hankel_1 - Cylindrical Hankel functions of the first kind"
   * - @ref __gnu_cxx::cyl_hankel_2 "cyl_hankel_2 - Cylindrical Hankel functions of the second kind"
   * - @ref __gnu_cxx::dawson "dawson - Dawson integrals"
   * - @ref __gnu_cxx::dilog "dilog - Dilogarithm functions"
   * - @ref __gnu_cxx::dirichlet_beta "dirichlet_beta - Dirichlet beta function"
   * - @ref __gnu_cxx::dirichlet_eta "dirichlet_eta - Dirichlet beta function"
   * - @ref __gnu_cxx::dirichlet_lambda "dirichlet_lambda - Dirichlet lambda function"
   * - @ref __gnu_cxx::double_factorial "double_factorial - "
   * - @ref __gnu_cxx::ellint_d "ellint_d - Legendre D elliptic integrals"
   * - @ref __gnu_cxx::ellint_rc "ellint_rc - Carlson elliptic functions R_C"
   * - @ref __gnu_cxx::ellint_rd "ellint_rd - Carlson elliptic functions R_D"
   * - @ref __gnu_cxx::ellint_rf "ellint_rf - Carlson elliptic functions R_F"
   * - @ref __gnu_cxx::ellint_rg "ellint_rg - Carlson elliptic functions R_G"
   * - @ref __gnu_cxx::ellint_rj "ellint_rj - Carlson elliptic functions R_J"
   * - @ref __gnu_cxx::ellnome "ellnome - Elliptic nome"
   * - @ref __gnu_cxx::expint "expint - Exponential integrals"
   * - @ref __gnu_cxx::factorial "factorial - Factorials"
   * - @ref __gnu_cxx::fermi_dirac "fermi_dirac - Fermi-Dirac integrals"
   * - @ref __gnu_cxx::fresnel_c "fresnel_c - Fresnel cosine integrals"
   * - @ref __gnu_cxx::fresnel_s "fresnel_s - Fresnel sine integrals"
   * - @ref __gnu_cxx::pgamma "pgamma - Regularized lower incomplete gamma functions"
   * - @ref __gnu_cxx::qgamma "qgamma - Regularized upper incomplete gamma functions"
   * - @ref __gnu_cxx::gegenbauer "gegenbauer - Gegenbauer polynomials"
   * - @ref __gnu_cxx::heuman_lambda "heuman_lambda - Heuman lambda functions"
   * - @ref __gnu_cxx::hurwitz_zeta "hurwitz_zeta - Hurwitz zeta functions"
   * - @ref __gnu_cxx::ibeta "ibeta - Regularized incomplete beta functions"
   * - @ref __gnu_cxx::jacobi "jacobi - Jacobi polynomials"
   * - @ref __gnu_cxx::jacobi_sn "jacobi_sn - Jacobi sine amplitude functions"
   * - @ref __gnu_cxx::jacobi_cn "jacobi_cn - Jacobi cosine amplitude functions"
   * - @ref __gnu_cxx::jacobi_dn "jacobi_dn - Jacobi delta amplitude functions"
   * - @ref __gnu_cxx::jacobi_zeta "jacobi_zeta - Jacobi zeta functions"
   * - @ref __gnu_cxx::lbincoef "lbincoef - Log binomial coefficients"
   * - @ref __gnu_cxx::ldouble_factorial "ldouble_factorial - Log double factorials"
   * - @ref __gnu_cxx::legendre_q "legendre_q - Legendre functions of the second kind"
   * - @ref __gnu_cxx::lfactorial "lfactorial - Log factorials"
   * - @ref __gnu_cxx::lgamma "lgamma - Log gamma for complex arguments"
   * - @ref __gnu_cxx::lpochhammer_lower "lpochhammer_lower - Log lower Pochhammer functions"
   * - @ref __gnu_cxx::lpochhammer "lpochhammer - Log upper Pochhammer functions"
   * - @ref __gnu_cxx::owens_t "owens_t - Owens T functions"
   * - @ref __gnu_cxx::pochhammer_lower "pochhammer_lower - Lower Pochhammer functions"
   * - @ref __gnu_cxx::pochhammer "pochhammer - Upper Pochhammer functions"
   * - @ref __gnu_cxx::psi "psi - Psi or digamma function"
   * - @ref __gnu_cxx::radpoly "radpoly - Radial polynomials"
   * - @ref __gnu_cxx::sinhc "sinhc - Hyperbolic sinus cardinal function"
   * - @ref __gnu_cxx::sinhc_pi "sinhc_pi - "
   * - @ref __gnu_cxx::sinc "sinc - Normalized sinus cardinal function"
   * - @ref __gnu_cxx::sincos "sincos - Sine + cosine function"
   * - @ref __gnu_cxx::sincos_pi "sincos_pi - Reperiodized sine + cosine function"
   * - @ref __gnu_cxx::sin_pi "sin_pi - Reperiodized sine function."
   * - @ref __gnu_cxx::sinh_pi "sinh_pi - Reperiodized hyperbolic sine function."
   * - @ref __gnu_cxx::sinc_pi "sinc_pi - Sinus cardinal function"
   * - @ref __gnu_cxx::sinhint "sinhint - Hyperbolic sine integral"
   * - @ref __gnu_cxx::sinint "sinint - Sine integral"
   * - @ref __gnu_cxx::sph_bessel_i "sph_bessel_i - Spherical regular modified Bessel functions"
   * - @ref __gnu_cxx::sph_bessel_k "sph_bessel_k - Spherical iregular modified Bessel functions"
   * - @ref __gnu_cxx::sph_hankel_1 "sph_hankel_1 - Spherical Hankel functions of the first kind"
   * - @ref __gnu_cxx::sph_hankel_2 "sph_hankel_2 - Spherical Hankel functions of the first kind"
   * - @ref __gnu_cxx::sph_harmonic "sph_harmonic - Spherical"
   * - @ref __gnu_cxx::tan_pi "tan_pi - Reperiodized tangent function."
   * - @ref __gnu_cxx::tanh_pi "tanh_pi - Reperiodized hyperbolic tangent function."
   * - @ref __gnu_cxx::tgamma "tgamma - Gamma for complex arguments"
   * - @ref __gnu_cxx::tgamma "tgamma - Upper incomplete gamma functions"
   * - @ref __gnu_cxx::tgamma_lower "tgamma_lower - Lower incomplete gamma functions"
   * - @ref __gnu_cxx::theta_1 "theta_1 - Exponential theta function 1"
   * - @ref __gnu_cxx::theta_2 "theta_2 - Exponential theta function 2"
   * - @ref __gnu_cxx::theta_3 "theta_3 - Exponential theta function 3"
   * - @ref __gnu_cxx::theta_4 "theta_4 - Exponential theta function 4"
   * - @ref __gnu_cxx::zernike "zernike - Zernike polynomials"
   *
   * @section general General Features
   *
   * @subsection promotion Argument Promotion
   * The arguments suppled to the non-suffixed functions will be promoted
   * according to the following rules:
   * 1. If any argument intended to be floating point is given an integral value
   * That integral value is promoted to double.
   * 2. All floating point arguments are promoted up to the largest floating
   *    point precision among them.
   *
   * @subsection NaN NaN Arguments
   * If any of the floating point arguments supplied to these functions is
   * invalid or NaN (std::numeric_limits<Tp>::quiet_NaN),
   * the value NaN is returned.
   *
   * @section impl Implementation
   *
   * We strive to implement the underlying math with type generic algorithms
   * to the greatest extent possible.  In practice, the functions are thin
   * wrappers that dispatch to function templates. Type dependence is
   * controlled with std::numeric_limits and functions thereof.
   *
   * We don't promote @c float to @c double or @c double to <tt>long double</tt>
   * reflexively.  The goal is for @c float functions to operate more quickly,
   * at the cost of @c float accuracy and possibly a smaller domain of validity.
   * Similaryly, <tt>long double</tt> should give you more dynamic range
   * and slightly more pecision than @c double on many systems.
   *
   * @section testing Testing
   *
   * These functions have been tested against equivalent implementations
   * from the <a href="http://www.gnu.org/software/gsl">
   * Gnu Scientific Library, GSL</a> and
   * <a href="http://www.boost.org/doc/libs/1_60_0/libs/math/doc/html/index.html>Boost</a>
   * and the ratio
   * @f[
   *   \frac{|f - f_{test}|}{|f_{test}|}
   * @f]
   * is generally found to be within 10^-15 for 64-bit double on linux-x86_64 systems
   * over most of the ranges of validity.
   * 
   * @todo Provide accuracy comparisons on a per-function basis for a small
   *       number of targets.
   *
   * @section bibliography General Bibliography
   *
   * @see Abramowitz and Stegun: Handbook of Mathematical Functions,
   * with Formulas, Graphs, and Mathematical Tables
   * Edited by Milton Abramowitz and Irene A. Stegun,
   * National Bureau of Standards  Applied Mathematics Series - 55
   * Issued June 1964, Tenth Printing, December 1972, with corrections
   * Electronic versions of A&S abound including both pdf and navigable html.
   * @see for example  http://people.math.sfu.ca/~cbm/aands/
   *
   * @see The old A&S has been redone as the
   * NIST Digital Library of Mathematical Functions: http://dlmf.nist.gov/
   * This version is far more navigable and includes more recent work.
   *
   * @see An Atlas of Functions: with Equator, the Atlas Function Calculator
   * 2nd Edition, by Oldham, Keith B., Myland, Jan, Spanier, Jerome
   *
   * @see Asymptotics and Special Functions by Frank W. J. Olver,
   * Academic Press, 1974
   *
   * @see Numerical Recipes in C, The Art of Scientific Computing,
   * by William H. Press, Second Ed., Saul A. Teukolsky,
   * William T. Vetterling, and Brian P. Flannery,
   * Cambridge University Press, 1992
   *
   * @see The Special Functions and Their Approximations: Volumes 1 and 2,
   * by Yudell L. Luke, Academic Press, 1969
   */

  /** @} */ // math_spec_func

  /**
   * @defgroup tr29124_math_spec_func C++17/IS29124 Mathematical Special Functions
   * @ingroup math_spec_func
   *
   * A collection of advanced mathematical special functions for C++17
   * and IS29124.
   * @{
   */

  // Associated Laguerre polynomials

  /**
   * Return the associated Laguerre polynomial @f$ L_n^m(x) @f$
   * of order @f$ n @f$, degree @f$ m @f$, and @c float argument @f$ x @f$.
   *
   * @see assoc_laguerre for more details.
   */
  inline float
  assoc_laguerref(unsigned int __n, unsigned int __m, float __x)
  { return __detail::__assoc_laguerre<float>(__n, __m, __x); }

  /**
   * Return the associated Laguerre polynomial @f$ L_n^m(x) @f$
   * of order @f$ n @f$, degree @f$ m @f$ and <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see assoc_laguerre for more details.
   */
  inline long double
  assoc_laguerrel(unsigned int __n, unsigned int __m, long double __x)
  { return __detail::__assoc_laguerre<long double>(__n, __m, __x); }

  /**
   * Return the associated Laguerre polynomial @f$ L_n^m(x) @f$
   * of nonnegative order @f$ n @f$, nonnegative degree @f$ m @f$
   * and real argument @f$ x @f$.
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
   * and @f$ x >= 0 @f$.
   * @see laguerre for details of the Laguerre function of degree @c n
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __n The order of the Laguerre function, <tt>__n >= 0</tt>.
   * @param __m The degree of the Laguerre function, <tt>__m >= 0</tt>.
   * @param __x The argument of the Laguerre function, <tt>__x >= 0</tt>.
   * @throw std::domain_error if <tt>__x < 0</tt>.
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
   * Return the associated Legendre function @f$ P_l^m(x) @f$
   * of degree @f$ l @f$, order @f$ m @f$, and @c float argument @f$ x @f$.
   *
   * @see assoc_legendre for more details.
   */
  inline float
  assoc_legendref(unsigned int __l, unsigned int __m, float __x)
  { return __detail::__assoc_legendre_p<float>(__l, __m, __x); }

  /**
   * Return the associated Legendre function @f$ P_l^m(x) @f$
   * of degree @f$ l @f$, order @f$ m @f$, and @c <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see assoc_legendre for more details.
   */
  inline long double
  assoc_legendrel(unsigned int __l, unsigned int __m, long double __x)
  { return __detail::__assoc_legendre_p<long double>(__l, __m, __x); }

  /**
   * Return the associated Legendre function @f$ P_l^m(x) @f$
   * of degree @f$ l @f$, order @f$ m @f$, and real argument @f$ x @f$.
   *
   * The associated Legendre function is derived from the Legendre function
   * @f$ P_l(x) @f$ by the Rodrigues formula:
   * @f[
   *   P_l^m(x) = (1 - x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)
   * @f]
   * @see legendre for details of the Legendre function of degree @c l
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __l  The degree <tt>__l >= 0</tt>.
   * @param  __m  The order <tt>__m <= l</tt>.
   * @param  __x  The argument, <tt>abs(__x) <= 1</tt>.
   * @throw std::domain_error if <tt>abs(__x) > 1</tt>.
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
   * Return the beta function, @f$ B(a,b) @f$, for @c float parameters @f$ a @f$, @f$ b @f$.
   *
   * @see beta for more details.
   */
  inline float
  betaf(float __a, float __b)
  { return __detail::__beta<float>(__a, __b); }

  /**
   * Return the beta function, @f$B(a,b)@f$, for long double
   * parameters @f$ a @f$, @f$ b @f$.
   *
   * @see beta for more details.
   */
  inline long double
  betal(long double __a, long double __b)
  { return __detail::__beta<long double>(__a, __b); }

  /**
   * Return the beta function, @f$B(a,b)@f$, for real parameters @f$ a @f$, @f$ b @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   * where @f$ a > 0 @f$ and @f$ b > 0 @f$
   *
   * @tparam _Tpa The floating-point type of the parameter @c __a.
   * @tparam _Tpb The floating-point type of the parameter @c __b.
   * @param __a The first argument of the beta function, <tt> __a > 0 </tt>.
   * @param __b The second argument of the beta function, <tt> __b > 0 </tt>.
   * @throw std::domain_error if <tt> __a < 0 </tt> or <tt> __b < 0 </tt>.
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
   * for @c float modulus @f$ k @f$.
   *
   * @see comp_ellint_1 for details.
   */
  inline float
  comp_ellint_1f(float __k)
  { return __detail::__comp_ellint_1<float>(__k); }

  /**
   * Return the complete elliptic integral of the first kind @f$ E(k) @f$
   * for <tt>long double</tt> modulus @f$ k @f$.
   *
   * @see comp_ellint_1 for details.
   */
  inline long double
  comp_ellint_1l(long double __k)
  { return __detail::__comp_ellint_1<long double>(__k); }

  /**
   * Return the complete elliptic integral of the first kind
   * @f$ K(k) @f$ for real modulus @f$ k @f$.
   *
   * The complete elliptic integral of the first kind is defined as
   * @f[
   *   K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   * 					     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   * first kind and the modulus @f$ |k| <= 1 @f$.
   * @see ellint_1 for details of the incomplete elliptic function
   * of the first kind.
   *
   * @tparam _Tp The floating-point type of the modulus @c __k.
   * @param  __k  The modulus, <tt> abs(__k) <= 1 </tt>
   * @throw std::domain_error if <tt> abs(__k) > 1 </tt>.
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
   * for @c float modulus @f$ k @f$.
   *
   * @see comp_ellint_2 for details.
   */
  inline float
  comp_ellint_2f(float __k)
  { return __detail::__comp_ellint_2<float>(__k); }

  /**
   * Return the complete elliptic integral of the second kind @f$ E(k) @f$
   * for <tt>long double</tt> modulus @f$ k @f$.
   *
   * @see comp_ellint_2 for details.
   */
  inline long double
  comp_ellint_2l(long double __k)
  { return __detail::__comp_ellint_2<long double>(__k); }

  /**
   * Return the complete elliptic integral of the second kind @f$ E(k) @f$
   * for real modulus @f$ k @f$.
   *
   * The complete elliptic integral of the second kind is defined as
   * @f[
   *   E(k) = E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   * @f]
   * where @f$ E(k,\phi) @f$ is the incomplete elliptic integral of the
   * second kind and the modulus @f$ |k| <= 1 @f$.
   * @see ellint_2 for details of the incomplete elliptic function
   * of the second kind.
   *
   * @tparam _Tp The floating-point type of the modulus @c __k.
   * @param  __k  The modulus, @c abs(__k) <= 1
   * @throw std::domain_error if @c abs(__k) > 1.
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
   * @f$ \Pi(k,\nu) @f$ for @c float modulus @f$ k @f$.
   *
   * @see comp_ellint_3 for details.
   */
  inline float
  comp_ellint_3f(float __k, float __nu)
  { return __detail::__comp_ellint_3<float>(__k, __nu); }

  /**
   * @brief Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) @f$ for <tt>long double</tt> modulus @f$ k @f$.
   *
   * @see comp_ellint_3 for details.
   */
  inline long double
  comp_ellint_3l(long double __k, long double __nu)
  { return __detail::__comp_ellint_3<long double>(__k, __nu); }

  /**
   * Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ for real modulus @f$ k @f$.
   *
   * The complete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu) = \Pi(k,\nu,\pi/2) = \int_0^{\pi/2}
   * 		     \frac{d\theta}
   * 		   {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   * where @f$ \Pi(k,\nu,\phi) @f$ is the incomplete elliptic integral of the
   * second kind and the modulus @f$ |k| <= 1 @f$.
   * @see ellint_3 for details of the incomplete elliptic function
   * of the third kind.
   *
   * @tparam _Tp The floating-point type of the modulus @c __k.
   * @tparam _Tpn The floating-point type of the argument @c __nu.
   * @param  __k  The modulus, @c abs(__k) <= 1
   * @param  __nu  The argument
   * @throw std::domain_error if @c abs(__k) > 1.
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
   * Return the regular modified Bessel function @f$ I_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_i for setails.
   */
  inline float
  cyl_bessel_if(float __nu, float __x)
  { return __detail::__cyl_bessel_i<float>(__nu, __x); }

  /**
   * Return the regular modified Bessel function @f$ I_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_i for setails.
   */
  inline long double
  cyl_bessel_il(long double __nu, long double __x)
  { return __detail::__cyl_bessel_i<long double>(__nu, __x); }

  /**
   * Return the regular modified Bessel function @f$ I_{\nu}(x) @f$
   * for real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The regular modified cylindrical Bessel function is:
   * @f[
   *  I_{\nu}(x) = i^{-\nu}J_\nu(ix) = \sum_{k=0}^{\infty}
   * 		\frac{(x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
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
   * Return the Bessel function of the first kind @f$ J_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_j for setails.
   */
  inline float
  cyl_bessel_jf(float __nu, float __x)
  { return __detail::__cyl_bessel_j<float>(__nu, __x); }

  /**
   * Return the Bessel function of the first kind @f$ J_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_j for setails.
   */
  inline long double
  cyl_bessel_jl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_j<long double>(__nu, __x); }

  /**
   * Return the Bessel function @f$ J_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The cylindrical Bessel function is:
   * @f[
   *    J_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
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
   * Return the irregular modified Bessel function @f$ K_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_k for setails.
   */
  inline float
  cyl_bessel_kf(float __nu, float __x)
  { return __detail::__cyl_bessel_k<float>(__nu, __x); }

  /**
   * Return the irregular modified Bessel function @f$ K_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_k for setails.
   */
  inline long double
  cyl_bessel_kl(long double __nu, long double __x)
  { return __detail::__cyl_bessel_k<long double>(__nu, __x); }

  /**
   * Return the irregular modified Bessel function @f$ K_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * The irregular modified Bessel function is defined by:
   * @f[
   * 	K_{\nu}(x) = \frac{\pi}{2}
   * 		     \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   * For negative argument we have simply:
   * @f[
   * 	K_{-\nu}(x) = K_{\nu}(x)
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
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
   * of <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * @see cyl_neumann for setails.
   */
  inline long double
  cyl_neumannl(long double __nu, long double __x)
  { return __detail::__cyl_neumann_n<long double>(__nu, __x); }

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The Neumann function is defined by:
   * @f[
   *    N_{\nu}(x) = \frac{J_{\nu}(x) \cos \nu\pi - J_{-\nu}(x)}
   *                      {\sin \nu\pi}
   * @f]
   * where @f$ x >= 0 @f$ and for integral order @f$ \nu = n @f$
   * a limit is taken: @f$ lim_{\nu \to n} @f$.
   *
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
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
   * for <tt>long double</tt> modulus @f$ k @f$ and angle @f$ \phi @f$.
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
   * @tparam _Tp The floating-point type of the modulus @c __k.
   * @tparam _Tpp The floating-point type of the angle @c __phi.
   * @param  __k  The modulus, <tt> abs(__k) <= 1 </tt>
   * @param  __phi  The integral limit argument in radians
   * @throw std::domain_error if <tt> abs(__k) > 1 </tt>.
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
   * @tparam _Tp The floating-point type of the modulus @c __k.
   * @tparam _Tpp The floating-point type of the angle @c __phi.
   * @param  __k  The modulus, <tt> abs(__k) <= 1 </tt>
   * @param  __phi  The integral limit argument in radians
   * @return  The elliptic function of the second kind.
   * @throw std::domain_error if <tt> abs(__k) > 1 </tt>.
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
   * @tparam _Tp The floating-point type of the modulus @c __k.
   * @tparam _Tpn The floating-point type of the argument @c __nu.
   * @tparam _Tpp The floating-point type of the angle @c __phi.
   * @param  __k  The modulus, <tt> abs(__k) <= 1 </tt>
   * @param  __nu  The second argument
   * @param  __phi  The integral limit argument in radians
   * @return  The elliptic function of the third kind.
   * @throw std::domain_error if <tt> abs(__k) > 1 </tt>.
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
   * Return the exponential integral @f$ Ei(x) @f$ for @c float argument @f$ x @f$.
   *
   * @see expint for details.
   */
  inline float
  expintf(float __x)
  { return __detail::__expint<float>(__x); }

  /**
   * Return the exponential integral @f$ Ei(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see expint for details.
   */
  inline long double
  expintl(long double __x)
  { return __detail::__expint<long double>(__x); }

  /**
   * Return the exponential integral @f$ Ei(x) @f$ for @c real argument @f$ x @f$.
   *
   * The exponential integral is given by
   * \f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * \f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
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
   * Return the Hermite polynomial @f$ H_n(x) @f$ of nonnegative order n
   * and float argument @f$ x @f$.
   *
   * @see hermite for details.
   */
  inline float
  hermitef(unsigned int __n, float __x)
  { return __detail::__poly_hermite<float>(__n, __x); }

  /**
   * Return the Hermite polynomial @f$ H_n(x) @f$ of nonnegative order n
   * and <tt>long double</tt> argument @f$ x @f$.
   *
   * @see hermite for details.
   */
  inline long double
  hermitel(unsigned int __n, long double __x)
  { return __detail::__poly_hermite<long double>(__n, __x); }

  /**
   * Return the Hermite polynomial @f$ H_n(x) @f$ of order n
   * and @c real argument @f$ x @f$.
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
   * @tparam _Tp The floating-point type of the argument @c __x.
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
   * Returns the Laguerre polynomial @f$ L_n(x) @f$ of nonnegative degree @c n
   * and @c float argument  @f$ x >= 0 @f$.
   *
   * @see laguerre for more details.
   */
  inline float
  laguerref(unsigned int __n, float __x)
  { return __detail::__laguerre<float>(__n, __x); }

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$ of nonnegative degree @c n
   * and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see laguerre for more details.
   */
  inline long double
  laguerrel(unsigned int __n, long double __x)
  { return __detail::__laguerre<long double>(__n, __x); }

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$
   * of nonnegative degree @f$ n @f$ and real argument @f$ x >= 0 @f$.
   *
   * The Laguerre polynomial is defined by:
   * @f[
   * 	 L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __n The nonnegative order
   * @param __x The argument <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
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
   * Return the Legendre polynomial @f$ P_l(x) @f$ of nonnegative
   * degree @f$ l @f$ and @c float argument @f$ |x| <= 0 @f$.
   *
   * @see legendre for more details.
   */
  inline float
  legendref(unsigned int __l, float __x)
  { return __detail::__poly_legendre_p<float>(__l, __x); }

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of nonnegative
   * degree @f$ l @f$ and <tt>long double</tt> argument @f$ |x| <= 0 @f$.
   *
   * @see legendre for more details.
   */
  inline long double
  legendrel(unsigned int __l, long double __x)
  { return __detail::__poly_legendre_p<long double>(__l, __x); }

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of nonnegative
   * degree @f$ l @f$ and real argument @f$ |x| <= 0 @f$.
   *
   * The Legendre function of order @f$ l @f$ and argument @f$ x @f$,
   * @f$ P_l(x) @f$, is defined by:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __l The degree @f$ l >= 0 @f$
   * @param __x The argument @c abs(__x) <= 1
   * @throw std::domain_error if @c abs(__x) > 1
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
   * for <tt>long double</tt> argument @f$ s @f$.
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
   * 	\zeta(s) = \sum_{k=1}^{\infty} k^{-s} \mbox{ for } s > 1
   * @f]
   * and
   * @f[
   * 	\zeta(s) = \frac{1}{1-2^{1-s}}\sum_{k=1}^{\infty}(-1)^{k-1}k^{-s}
   *              \mbox{ for } 0 <= s <= 1
   * @f]
   * For s < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = 2^s \pi^{s-1} \sin(\frac{\pi s}{2}) \Gamma(1-s) \zeta(1-s)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __s.
   * @param __s The argument <tt> s != 1 </tt>
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
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel for more details.
   */
  inline float
  sph_besself(unsigned int __n, float __x)
  { return __detail::__sph_bessel<float>(__n, __x); }

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel for more details.
   */
  inline long double
  sph_bessell(unsigned int __n, long double __x)
  { return __detail::__sph_bessel<long double>(__n, __x); }

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  j_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} J_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __n  The integral order <tt> n >= 0 </tt>
   * @param  __x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
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
   * Return the spherical Legendre function of nonnegative integral
   * degree @f$ l @f$ and order @f$ m @f$ and float angle @f$ \theta @f$ in radians.
   *
   * @see sph_legendre for details.
   */
  inline float
  sph_legendref(unsigned int __l, unsigned int __m, float __theta)
  { return __detail::__sph_legendre<float>(__l, __m, __theta); }

  /**
   * Return the spherical Legendre function of nonnegative integral
   * degree @f$ l @f$ and order @f$ m @f$ and <tt>long double</tt> angle @f$ \theta @f$
   * in radians.
   *
   * @see sph_legendre for details.
   */
  inline long double
  sph_legendrel(unsigned int __l, unsigned int __m, long double __theta)
  { return __detail::__sph_legendre<long double>(__l, __m, __theta); }

  /**
   * Return the spherical Legendre function of nonnegative integral
   * degree @f$ l @f$ and order @f$ m @f$ and real angle @f$ \theta @f$ in radians.
   *
   * The spherical Legendre function is defined by
   * @f[
   *  Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   *                              \frac{(l-m)!}{(l+m)!}]
   *                   P_l^m(\cos\theta) \exp^{im\phi}
   * @f]
   *
   * @tparam _Tp The floating-point type of the angle @c __theta.
   * @param __l The order <tt> __l >= 0 </tt>
   * @param __m The degree <tt> __m >= 0 </tt> and <tt> __m <= __l </tt>
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
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_neumann for details.
   */
  inline float
  sph_neumannf(unsigned int __n, float __x)
  { return __detail::__sph_neumann<float>(__n, __x); }

  /**
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and <tt>long double</tt> @f$ x >= 0 @f$.
   *
   * @see sph_neumann for details.
   */
  inline long double
  sph_neumannl(unsigned int __n, long double __x)
  { return __detail::__sph_neumann<long double>(__n, __x); }

  /**
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and real argument @f$ x >= 0 @f$.
   *
   * The spherical Neumann function is defined by
   * @f[
   *    n_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} N_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __n  The integral order <tt> n >= 0 </tt>
   * @param  __x  The real argument <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sph_neumann(unsigned int __n, _Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return __detail::__sph_neumann<__type>(__n, __x);
    }

  /** @} */ // tr29124_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @defgroup gnu_math_spec_func GNU Extended Mathematical Special Functions
   * @ingroup math_spec_func
   *
   * An extended collection of advanced mathematical special functions for GNU.
   * @{
   */

  // Confluent hypergeometric functions

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of @c float numeratorial parameter @f$ a @f$, denominatorial parameter @f$ c @f$,
   * and argument @f$ x @f$.
   *
   * @see conf_hyperg for details.
   */
  inline float
  conf_hypergf(float __a, float __c, float __x)
  { return std::__detail::__conf_hyperg<float>(__a, __c, __x); }

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of <tt>long double</tt> numeratorial parameter @f$ a @f$,
   * denominatorial parameter @f$ c @f$, and argument @f$ x @f$.
   *
   * @see conf_hyperg for details.
   */
  inline long double
  conf_hypergl(long double __a, long double __c, long double __x)
  { return std::__detail::__conf_hyperg<long double>(__a, __c, __x); }

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of real numeratorial parameter @f$ a @f$, denominatorial parameter @f$ c @f$,
   * and argument @f$ x @f$.
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
   * of @ float numeratorial parameters @f$ a @f$ and @f$ b @f$,
   * denominatorial parameter @f$ c @f$, and argument @f$ x @f$.
   *
   * @see hyperg for details.
   */
  inline float
  hypergf(float __a, float __b, float __c, float __x)
  { return std::__detail::__hyperg<float>(__a, __b, __c, __x); }

  /**
   * Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * of <tt>long double</tt> numeratorial parameters @f$ a @f$ and @f$ b @f$,
   * denominatorial parameter @f$ c @f$, and argument @f$ x @f$.
   *
   * @see hyperg for details.
   */
  inline long double
  hypergl(long double __a, long double __b, long double __c, long double __x)
  { return std::__detail::__hyperg<long double>(__a, __b, __c, __x); }

  /**
   * Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * of real numeratorial parameters @f$ a @f$ and @f$ b @f$,
   * denominatorial parameter @f$ c @f$, and argument @f$ x @f$.
   *
   * The hypergeometric function is defined by
   * @f[
   *    {}_2F_1(a,b;c;x) = \sum_{n=0}^{\infty} \frac{(a)_n (b)_n x^n}{(c)_n n!}
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
   * of @c float numeratorial parameter @f$ c @f$ and argument @f$ x @f$.
   *
   * @see conf_hyperg_lim for details.
   */
  inline float
  conf_hyperg_limf(float __c, float __x)
  { return std::__detail::__conf_hyperg_lim<float>(__c, __x); }

  /**
   * Return the confluent hypergeometric limit function @f$ {}_0F_1(;c;x) @f$
   * of <tt>long double</tt> numeratorial parameter @f$ c @f$ and argument @f$ x @f$.
   *
   * @see conf_hyperg_lim for details.
   */
  inline long double
  conf_hyperg_liml(long double __c, long double __x)
  { return std::__detail::__conf_hyperg_lim<long double>(__c, __x); }

  /**
   * Return the confluent hypergeometric limit function @f$ {}_0F_1(;c;x) @f$
   * of real numeratorial parameter @f$ c @f$ and argument @f$ x @f$.
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

  // Sinus cardinal functions

  /**
   * Return the sinus cardinal function @f$ sinc_\pi(x) @f$
   * for @c float argument @c __x.
   *
   * @see sinc_pi for details.
   */
  inline float
  sincf(float __x)
  { return std::__detail::__sinc<float>(__x); }

  /**
   * Return the sinus cardinal function @f$ sinc_\pi(x) @f$
   * for <tt>long double</tt> argument @c __x.
   *
   * @see sinc_pi for details.
   */
  inline long double
  sincl(long double __x)
  { return std::__detail::__sinc<long double>(__x); }

  /**
   * Return the sinus cardinal function @f$ sinc_\pi(x) @f$
   * for real argument @c __x.
   * The sinus cardinal function is defined by:
   * @f[
   *    sinc(x) = \frac{sin(x)}{x}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sinc(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sinc<__type>(__x);
    }

  // Normalized sinus cardinal functions

  /**
   * Return the reperiodized sinus cardinal function @f$ sinc(x) @f$
   * for @c float argument @c __x.
   *
   * @see sinc for details.
   */
  inline float
  sinc_pif(float __x)
  { return std::__detail::__sinc_pi<float>(__x); }

  /**
   * Return the reperiodized sinus cardinal function @f$ sinc(x) @f$
   * for <tt>long double</tt> argument @c __x.
   *
   * @see sinc for details.
   */
  inline long double
  sinc_pil(long double __x)
  { return std::__detail::__sinc_pi<long double>(__x); }

  /**
   * Return the reperiodized sinus cardinal function @f$ sinc(x) @f$
   * for real argument @c __x.
   * The normalized sinus cardinal function is defined by:
   * @f[
   *    sinc_\pi(x) = \frac{sin(\pi x)}{\pi x}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sinc_pi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sinc_pi<__type>(__x);
    }

  // Logarithmic integrals

  /**
   * Return the logarithmic integral of argument @f$ x @f$.
   *
   * @see logint for details.
   */
  inline float
  logintf(float __x)
  { return std::__detail::__logint<float>(__x); }

  /**
   * Return the logarithmic integral of argument @f$ x @f$.
   *
   * @see logint for details.
   */
  inline long double
  logintl(long double __x)
  { return std::__detail::__logint<long double>(__x); }

  /**
   * Return the logarithmic integral of argument @f$ x @f$.
   *
   * The logarithmic integral is defined by
   * @f[
   *    li(x) = \int_0^x \frac{dt}{ln(t)}
   * @f]
   *
   * @param __x The real upper integration limit
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    logint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__logint<__type>(__x);
    }

  // Sine integrals

  /**
   * Return the sine integral @f$ Si(x) @f$ of @c float argument @f$ x @f$.
   *
   * @see sinint for details.
   */
  inline float
  sinintf(float __x)
  { return std::__detail::__sincosint<float>(__x).first; }

  /**
   * Return the sine integral @f$ Si(x) @f$ of <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see sinint for details.
   */
  inline long double
  sinintl(long double __x)
  { return std::__detail::__sincosint<long double>(__x).first; }

  /**
   * Return the sine integral @f$ Si(x) @f$ of real argument @f$ x @f$.
   *
   * The sine integral is defined by
   * @f[
   *    Si(x) = \int_0^x \frac{sin(t)}{t}dt
   * @f]
   *
   * @param __x The real upper integration limit
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sinint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sincosint<__type>(__x).first;
    }

  // Cosine integrals

  /**
   * Return the cosine integral @f$ Ci(x) @f$ of @c float argument @f$ x @f$.
   *
   * @see cosint for details.
   */
  inline float
  cosintf(float __x)
  { return std::__detail::__sincosint<float>(__x).second; }

  /**
   * Return the cosine integral @f$ Ci(x) @f$ of <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see cosint for details.
   */
  inline long double
  cosintl(long double __x)
  { return std::__detail::__sincosint<long double>(__x).second; }

  /**
   * Return the cosine integral @f$ Ci(x) @f$ of real argument @f$ x @f$.
   *
   * The cosine integral is defined by
   * @f[
   *    Ci(x) = -\int_x^\infty \frac{cos(t)}{t}dt
   *     = \gamma_E + ln(x) + \int_0^x \frac{cos(t)-1}{t}dt
   * @f]
   *
   * @param __x The real upper integration limit
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    cosint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sincosint<__type>(__x).second;
    }

  // Hyperbolic sine integrals

  /**
   * Return the hyperbolic sine integral of @c float argument @f$ x @f$.
   *
   * @see sinhint for details.
   */
  inline float
  sinhintf(float __x)
  { return std::__detail::__sinhint<float>(__x); }

  /**
   * Return the hyperbolic sine integral @f$ Shi(x) @f$ of <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see sinhint for details.
   */
  inline long double
  sinhintl(long double __x)
  { return std::__detail::__sinhint<long double>(__x); }

  /**
   * Return the hyperbolic sine integral @f$ Shi(x) @f$
   * of real argument @f$ x @f$.
   *
   * The hyperbolic sine integral is defined by
   * @f[
   *    Shi(x) = \int_0^x \frac{\sinh(t)}{t}dt
   * @f]
   *
   * @tparam _Tp The type of the real argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sinhint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sinhint<__type>(__x);
    }

  // Hyperbolic cosine integrals

  /**
   * Return the hyperbolic cosine integral of @c float argument @f$ x @f$.
   *
   * @see coshint for details.
   */
  inline float
  coshintf(float __x)
  { return std::__detail::__coshint<float>(__x); }

  /**
   * Return the hyperbolic cosine integral @f$ Chi(x) @f$
   * of <tt>long double</tt> argument @f$ x @f$.
   *
   * @see coshint for details.
   */
  inline long double
  coshintl(long double __x)
  { return std::__detail::__coshint<long double>(__x); }

  /**
   * Return the hyperbolic cosine integral @f$ Chi(x) @f$
   * of real argument @f$ x @f$.
   *
   * The hyperbolic cosine integral is defined by
   * @f[
   *    Chi(x) = -\int_x^\infty \frac{cosh(t)}{t}dt
   *     = \gamma_E + ln(x) + \int_0^x \frac{cosh(t)-1}{t}dt
   * @f]
   *
   * @tparam _Tp The type of the real argument
   * @param __x The real argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    coshint(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__coshint<__type>(__x);
    }

  // Slots for Jacobi elliptic function tuple.
  enum
  {
    _GLIBCXX_JACOBI_SN,
    _GLIBCXX_JACOBI_CN,
    _GLIBCXX_JACOBI_DN
  };

  // Jacobi elliptic sine amplitude functions.

  /**
   * Return the Jacobi elliptic @f$ sn(k,u) @f$ integral
   * of @c float modulus @f$ k @f$ and argument @f$ u @f$.
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
   * Return the Jacobi elliptic @f$ sn(k,u) @f$ integral
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ u @f$.
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
   * Return the Jacobi elliptic @f$ sn(k,u) @f$ integral
   * of real modulus @f$ k @f$ and argument @f$ u @f$.
   *
   * The Jacobi elliptic @c sn integral is defined by
   * @f[
   *    \sin(\phi) = sn(k, F(k,\phi))
   * @f]
   * where @f$ F(k,\phi) @f$ is the elliptic integral of the first kind.
   *
   * @tparam _Kp The type of the real modulus
   * @tparam _Up The type of the real argument
   * @param __k The real modulus
   * @param __u The real argument
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_sn(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::get<_GLIBCXX_JACOBI_SN>
		(std::__detail::__jacobi_sncndn<__type>(__k, __u));
    }

  // Jacobi elliptic cosine amplitude functions.

  /**
   * Return the Jacobi elliptic @f$ cn(k,u) @f$ integral
   * of @c float modulus @f$ k @f$ and argument @f$ u @f$.
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
   * Return the Jacobi elliptic @f$ cn(k,u) @f$ integral
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ u @f$.
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
   * Return the Jacobi elliptic @f$ cn(k,u) @f$ integral
   * of real modulus @f$ k @f$ and argument @f$ u @f$.
   *
   * The Jacobi elliptic @c cn integral is defined by
   * @f[
   *    \cos(\phi) = cn(k, F(k,\phi))
   * @f]
   * where @f$ F(k,\phi) @f$ is the elliptic integral of the first kind.
   *
   * @tparam _Kp The type of the real modulus
   * @tparam _Up The type of the real argument
   * @param __k The real modulus
   * @param __u The real argument
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_cn(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::get<_GLIBCXX_JACOBI_CN>
		(std::__detail::__jacobi_sncndn<__type>(__k, __u));
    }

  // Jacobi elliptic delta amplitude functions.

  /**
   * Return the Jacobi elliptic @f$ dn(k,u) @f$ integral
   * of @c float modulus @f$ k @f$ and argument @f$ u @f$.
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
   * Return the Jacobi elliptic @f$ dn(k,u) @f$ integral
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ u @f$.
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
   * Return the Jacobi elliptic @f$ dn(k,u) @f$ integral
   * of real modulus @f$ k @f$ and argument @f$ u @f$.
   *
   * The Jacobi elliptic @c dn integral is defined by
   * @f[
   *    \sqrt{1 - k^2\sin(\phi)} = dn(k, F(k,\phi))
   * @f]
   * where @f$ F(k,\phi) @f$ is the elliptic integral of the first kind.
   *
   * @tparam _Kp The type of the real modulus
   * @tparam _Up The type of the real argument
   * @param __k The real modulus
   * @param __u The real argument
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_dn(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::get<_GLIBCXX_JACOBI_DN>
		(std::__detail::__jacobi_sncndn<__type>(__k, __u));
    }

  // Chebyshev polynomials of the first kind

  /**
   * Return the Chebyshev polynomials of the first kind @f$ T_n(x) @f$
   * of non-negative order @f$ n @f$ and @c float argument @f$ x @f$.
   *
   * @see chebyshev_t for details.
   */
  inline float
  chebyshev_tf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_t<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the first kind @f$ T_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * @see chebyshev_t for details.
   */
  inline long double
  chebyshev_tl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_t<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomial of the first kind @f$ T_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the first kind is defined by:
   * @f[
   *    T_n(x) = \cos(n \theta)
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    chebyshev_t(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__chebyshev_t<__type>(__n, __x);
    }

  // Chebyshev polynomials of the second kind

  /**
   * Return the Chebyshev polynomials of the second kind @f$ U_n(x) @f$
   * of non-negative order @f$ n @f$ and @c float argument @f$ x @f$.
   *
   * @see chebyshev_u for details.
   */
  inline float
  chebyshev_uf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_u<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the second kind  @f$ U_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * @see chebyshev_u for details.
   */
  inline long double
  chebyshev_ul(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_u<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomial of the second kind @f$ U_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the second kind is defined by:
   * @f[
   *    U_n(x) = \frac{\sin \left[(n+1)\theta \right]}{\sin(\theta)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    chebyshev_u(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__chebyshev_u<__type>(__n, __x);
    }

  // Chebyshev polynomials of the third kind

  /**
   * Return the Chebyshev polynomials of the third kind @f$ V_n(x) @f$
   * of non-negative order @f$ n @f$ and @c float argument @f$ x @f$.
   *
   * @see chebyshev_v for details.
   */
  inline float
  chebyshev_vf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_v<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the third kind @f$ V_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * @see chebyshev_v for details.
   */
  inline long double
  chebyshev_vl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_v<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomial of the third kind @f$ V_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the third kind is defined by:
   * @f[
   *    V_n(x) = \frac{\cos \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\cos \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    chebyshev_v(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__chebyshev_v<__type>(__n, __x);
    }

  // Chebyshev polynomials of the fourth kind

  /**
   * Return the Chebyshev polynomials of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @f$ n @f$ and @c float argument @f$ x @f$.
   *
   * @see chebyshev_w for details.
   */
  inline float
  chebyshev_wf(unsigned int __n, float __x)
  { return std::__detail::__chebyshev_w<float>(__n, __x); }

  /**
   * Return the Chebyshev polynomials of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * @see chebyshev_w for details.
   */
  inline long double
  chebyshev_wl(unsigned int __n, long double __x)
  { return std::__detail::__chebyshev_w<long double>(__n, __x); }

  /**
   * Return the Chebyshev polynomial of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the fourth kind is defined by:
   * @f[
   *    W_n(x) = \frac{\sin \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\sin \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    chebyshev_w(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__chebyshev_w<__type>(__n, __x);
    }

  // Jacobi polynomials

  /**
   * Return the Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$
   * of degree @f$ n @f$ and @c float orders @f$ \alpha, \beta > -1 @f$
   * and argument @f$ x @f$.
   *
   * @see jacobi for details.
   */
  inline float
  jacobif(unsigned __n, float __alpha, float __beta, float __x)
  { return std::__detail::__poly_jacobi<float>(__n, __alpha, __beta, __x); }

  /**
   * Return the Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$
   * of degree @f$ n @f$ and @c <tt>long double</tt> orders @f$ \alpha, \beta > -1 @f$
   * and argument @f$ x @f$.
   *
   * @see jacobi for details.
   */
  inline long double
  jacobil(unsigned __n, long double __alpha, long double __beta, long double __x)
  { return std::__detail::__poly_jacobi<long double>(__n, __alpha, __beta, __x); }

  /**
   * Return the Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$
   * of degree @f$ n @f$ and @c float orders @f$ \alpha, \beta > -1 @f$
   * and argument @f$ x @f$.
   *
   * The Jacobi polynomials are generated by a three-term recursion relation:
   * @f[
   *   2 n(\alpha + \beta + n) (\alpha + \beta + 2n - 2)
   *         P^{(\alpha, \beta)}_{n}(x)
   *     = (\alpha + \beta + 2n - 1)
   *       ((\alpha^2 - \beta^2)
   *        + x(\alpha + \beta + 2n - 2)(\alpha + \beta + 2n))
   *         P^{(\alpha, \beta)}_{n-1}(x)
   *     - 2 (\alpha + n - 1)(\beta + n - 1)(\alpha + \beta + 2n)
   *         P^{(\alpha, \beta)}_{n-2}(x)
   * @f]
   * where @f$ P_0^{(\alpha,\beta)}(x) = 1 @f$ and
   * @f$ P_1^{(\alpha,\beta)}(x)
   *      = ((\alpha-\beta) + (2 + (\alpha+\beta)) * x) / 2 @f$.
   *
   * @tparam _Talpha The real type of the order @f$ \alpha @f$
   * @tparam _Tbeta The real type of the order @f$ \beta @f$
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral degree
   * @param __alpha The real order
   * @param __beta The real order
   * @param __x The real argument
   */
  template<typename _Talpha, typename _Tbeta, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Talpha, _Tbeta, _Tp>
    jacobi(unsigned __n, _Talpha __alpha, _Tbeta __beta, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Talpha, _Tbeta, _Tp>;
      return std::__detail::__poly_jacobi<__type>(__n, __alpha, __beta, __x);
    }

  // Gegenbauer polynomials

  /**
   * Return the Gegenbauer polynomial @f$ C_n^{\alpha}(x) @f$ of degree @c n
   * and @c float order @f$ \alpha > -1/2, \alpha \neq 0 @f$ and argument @f$ x @f$.
   *
   * @see gegenbauer for details.
   */
  inline float
  gegenbauerf(unsigned int __n, float __alpha, float __x)
  { return std::__detail::__gegenbauer_poly<float>(__n, __alpha, __x); }

  /**
   * Return the Gegenbauer polynomial @f$ C_n^{\alpha}(x) @f$ of degree @c n
   * and <tt>long double</tt> order @f$ \alpha > -1/2, \alpha \neq 0 @f$
   * and argument @f$ x @f$.
   *
   * @see gegenbauer for details.
   */
  inline long double
  gegenbauerl(unsigned int __n, long double __alpha, long double __x)
  { return std::__detail::__gegenbauer_poly<long double>(__n, __alpha, __x); }

  /**
   * Return the Gegenbauer polynomial @f$ C_n^{\alpha}(x) @f$ of degree @c n
   * and real order @f$ \alpha > -1/2, \alpha \neq 0 @f$ and argument @f$ x @f$.
   *
   * The Gegenbauer polynomials are generated by a three-term recursion relation:
   * @f[
   *    C_n^{\alpha}(x) = \frac{1}{n}\left[ 2x(n+\alpha-1)C_{n-1}^{\alpha}(x)
   *                   - (n+2\alpha-2)C_{n-2}^{\alpha}(x) \right]
   * @f]
   * and @f$ C_0^{\alpha}(x) = 1 @f$, @f$ C_1^{\alpha}(x) = 2\alpha x @f$.
   *
   * @tparam _Talpha The real type of the order
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral degree
   * @param __alpha The real order
   * @param __x The real argument
   */
  template<typename _Talpha, typename _Tp>
    inline typename __gnu_cxx::__promote_fp_t<_Talpha, _Tp>
    gegenbauer(unsigned int __n, _Talpha __alpha, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Talpha, _Tp>;
      return std::__detail::__gegenbauer_poly<__type>(__n, __alpha, __x);
    }

  // Zernike polynomials

  /**
   * Return the Zernicke polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative degree @f$ n @f$, signed order @f$ m @f$,
   * and real radial argument @f$ \rho @f$ and azimuthal angle @f$ \phi @f$.
   *
   * @see zernike for details.
   */
  inline float
  zernikef(unsigned int __n, int __m, float __rho, float __phi)
  { return std::__detail::__zernike<float>(__n, __m, __rho, __phi); }

  /**
   * Return the Zernicke polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative degree @f$ n @f$, signed order @f$ m @f$,
   * and real radial argument @f$ \rho @f$ and azimuthal angle @f$ \phi @f$.
   *
   * @see zernike for details.
   */
  inline long double
  zernikel(unsigned int __n, int __m, long double __rho, long double __phi)
  { return std::__detail::__zernike<long double>(__n, __m, __rho, __phi); }

  /**
   * Return the Zernicke polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative degree @f$ n @f$, signed order @f$ m @f$,
   * and real radial argument @f$ \rho @f$ and azimuthal angle @f$ \phi @f$.
   *
   * The even Zernicke polynomials are defined by:
   * @f[
   *    Z_n^m(\rho,\phi) = R_n^m(\rho)\cos(m\phi)
   * @f]
   * and the odd Zernicke polynomials are defined by:
   * @f[
   *    Z_n^{-m}(\rho,\phi) = R_n^m(\rho)\sin(m\phi)
   * @f]
   * for non-negative degree @f$ m @f$ and @f$ m <= n @f$
   * and where @f$ R_n^m(\rho) @f$ is the radial polynomial (@see radpoly).
   *
   * @tparam _Trho The real type of the radial coordinate
   * @tparam _Tphi The real type of the azimuthal angle
   * @param __n The non-negative degree.
   * @param __m The (signed) azimuthal order
   * @param __rho The radial coordinate
   * @param __phi The azimuthal angle
   */
  template<typename _Trho, typename _Tphi>
    inline __gnu_cxx::__promote_fp_t<_Trho, _Tphi>
    zernike(unsigned int __n, int __m, _Trho __rho, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Trho, _Tphi>;
      return std::__detail::__zernike<__type>(__n, __m, __rho, __phi);
    }

  // Radial polynomials

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @f$ n @f$, order @f$ m <= n @f$, and @c float radial
   * argument @f$ \rho @f$.
   *
   * @see radpoly for details.
   */
  inline float
  radpolyf(unsigned int __n, unsigned int __m, float __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @f$ n @f$, order @f$ m <= n @f$, and <tt>long double</tt> radial
   * argument @f$ \rho @f$.
   *
   * @see radpoly for details.
   */
  inline long double
  radpolyl(unsigned int __n, unsigned int __m, long double __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @f$ n @f$, order @f$ m <= n @f$, and real radial
   * argument @f$ \rho @f$.
   *
   * The radial polynomials are defined by 
   * @f[
   *     R_n^m(\rho) = \sum_{k=0}^{\frac{n-m}{2}}
   *       \frac{(-1)^k(n-k)!}{k!(\frac{n+m}{2}-k)!(\frac{n-m}{2}-k)!}
   *       \rho^{n-2k}
   * @f]
   * for @f$ n - m @f$ even and identically 0 for @f$ n - m @f$ odd.
   * The radial polynomials can be related to the Jacobi polynomials:
   * @f[
   *    R_n^m(\rho) = 
   * @f]
   * @see jacobi for details on the Jacobi polynomials.
   *
   * @tparam _Tp The real type of the radial coordinate
   * @param __n The non-negative degree.
   * @param __m The non-negative azimuthal order
   * @param __rho The radial argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    radpoly(unsigned int __n, unsigned int __m, _Tp __rho)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__poly_radial_jacobi<__type>(__n, __m, __rho);
    }

  // Unnormalized hyperbolic sinus cardinal functions

  /**
   * Return the hyperbolic sinus cardinal function @f$ sinhc_\pi(x) @f$
   * for @c float argument @c __x.
   *
   * @see sinhc_pi for details.
   */
  inline float
  sinhc_pif(float __x)
  { return std::__detail::__sinhc_pi<float>(__x); }

  /**
   * Return the hyperbolic sinus cardinal function @f$ sinhc_\pi(x) @f$
   * for <tt>long double</tt> argument @c __x.
   *
   * @see sinhc_pi for details.
   */
  inline long double
  sinhc_pil(long double __x)
  { return std::__detail::__sinhc_pi<long double>(__x); }

  /**
   * Return the hyperbolic sinus cardinal function @f$ sinhc_\pi(x) @f$
   * for real argument @c __x.
   * The sinus cardinal function is defined by:
   * @f[
   *    sinhc_\pi(x) = \frac{\sinh(x)}{x}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sinhc_pi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sinhc_pi<__type>(__x);
    }

  // Normalized hyperbolic sinus cardinal functions

  /**
   * Return the normalized hyperbolic sinus cardinal function @f$ sinhc(x) @f$
   * for @c float argument @c __x.
   *
   * @see sinhc for details.
   */
  inline float
  sinhcf(float __x)
  { return std::__detail::__sinhc<float>(__x); }

  /**
   * Return the normalized hyperbolic sinus cardinal function @f$ sinhc(x) @f$
   * for <tt>long double</tt> argument @c __x.
   *
   * @see sinhc for details.
   */
  inline long double
  sinhcl(long double __x)
  { return std::__detail::__sinhc<long double>(__x); }

  /**
   * Return the normalized hyperbolic sinus cardinal function @f$ sinhc(x) @f$
   * for real argument @c __x.
   * The normalized hyperbolic sinus cardinal function is defined by:
   * @f[
   *    sinhc(x) = \frac{\sinh(\pi x)}{\pi x}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sinhc(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sinhc<__type>(__x);
    }

  // Cylindrical Hankel functions of the first kind

  /**
   * Return the cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_1 for details.
   */
  inline std::complex<float>
  cyl_hankel_1f(float __nu, float __z)
  { return std::__detail::__cyl_hankel_1<float>(__nu, __z); }

  /**
   * Return the cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_1 for details.
   */
  inline std::complex<long double>
  cyl_hankel_1l(long double __nu, long double __z)
  { return std::__detail::__cyl_hankel_1<long double>(__nu, __z); }

  /**
   * Return the cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_n(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The cylindrical Hankel function of the first kind is defined by:
   * @f[
   *    H^{(1)}_\nu(x) = \left(\frac{\pi}{2x} \right) ^{1/2}
   *       \left[ J_{n+1/2}(x) + iN_{n+1/2}(x) \right]
   * @f]
   * where @f$ J_\nu(x) @f$ and @f$ N_\nu(x) @f$ are the cylindrical Bessel
   * and Neumann functions respectively (@see cyl_bessel and cyl_neumann).
   *
   * @tparam _Tp The real type of the argument
   * @param __nu The real order
   * @param __z The real argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tpnu, _Tp>>
    cyl_hankel_1(_Tpnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_1<__type>(__nu, __z);
    }

  // Cylindrical Hankel functions of the second kind

  /**
   * Return the cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_2 for details.
   */
  inline std::complex<float>
  cyl_hankel_2f(float __nu, float __z)
  { return std::__detail::__cyl_hankel_2<float>(__nu, __z); }

  /**
   * Return the cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_2 for details.
   */
  inline std::complex<long double>
  cyl_hankel_2l(long double __nu, long double __z)
  { return std::__detail::__cyl_hankel_2<long double>(__nu, __z); }

  /**
   * Return the cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_n(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The cylindrical Hankel function of the second kind is defined by:
   * @f[
   *    H^{(2)}_\nu(x) = \left(\frac{\pi}{2x} \right) ^{1/2}
   *       \left[ J_{n+1/2}(x) - iN_{n+1/2}(x) \right]
   * @f]
   * where @f$ J_\nu(x) @f$ and @f$ N_\nu(x) @f$ are the cylindrical Bessel
   * and Neumann functions respectively (@see cyl_bessel and cyl_neumann).
   *
   * @tparam _Tp The real type of the argument
   * @param __nu The real order
   * @param __z The real argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tpnu, _Tp>>
    cyl_hankel_2(_Tpnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_2<__type>(__nu, __z);
    }

  // Spherical Hankel functions of the first kind

  /**
   * Return the spherical Hankel function of the first kind @f$ h^{(1)}_n(x) @f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_1 for details.
   */
  inline std::complex<float>
  sph_hankel_1f(unsigned int __n, float __z)
  { return std::__detail::__sph_hankel_1<float>(__n, __z); }

  /**
   * Return the spherical Hankel function of the first kind @f$ h^{(1)}_n(x) @f$
   * of nonnegative order n and @c <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_1 for details.
   */
  inline std::complex<long double>
  sph_hankel_1l(unsigned int __n, long double __z)
  { return std::__detail::__sph_hankel_1<long double>(__n, __z); }

  /**
   * Return the spherical Hankel function of the first kind @f$ h^{(1)}_n(x) @f$
   * of nonnegative order @f$ n @f$ and real argument @f$ x >= 0 @f$.
   *
   * The spherical Hankel function of the first kind is defined by:
   * @f[
   *    h^{(1)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(1)}_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative order
   * @param __z The real argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    sph_hankel_1(unsigned int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sph_hankel_1<__type>(__n, __z);
    }

  // Spherical Hankel functions of the second kind

  /**
   * Return the spherical Hankel function of the second kind @f$ h^{(2)}_n(x)@f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_2 for details.
   */
  inline std::complex<float>
  sph_hankel_2f(unsigned int __n, float __z)
  { return std::__detail::__sph_hankel_2<float>(__n, __z); }

  /**
   * Return the spherical Hankel function of the second kind @f$ h^{(2)}_n(x)@f$
   * of nonnegative order n and @c <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_2 for details.
   */
  inline std::complex<long double>
  sph_hankel_2l(unsigned int __n, long double __z)
  { return std::__detail::__sph_hankel_2<long double>(__n, __z); }

  /**
   * Return the spherical Hankel function of the second kind @f$ h^{(2)}_n(x)@f$
   * of nonnegative order @f$ n @f$ and real argument @f$ x >= 0 @f$.
   *
   * The spherical Hankel function of the second kind is defined by:
   * @f[
   *    h^{(2)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(2)}_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative order
   * @param __z The real argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    sph_hankel_2(unsigned int __n, _Tp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sph_hankel_2<__type>(__n, __z);
    }

  // Modified spherical Bessel functions of the first kind

  /**
   * Return the regular modified spherical Bessel function @f$ i_n(x) @f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_i for details.
   */
  inline float
  sph_bessel_if(unsigned int __n, float __x)
  {
    float __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<float>(__n, __x,
        				  __i_n, __k_n, __ip_n, __kp_n);
    return __i_n;
  }

  /**
   * Return the regular modified spherical Bessel function @f$ i_n(x) @f$
   * of nonnegative order n and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_i for details.
   */
  inline long double
  sph_bessel_il(unsigned int __n, long double __x)
  {
    long double __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<long double>(__n, __x,
        					__i_n, __k_n, __ip_n, __kp_n);
    return __i_n;
  }

  /**
   * Return the regular modified spherical Bessel function @f$ i_n(x) @f$
   * of nonnegative order n and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  i_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} I_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __n  The integral order <tt> n >= 0 </tt>
   * @param  __x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sph_bessel_i(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      __type __i_n, __k_n, __ip_n, __kp_n;
      std::__detail::__sph_bessel_ik<__type>(__n, __x,
        				     __i_n, __k_n, __ip_n, __kp_n);
      return __i_n;
    }

  // Modified spherical Bessel functions of the second kind

  /**
   * Return the irregular modified spherical Bessel function @f$ k_n(x) @f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_k for more details.
   */
  inline float
  sph_bessel_kf(unsigned int __n, float __x)
  {
    float __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<float>(__n, __x,
        				  __i_n, __k_n, __ip_n, __kp_n);
    return __k_n;
  }

  /**
   * Return the irregular modified spherical Bessel function @f$ k_n(x) @f$
   * of nonnegative order n and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_k for more details.
   */
  inline long double
  sph_bessel_kl(unsigned int __n, long double __x)
  {
    long double __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<long double>(__n, __x,
        					__i_n, __k_n, __ip_n, __kp_n);
    return __k_n;
  }

  /**
   * Return the irregular modified spherical Bessel function @f$ k_n(x) @f$
   * of nonnegative order n and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  k_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} K_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __n  The integral order <tt> n >= 0 </tt>
   * @param  __x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    sph_bessel_k(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      __type __i_n, __k_n, __ip_n, __kp_n;
      std::__detail::__sph_bessel_ik<__type>(__n, __x,
        				     __i_n, __k_n, __ip_n, __kp_n);
      return __k_n;
    }

  // Airy functions of the first kind

  /**
   * Return the Airy function @f$ Ai(x) @f$ for @c float argument @f$ x @f$.
   *
   * @see airy_ai for details.
   */
  inline float
  airy_aif(float __x)
  {
    float __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<float>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Ai;
  }

  /**
   * Return the Airy function @f$ Ai(x) @f$ for <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see airy_ai for details.
   */
  inline long double
  airy_ail(long double __x)
  {
    long double __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<long double>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Ai;
  }

  /**
   * Return the Airy function @f$ Ai(x) @f$ of real argument @f$ x @f$.
   *
   * The Airy function is defined by:
   * @f[
   *    Ai(x) = \frac{1}{\pi}\int_0^\infty
   *      \cos \left(\frac{t^3}{3} + xt \right)dt
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    airy_ai(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      __type __Ai, __Bi, __Aip, __Bip;
      std::__detail::__airy<__type>(__x, __Ai, __Bi, __Aip, __Bip);
      return __Ai;
    }

  /**
   * Return the Airy function @f$ Ai(x) @f$ of complex argument @f$ x @f$.
   *
   * The Airy function is defined by:
   * @f[
   *    Ai(x) = \frac{1}{\pi}\int_0^\infty
   *      \cos \left(\frac{t^3}{3} + xt \right)dt
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    airy_ai(std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__airy_ai<__type>(__x);
    }

  // Airy functions of the second kind

  /**
   * Return the Airy function @f$ Bi(x) @f$ for @c float argument @f$ x @f$.
   *
   * @see airy_bi for details.
   */
  inline float
  airy_bif(float __x)
  {
    float __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<float>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Bi;
  }

  /**
   * Return the Airy function @f$ Bi(x) @f$ for <tt>long double</tt>
   * argument @f$ x @f$.
   *
   * @see airy_bi for details.
   */
  inline long double
  airy_bil(long double __x)
  {
    long double __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<long double>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Bi;
  }

  /**
   * Return the Airy function @f$ Bi(x) @f$ of real argument @f$ x @f$.
   *
   * The Airy function is defined by:
   * @f[
   *    Bi(x) = \frac{1}{\pi}\int_0^\infty \left[
   *           \exp \left(-\frac{t^3}{3} + xt \right)
   *          + \sin \left(\frac{t^3}{3} + xt \right) \right] dt
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    airy_bi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      __type __Ai, __Bi, __Aip, __Bip;
      std::__detail::__airy<__type>(__x, __Ai, __Bi, __Aip, __Bip);
      return __Bi;
    }

  /**
   * Return the Airy function @f$ Bi(x) @f$ of complex argument @f$ x @f$.
   *
   * The Airy function is defined by:
   * @f[
   *    Bi(x) = \frac{1}{\pi}\int_0^\infty \left[ 
   *           \exp \left(-\frac{t^3}{3} + xt \right)
   *          + \sin \left(\frac{t^3}{3} + xt \right) \right] dt
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __x The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    airy_bi(std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__airy_bi<__type>(__x);
    }

  // Log Gamma function for complex argument.

  /**
   * Return the logarithm of the gamma function for
   * <tt> std::complex<float> </tt> argument.
   *
   * @see lgamma for details.
   */
  inline std::complex<float>
  lgammaf(std::complex<float> __a)
  { return std::__detail::__log_gamma<std::complex<float>>(__a); }

  /**
   * Return the logarithm of the gamma function for
   * <tt> std::complex<long double> </tt> argument.
   *
   * @see lgamma for details.
   */
  inline std::complex<long double>
  lgammal(std::complex<long double> __a)
  { return std::__detail::__log_gamma<std::complex<long double>>(__a); }

  /**
   * Return the logarithm of the gamma function for complex argument.
   */
  template<typename _Ta>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Ta>>
    lgamma(std::complex<_Ta> __a)
    {
      using __type = std::complex<__gnu_cxx::__promote_fp_t<_Ta>>;
      return std::__detail::__log_gamma<__type>(__a);
    }

  // Gamma function for complex argument.

  /**
   * Return the gamma function for <tt> std::complex<float> </tt> argument.
   *
   * @see lgamma for details.
   */
  inline std::complex<float>
  tgammaf(std::complex<float> __a)
  { return std::__detail::__gamma<std::complex<float>>(__a); }

  /**
   * Return the gamma function for <tt> std::complex<long double> </tt>
   * argument.
   *
   * @see lgamma for details.
   */
  inline std::complex<long double>
  tgammal(std::complex<long double> __a)
  { return std::__detail::__gamma<std::complex<long double>>(__a); }

  /**
   * Return the gamma function for complex argument.
   */
  template<typename _Ta>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Ta>>
    tgamma(std::complex<_Ta> __a)
    {
      using __type = std::complex<__gnu_cxx::__promote_fp_t<_Ta>>;
      return std::__detail::__gamma<__type>(__a);
    }

  // Upper incomplete gamma functions

  /**
   * Return the upper incomplete gamma function @f$ \Gamma(a,x) @f$
   * for @c float argument.
   *
   * @see tgamma for details.
   */
  inline float
  tgammaf(float __a, float __x)
  { return std::__detail::__tgamma<float>(__a, __x); }

  /**
   * Return the upper incomplete gamma function @f$ \Gamma(a,x) @f$
   * for <tt>long double</tt> argument.
   *
   * @see tgamma for details.
   */
  inline long double
  tgammal(long double __a, long double __x)
  { return std::__detail::__tgamma<long double>(__a, __x); }

  /**
   * Return the upper incomplete gamma function @f$ \Gamma(a,x) @f$.
   * The (upper) incomplete gamma function is defined by
   * @f[
   *   \Gamma(a,x) = \int_x^\infty t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Ta, _Tp>
    tgamma(_Ta __a, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ta, _Tp>;
      return std::__detail::__tgamma<__type>(__a, __x);
    }

  // Lower incomplete gamma functions

  /**
   * Return the lower incomplete gamma function @f$ \gamma(a,x) @f$
   * for @c float argument.
   *
   * @see tgamma_lower for details.
   */
  inline float
  tgamma_lowerf(float __a, float __x)
  { return std::__detail::__tgamma_lower<float>(__a, __x); }

  /**
   * Return the lower incomplete gamma function @f$ \gamma(a,x) @f$
   * for <tt>long double</tt> argument.
   *
   * @see tgamma_lower for details.
   */
  inline long double
  tgamma_lowerl(long double __a, long double __x)
  { return std::__detail::__tgamma_lower<long double>(__a, __x); }

  /**
   * Return the lower incomplete gamma function @f$ \gamma(a,x) @f$.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \gamma(a,x) = \int_0^x t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Ta, _Tp>
    tgamma_lower(_Ta __a, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ta, _Tp>;
      return std::__detail::__tgamma_lower<__type>(__a, __x);
    }

  // Dilogarithm functions

  /**
   * Return the dilogarithm function @f$ \psi(z) @f$
   * for @c float argument.
   *
   * @see dilog for details.
   */
  inline float
  dilogf(float __x)
  { return std::__detail::__dilog<float>(__x); }

  /**
   * Return the dilogarithm function @f$ \psi(z) @f$
   * for <tt>long double</tt> argument.
   *
   * @see dilog for details.
   */
  inline long double
  dilogl(long double __x)
  { return std::__detail::__dilog<long double>(__x); }

  /**
   * Return the dilogarithm function @f$ \psi(z) @f$
   * for real argument.
   *
   * The dilogarithm is defined by:
   * @f[
   *    Li_2(x) = \sum_{k=1}^{\infty}\frac{x^k}{k^2}
   * @f]
   *
   * @param __x The argument.
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    dilog(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__dilog<__type>(__x);
    }

  // Complete Carlson elliptic R_F functions

  /**
   * Return the complete Carlson elliptic function @f$ R_F(x,y,z) @f$
   * for @c float arguments.
   *
   * @see comp_ellint_rf for details.
   */
  inline float
  comp_ellint_rf(float __x, float __y)
  { return std::__detail::__comp_ellint_rf<float>(__x, __y); }

  /**
   * Return the complete Carlson elliptic function @f$ R_F(x,y) @f$
   * for <tt>long double</tt> arguments.
   *
   * @see comp_ellint_rf for details.
   */
  inline long double
  comp_ellint_rf(long double __x, long double __y)
  { return std::__detail::__comp_ellint_rf<long double>(__x, __y); }

  /**
   * Return the complete Carlson elliptic function @f$ R_F(x,y) @f$
   * for real arguments.
   *
   * The complete Carlson elliptic function of the first kind is defined by:
   * @f[
   *    R_F(x,y) = R_F(x,y,y) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)}
   * @f]
   *
   * @param  __x  The first argument.
   * @param  __y  The second argument.
   */
  template<typename _Tx, typename _Ty>
    inline __gnu_cxx::__promote_fp_t<_Tx, _Ty>
    comp_ellint_rf(_Tx __x, _Ty __y)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tx, _Ty>;
      return std::__detail::__comp_ellint_rf<__type>(__x, __y);
    }

  // Carlson elliptic R_F functions

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * of the first kind for @c float arguments.
   *
   * @see ellint_rf for details.
   */
  inline float
  ellint_rff(float __x, float __y, float __z)
  { return std::__detail::__ellint_rf<float>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * of the first kind for <tt>long double</tt> arguments.
   *
   * @see ellint_rf for details.
   */
  inline long double
  ellint_rfl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rf<long double>(__x, __y, __z); }

  /**
   * Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * of the first kind for real arguments.
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
    inline __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp>
    ellint_rf(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp, _Up>
    ellint_rc(_Tp __x, _Up __y)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp, _Wp>
    ellint_rj(_Tp __x, _Up __y, _Vp __z, _Wp __p)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp, _Wp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp>
    ellint_rd(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp>;
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
   * Return the complete Carlson elliptic function @f$ R_G(x,y) @f$
   * for real arguments.
   *
   * The complete Carlson elliptic function is defined by:
   * @f[
   *    R_G(x,y) = R_G(x,y,y) = \frac{1}{4} \int_0^\infty
   *            dt t (t + x)^{-1/2}(t + y)^{-1}
   *         (\frac{x}{t + x} + \frac{2y}{t + y})
   * @f]
   *
   * @param  __x  The first argument.
   * @param  __y  The second argument.
   */
  template<typename _Tx, typename _Ty>
    inline __gnu_cxx::__promote_fp_t<_Tx, _Ty>
    comp_ellint_rg(_Tx __x, _Ty __y)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tx, _Ty>;
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
   * @f]
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
    inline __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp>
    ellint_rg(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up, _Vp>;
      return std::__detail::__ellint_rg<__type>(__x, __y, __z);
    }

  // Hurwitz zeta functions

  /**
   * Return the Hurwitz zeta function of @c float argument @f$ s @f$,
   * and parameter @f$ a @f$.
   *
   * @see hurwitz_zeta for details.
   */
  inline float
  hurwitz_zetaf(float __s, float __a)
  { return std::__detail::__hurwitz_zeta<float>(__s, __a); }

  /**
   * Return the Hurwitz zeta function of <tt>long double</tt> argument @f$ s @f$,
   * and parameter @f$ a @f$.
   *
   * @see hurwitz_zeta for details.
   */
  inline long double
  hurwitz_zetal(long double __s, long double __a)
  { return std::__detail::__hurwitz_zeta<long double>(__s, __a); }

  /**
   * Return the Hurwitz zeta function of real argument @f$ s @f$, and parameter @f$ a @f$.
   *
   * The the Hurwitz zeta function is defined by
   * @f[
   *    \zeta(s, a) = \sum_{n=0}^{\infty}\frac{1}{(a + n)^s}
   * @f]
   *
   * @param __s The argument
   * @param __a The parameter
   */
  template<typename _Tp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Up>
    hurwitz_zeta(_Tp __s, _Up __a)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up>;
      return std::__detail::__hurwitz_zeta<__type>(__s, __a);
    }

  /**
   * Return the Hurwitz zeta function of real argument @f$ s @f$,
   * and complex parameter @f$ a @f$.
   *
   * @see hurwitz_zeta for details.
   */
  template<typename _Tp, typename _Up>
    std::complex<_Tp>
    hurwitz_zeta(_Tp __s, std::complex<_Up> __a)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Up>;
      return std::__detail::__hurwitz_zeta<__type>(__s, __a);
    }

  // Digamma or psi functions

  /**
   * Return the psi or digamma function of @c float argument @f$ x @f$.
   *
   * @see psi for details.
   */
  inline float
  psif(float __x)
  { return std::__detail::__psi<float>(__x); }

  /**
   * Return the psi or digamma function of <tt>long double</tt> argument
   * @f$ x @f$.
   *
   * @see psi for details.
   */
  inline long double
  psil(long double __x)
  { return std::__detail::__psi<long double>(__x); }

  /**
   * Return the psi or digamma function of argument @f$ x @f$.
   *
   * The the psi or digamma function is defined by
   * @f[
   *    \psi(x) = \frac{d}{dx}log\left(\Gamma(x)\right)
   *            = \frac{\Gamma'(x)}{\Gamma(x)}
   * @f]
   *
   * @param __x The parameter
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    psi(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__psi<__type>(__x);
    }

  // Incomplete beta functions

  /**
   * Return the regularized incomplete beta function of parameters @f$ a @f$, @f$ b @f$,
   * and argument @f$ x @f$.
   *
   * See ibeta for details.
   */
  inline float
  ibetaf(float __a, float __b, float __x)
  { return std::__detail::__beta_inc<float>(__a, __b, __x); }

  /**
   * Return the regularized incomplete beta function of parameters @f$ a @f$, @f$ b @f$,
   * and argument @f$ x @f$.
   *
   * See ibeta for details.
   */
  inline long double
  ibetal(long double __a, long double __b, long double __x)
  { return std::__detail::__beta_inc<long double>(__a, __b, __x); }

  /**
   * Return the regularized incomplete beta function of parameters @f$ a @f$, @f$ b @f$,
   * and argument @f$ x @f$.
   *
   * The regularized incomplete beta function is defined by
   * @f[
   *    I_x(a, b) = \frac{B_x(a,b)}{B(a,b)}
   * @f]
   * where
   * @f[
   *   B_x(a,b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt
   * @f]
   * is the non-regularized incomplete beta function and @f$ B(a,b) @f$
   * is the usual beta function.
   *
   * @param __a The first parameter
   * @param __b The second parameter
   * @param __x The argument
   */
  template<typename _Ta, typename _Tb, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Ta, _Tb, _Tp>
    ibeta(_Ta __a, _Tb __b, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ta, _Tb, _Tp>;
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
   * of parameters @f$ a @f$, @f$ b @f$, and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Ta, _Tb, _Tp>
    ibetac(_Ta __a, _Tb __b, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ta, _Tb, _Tp>;
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
   * Return the Fresnel sine integral of argument @f$ x @f$.
   *
   * The Fresnel sine integral is defined by
   * @f[
   *    S(x) = \int_0^x \sin(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    fresnel_s(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
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
   * Return the Fresnel cosine integral of argument @f$ x @f$.
   *
   * The Fresnel cosine integral is defined by
   * @f[
   *    C(x) = \int_0^x \cos(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param __x The argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    fresnel_c(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::real(std::__detail::__fresnel<__type>(__x));
    }

  // Dawson integral

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for @c float argument @f$ x @f$.
   *
   * @see dawson for details.
   */
  inline float
  dawsonf(float __x)
  { return std::__detail::__dawson<float>(__x); }

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see dawson for details.
   */
  inline long double
  dawsonl(long double __x)
  { return std::__detail::__dawson<long double>(__x); }

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for real argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    dawson(_Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__dawson<__type>(__x);
    }

  //  Exponential integrals

  /**
   * Return the exponential integral @f$ E_n(x) @f$ for integral
   * order @f$ n @f$ and @c float argument @f$ x @f$.
   *
   * @see expint for details.
   */
  inline float
  expintf(unsigned int __n, float __x)
  { return std::__detail::__expint<float>(__n, __x); }

  /**
   * Return the exponential integral @f$ E_n(x) @f$ for integral
   * order @f$ n @f$ and <tt>long double</tt> argument @f$ x @f$.
   *
   * @see expint for details.
   */
  inline long double
  expintl(unsigned int __n, long double __x)
  { return std::__detail::__expint<long double>(__n, __x); }

  /**
   * Return the exponential integral @f$ E_n(x) @f$ of integral
   * order @f$ n @f$ and real argument @f$ x @f$.
   * The exponential integral is defined by:
   * @f[
   *    E_n(x) = \int_1^\infty \frac{e^{-tx}}{t^n}dt
   * @f]
   * In particular
   * @f[
   *    E_1(x) = \int_1^\infty \frac{e^{-tx}}{t}dt = -Ei(-x)
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __n The integral order
   * @param __x The real argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    expint(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__expint<__type>(__n, __x);
    }

  //  Log upper Pochhammer symbol

  inline float
  lpochhammerf(float __a, float __n)
  { return std::__detail::__log_pochhammer<float>(__a, __n); }

  inline long double
  lpochhammerl(long double __a, long double __n)
  { return std::__detail::__log_pochhammer<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Tn>
    lpochhammer(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Tn>;
      return std::__detail::__log_pochhammer<__type>(__a, __n);
    }

  //  Log lower Pochhammer symbol

  inline float
  lpochhammer_lowerf(float __a, float __n)
  { return std::__detail::__log_pochhammer_lower<float>(__a, __n); }

  inline long double
  lpochhammer_lowerl(long double __a, long double __n)
  { return std::__detail::__log_pochhammer_lower<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Tn>
    lpochhammer_lower(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Tn>;
      return std::__detail::__log_pochhammer_lower<__type>(__a, __n);
    }

  //  Upper Pochhammer symbols

  inline float
  pochhammerf(float __a, float __n)
  { return std::__detail::__pochhammer<float>(__a, __n); }

  inline long double
  pochhammerl(long double __a, long double __n)
  { return std::__detail::__pochhammer<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Tn>
    pochhammer(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Tn>;
      return std::__detail::__pochhammer<__type>(__a, __n);
    }

  //  Lower Pochhammer symbols

  inline float
  pochhammer_lowerf(float __a, float __n)
  { return std::__detail::__pochhammer_lower<float>(__a, __n); }

  inline long double
  pochhammer_lowerl(long double __a, long double __n)
  { return std::__detail::__pochhammer_lower<long double>(__a, __n); }

  /**
   * 
   */
  template<typename _Tp, typename _Tn>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Tn>
    pochhammer_lower(_Tp __a, _Tn __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Tn>;
      return std::__detail::__pochhammer_lower<__type>(__a, __n);
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    factorial(unsigned int __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    double_factorial(int __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    lfactorial(unsigned int __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    ldouble_factorial(int __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    bincoef(unsigned int __n, unsigned int __k)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tp>
    lbincoef(unsigned int __n, unsigned int __k)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__log_bincoef<__type>(__n, __k);
    }

  // Bernoulli numbers

  /**
   * Return the Bernoulli number of integer order @f$ n @f$ as a @c float.
   *
   * @see bernoulli for details.
   */
  inline float
  bernoullif(unsigned int __n)
  { return std::__detail::__bernoulli<float>(__n); }

  /**
   * Return the Bernoulli number of integer order @f$ n @f$ as a
   * <tt>long double</tt>.
   *
   * @see bernoulli for details.
   */
  inline long double
  bernoullil(unsigned int __n)
  { return std::__detail::__bernoulli<long double>(__n); }

  /**
   * Return the Bernoulli number of integer order @f$ n @f$.
   *
   * The Bernoulli numbers are defined by
   * @f[
   *    
   * @f]
   *
   * @param __n The order.
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    bernoulli(unsigned int __n)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__bernoulli<__type>(__n);
    }

  // Legendre functions of the second kind

  /**
   * Return the Legendre function of the second kind @f$ Q_l(x) @f$
   * of nonnegative degree @f$ l @f$ and @c float argument.
   *
   * @see legendre_q for details.
   */
  inline float
  legendre_qf(unsigned int __n, float __x)
  { return std::__detail::__legendre_q<float>(__n, __x); }

  /**
   * Return the Legendre function of the second kind @f$ Q_l(x) @f$
   * of nonnegative degree @f$ l @f$ and <tt>long double</tt> argument.
   *
   * @see legendre_q for details.
   */
  inline long double
  legendre_ql(unsigned int __n, long double __x)
  { return std::__detail::__legendre_q<long double>(__n, __x); }

  /**
   * Return the Legendre function of the second kind @f$ Q_l(x) @f$ of
   * nonnegative degree @f$ l @f$ and real argument @f$ |x| <= 0 @f$.
   *
   * The Legendre function of the second kind of order @f$ l @f$
   * and argument @f$ x @f$, @f$ Q_l(x) @f$, is defined by:
   * @f[
   *   Q_l(x) = \frac{1}{2} \log{\frac{x+1}{x-1}} P_l(x)
   *           - \sum_{k=0}^{l-1}\frac{(l+k)!}{(l-k)!(k!)^2s^k}
   *             \left[\psi(l+1) - \psi(k+1)\right](x-1)^k
   * @f]
   * where @f$ P_l(x) @f$ is the Legendre polynomial of degree @f$ l @f$
   * and @f$ \psi(x) @f$ is the psi or dilogarithm function.
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __l The degree @f$ l >= 0 @f$
   * @param __x The argument @c abs(__x) <= 1
   * @throw std::domain_error if @c abs(__x) > 1
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    legendre_q(unsigned int __n, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__legendre_q<__type>(__n, __x);
    }

  // Scaled lower incomplete gamma

  inline float
  pgammaf(float __a, float __x)
  { return std::__detail::__pgamma<float>(__a, __x); }

  inline long double
  pgammal(long double __a, long double __x)
  { return std::__detail::__pgamma<long double>(__a, __x); }

  /**
   * 
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Ta, _Tp>
    pgamma(_Ta __a, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ta, _Tp>;
      return std::__detail::__pgamma<__type>(__a, __x);
    }

  // Scaled upper incomplete gamma

  inline float
  qgammaf(float __a, float __x)
  { return std::__detail::__qgamma<float>(__a, __x); }

  inline long double
  qgammal(long double __a, long double __x)
  { return std::__detail::__qgamma<long double>(__a, __x); }

  /**
   * 
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Ta, _Tp>
    qgamma(_Ta __a,_Tp  __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ta, _Tp>;
      return std::__detail::__qgamma<__type>(__a, __x);
    }

  // Jacobi zeta functions.

  inline float
  jacobi_zetaf(float __k, float __phi)
  { return std::__detail::__jacobi_zeta<float>(__k, __phi); }

  inline long double
  jacobi_zetal(long double __k, long double __phi)
  { return std::__detail::__jacobi_zeta<long double>(__k, __phi); }

  /**
   * Return the Jacobi zeta function of @f$ k @f$ and @f$ \phi @f$.
   *
   * The Jacobi zeta function is defined by
   * @f[
   *    Z(m,\phi) = E(m,\phi) - \frac{E(m)F(m,\phi)}{K(m)}
   * @f]
   * where @f$ E(m,\phi) @f$ is the elliptic function of the second kind,
   * @f$ E(m) @f$ is the complete ellitic function of the second kind,
   * and @f$ F(m,\phi) @f$ is the elliptic function of the first kind.
   *
   * @tparam _Tk the real type of the modulus
   * @tparam _Tphi the real type of the angle limit
   * @param __k The modulus
   * @param __phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline __gnu_cxx::__promote_fp_t<_Tk, _Tphi>
    jacobi_zeta(_Tk __k, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tphi>;
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
   * Return the Heuman lambda function @f$ \Lambda(k,\phi) @f$
   * of modulus @f$ k @f$ and angular limit @f$ \phi @f$.
   *
   * The complete Heuman lambda function is defined by
   * @f[
   *    \Lambda(k,\phi) = \frac{F(1-m,\phi)}{K(1-m)}
   *    + \frac{2}{\pi} K(m) Z(1-m,\phi)
   * @f]
   * where @f$ m = k^2 @f$, @f$ K(k) @f$ is the complete elliptic function
   * of the first kind, and @f$ Z(k,phi) @f$ is the Jacobi zeta function.
   *
   * @tparam _Tk the floating-point type of the modulus
   * @tparam _Tphi the floating-point type of the angular limit argument
   * @param __k The modulus
   * @param __phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline __gnu_cxx::__promote_fp_t<_Tk, _Tphi>
    heuman_lambda(_Tk __k, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tphi>;
      return std::__detail::__heuman_lambda<__type>(__k, __phi);
    }

  // Complete Legendre elliptic integral D.

  /**
   * Return the complete Legendre elliptic integral @f$ D(k) @f$
   * of @c float modulus @f$ k @f$.
   *
   * @see comp_ellint_d for details.
   */
  inline float
  comp_ellint_df(float __k)
  { return std::__detail::__comp_ellint_d<float>(__k); }

  /**
   * Return the complete Legendre elliptic integral @f$ D(k) @f$
   * of <tt>long double</tt> modulus @f$ k @f$.
   *
   * @see comp_ellint_d for details.
   */
  inline long double
  comp_ellint_dl(long double __k)
  { return std::__detail::__comp_ellint_d<long double>(__k); }

  /**
   * Return the complete Legendre elliptic integral @f$ D(k) @f$
   * of real modulus @f$ k @f$.
   *
   * The complete Legendre elliptic integral D is defined by
   * @f[
   *    D(k) = \int_0^{\pi/2} \frac{\sin^2\theta d\theta}{\sqrt{1-k^2sin2\theta}}
   * @f]
   *
   * @tparam _Tk The type of the modulus @c k
   * @param __k The modulus <tt>-1 <= __k <= +1</tt>
   */
  template<typename _Tk>
    inline __gnu_cxx::__promote_fp_t<_Tk>
    comp_ellint_d(_Tk __k)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk>;
      return std::__detail::__comp_ellint_d<__type>(__k);
    }

  // Legendre elliptic integrals D.

  /**
   * Return the incomplete Legendre elliptic integral @f$ D(k, \phi) @f$
   * of @c float modulus @f$ k @f$ and angular limit @f$ \phi @f$.
   *
   * @see ellint_d for details.
   */
  inline float
  ellint_df(float __k, float __phi)
  { return std::__detail::__ellint_d<float>(__k, __phi); }

  /**
   * Return the incomplete Legendre elliptic integral @f$ D(k, \phi) @f$
   * of <tt>long double</tt> modulus @f$ k @f$ and angular limit @f$ \phi @f$.
   *
   * @see ellint_d for details.
   */
  inline long double
  ellint_dl(long double __k, long double __phi)
  { return std::__detail::__ellint_d<long double>(__k, __phi); }

  /**
   * Return the incomplete Legendre elliptic integral @f$ D(k,\phi) @f$
   * of real modulus @f$ k @f$ and angular limit @f$ \phi @f$.
   *
   * The Legendre elliptic integral D is defined by
   * @f[
   *    D(k,\phi) = \int_0^\phi
   *         \frac{\sin^2\theta d\theta}{\sqrt{1-k^2sin^2\theta}}
   * @f]
   *
   * @param __k The modulus <tt>-1 <= __k <= +1</tt>
   * @param __phi The angle
   */
  template<typename _Tk, typename _Tphi>
    inline __gnu_cxx::__promote_fp_t<_Tk, _Tphi>
    ellint_d(_Tk __k, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tphi>;
      return std::__detail::__ellint_d<__type>(__k, __phi);
    }

  // Bulirsch elliptic integrals of the first kind.

  /**
   * Return the Bulirsch elliptic integral @f$ el1(x,k_c) @f$
   * of the first kind of @c float tangent limit @f$ x @f$
   * and complementary modulus @f$ k_c @f$.
   *
   * @see ellint_el1 for details.
   */
  inline float
  ellint_el1f(float __x, float __k_c)
  { return std::__detail::__ellint_el1<float>(__x, __k_c); }

  /**
   * Return the Bulirsch elliptic integral @f$ el1(x,k_c) @f$
   * of the first kind of real tangent limit @f$ x @f$
   * and complementary modulus @f$ k_c @f$.
   *
   * @see ellint_el1 for details.
   */
  inline long double
  ellint_el1l(long double __x, long double __k_c)
  { return std::__detail::__ellint_el1<long double>(__x, __k_c); }

  /**
   * Return the Bulirsch elliptic integral @f$ el1(x,k_c) @f$
   * of the first kind of real tangent limit @f$ x @f$
   * and complementary modulus @f$ k_c @f$.
   *
   * The Bulirsch elliptic integral of the first kind is defined by
   * @f[
   *    el1(x,k_c) = el2(x,k_c,1,1) = \int_0^{\arctan x} \frac{1+1\tan^2\theta}
   *           {\sqrt{(1+\tan^2\theta)(1+k_c^2\tan^2\theta)}}d\theta
   * @f]
   *
   * @param __x The tangent of the angular integration limit
   * @param __k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   */
  template<typename _Tp, typename _Tk>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Tk>
    ellint_el1(_Tp __x, _Tk __k_c)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Tk>;
      return std::__detail::__ellint_el1<__type>(__x, __k_c);
    }

  // Bulirsch elliptic integrals of the second kind.

  /**
   * Return the Bulirsch elliptic integral of the second kind
   * @f$ el2(x,k_c,a,b) @f$.
   *
   * @see ellint_el2 for details.
   */
  inline float
  ellint_el2f(float __x, float __k_c, float __a, float __b)
  { return std::__detail::__ellint_el2<float>(__x, __k_c, __a, __b); }

  /**
   * Return the Bulirsch elliptic integral of the second kind
   * @f$ el2(x,k_c,a,b) @f$.
   *
   * @see ellint_el2 for details.
   */
  inline long double
  ellint_el2l(long double __x, long double __k_c,
	      long double __a, long double __b)
  { return std::__detail::__ellint_el2<long double>(__x, __k_c, __a, __b); }

  /**
   * Return the Bulirsch elliptic integral of the second kind
   * @f$ el2(x,k_c,a,b) @f$.
   *
   * The Bulirsch elliptic integral of the second kind is defined by
   * @f[
   *    el2(x,k_c,a,b) = \int_0^{\arctan x} \frac{a+b\tan^2\theta}
   *           {\sqrt{(1+\tan^2\theta)(1+k_c^2\tan^2\theta)}}d\theta
   * @f]
   *
   * @param __x The tangent of the angular integration limit
   * @param __k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   * @param __a The  parameter
   * @param __b The  parameter
   */
  template<typename _Tp, typename _Tk, typename _Ta, typename _Tb>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Tk, _Ta, _Tb>
    ellint_el2(_Tp __x, _Tk __k_c, _Ta __a, _Tb __b)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Tk, _Ta, _Tb>;
      return std::__detail::__ellint_el2<__type>(__x, __k_c, __a, __b);
    }

  // Bulirsch elliptic integrals of the third kind.

  /**
   * Return the Bulirsch elliptic integral of the third kind
   * @f$ el3(x,k_c,p) @f$ of @c float tangent limit @f$ x @f$,
   * complementary modulus @f$ k_c @f$, and parameter @f$ p @f$.
   *
   * @see ellint_el3 for details.
   */
  inline float
  ellint_el3f(float __x, float __k_c, float __p)
  { return std::__detail::__ellint_el3<float>(__x, __k_c, __p); }

  /**
   * Return the Bulirsch elliptic integral of the third kind
   * @f$ el3(x,k_c,p) @f$ of <tt>long double</tt> tangent limit @f$ x @f$,
   * complementary modulus @f$ k_c @f$, and parameter @f$ p @f$.
   *
   * @see ellint_el3 for details.
   */
  inline long double
  ellint_el3l(long double __x, long double __k_c, long double __p)
  { return std::__detail::__ellint_el3<long double>(__x, __k_c, __p); }

  /**
   * Return the Bulirsch elliptic integral of the third kind
   * @f$ el3(x,k_c,p) @f$ of real tangent limit @f$ x @f$,
   * complementary modulus @f$ k_c @f$, and parameter @f$ p @f$.
   *
   * The Bulirsch elliptic integral of the third kind is defined by
   * @f[
   *    el3(x,k_c,p) = \int_0^{\arctan x} \frac{d\theta}
   *          {(cos^2\theta+p\sin^2\theta)\sqrt{cos^2\theta+k_c^2\sin^2\theta}}
   * @f]
   *
   * @param __x The tangent of the angular integration limit
   * @param __k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   * @param __p The paramenter
   */
  template<typename _Tx, typename _Tk, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tx, _Tk, _Tp>
    ellint_el3(_Tx __x, _Tk __k_c, _Tp __p)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tx, _Tk, _Tp>;
      return std::__detail::__ellint_el3<__type>(__x, __k_c, __p);
    }

  // Bulirsch complete elliptic integrals.

  /**
   * Return the Bulirsch complete elliptic integral @f$ cel(k_c,p,a,b) @f$
   * of real complementary modulus @f$ k_c @f$, and parameters @f$ p @f$,
   * @f$ a @f$, and @f$ b @f$.
   *
   * @see ellint_cel for details.
   */
  inline float
  ellint_celf(float __k_c, float __p, float __a, float __b)
  { return std::__detail::__ellint_cel<float>(__k_c, __p, __a, __b); }

  /**
   * Return the Bulirsch complete elliptic integral @f$ cel(k_c,p,a,b) @f$.
   *
   * @see ellint_cel for details.
   */
  inline long double
  ellint_cell(long double __k_c, long double __p,
	      long double __a, long double __b)
  { return std::__detail::__ellint_cel<long double>(__k_c, __p, __a, __b); }

  /**
   * Return the Bulirsch complete elliptic integral @f$ cel(k_c,p,a,b) @f$
   * of real complementary modulus @f$ k_c @f$, and parameters @f$ p @f$,
   * @f$ a @f$, and @f$ b @f$.
   *
   * The Bulirsch complete elliptic integral is defined by
   * @f[
   *    cel(k_c,p,a,b)=\int_0^{\pi/2}
   *        \frac{a\cos^2\theta + b\sin^2\theta}{cos^2\theta + p\sin^2\theta}
   *        \frac{d\theta}{\sqrt{cos^2\theta + k_c^2\sin^2\theta}}
   * @f]
   *
   * @param __k_c The complementary modulus @f$ k_c = \sqrt{1 - k^2} @f$
   * @param __p The  parameter
   * @param __a The  parameter
   * @param __b The  parameter
   */
  template<typename _Tk, typename _Tp, typename _Ta, typename _Tb>
    inline __gnu_cxx::__promote_fp_t<_Tk, _Tp, _Ta, _Tb>
    ellint_cel(_Tk __k_c, _Tp __p, _Ta __a, _Tb __b)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tp, _Ta, _Tb>;
      return std::__detail::__ellint_cel<__type>(__k_c, __p, __a, __b);
    }

  // Cylindrical Hankel functions of the first kind.

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of <tt>std::complex<float></tt> order @f$ \nu @f$
   * and argument @f$ x @f$.
   *
   * @see cyl_hankel_1 for more details.
   */
  inline std::complex<float>
  cyl_hankel_1f(std::complex<float> __nu, std::complex<float> __x)
  { return std::__detail::__cyl_hankel_1<float>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of <tt>std::complex<long double></tt>
   * order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * @see cyl_hankel_1 for more details.
   */
  inline std::complex<long double>
  cyl_hankel_1l(std::complex<long double> __nu, std::complex<long double> __x)
  { return std::__detail::__cyl_hankel_1<long double>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of complex order @f$ \nu @f$
   * and argument @f$ x @f$.
   *
   * The cylindrical Hankel function of the first kind is defined by
   * @f[
   *    H^{(1)}_\nu(x) = J_\nu(x) + i N_\nu(x)
   * @f]
   *
   * @tparam _Tpnu The complex type of the order
   * @tparam _Tp The complex type of the argument
   * @param __nu The complex order
   * @param __x The complex argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tpnu, _Tp>>
    cyl_hankel_1(std::complex<_Tpnu> __nu, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_1<__type>(__nu, __x);
    }

  // Cylindrical Hankel functions of the second kind.

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of <tt>std::complex<float></tt> order @f$ \nu @f$
   * and argument @f$ x @f$.
   *
   * @see cyl_hankel_2 for more details.
   */
  inline std::complex<float>
  cyl_hankel_2f(std::complex<float> __nu, std::complex<float> __x)
  { return std::__detail::__cyl_hankel_2<float>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of <tt>std::complex<long double></tt>
   * order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * @see cyl_hankel_2 for more details.
   */
  inline std::complex<long double>
  cyl_hankel_2l(std::complex<long double> __nu, std::complex<long double> __x)
  { return std::__detail::__cyl_hankel_2<long double>(__nu, __x); }

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of complex order @f$ \nu @f$
   * and argument @f$ x @f$.
   *
   * The cylindrical Hankel function of the second kind is defined by
   * @f[
   *    H^{(2)}_\nu(x) = J_\nu(x) - i N_\nu(x)
   * @f]
   *
   * @tparam _Tpnu The complex type of the order
   * @tparam _Tp The complex type of the argument
   * @param __nu The complex order
   * @param __x The complex argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tpnu, _Tp>>
    cyl_hankel_2(std::complex<_Tpnu> __nu, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__cyl_hankel_2<__type>(__nu, __x);
    }

  // Spherical Hankel functions of the first kind.

  /**
   * Return the complex spherical Hankel function of the first kind
   * @f$ h^{(1)}_n(x) @f$ of non-negative integral @f$ n @f$
   * and <tt>std::complex<float></tt> argument @f$ x @f$.
   *
   * @see sph_hankel_1 for more details.
   */
  inline std::complex<float>
  sph_hankel_1f(unsigned int __n, std::complex<float> __x)
  { return std::__detail::__sph_hankel_1<float>(__n, __x); }

  /**
   * Return the complex spherical Hankel function of the first kind
   * @f$ h^{(1)}_n(x) @f$ of non-negative integral @f$ n @f$
   * and <tt>std::complex<long double></tt> argument @f$ x @f$.
   *
   * @see sph_hankel_1 for more details.
   */
  inline std::complex<long double>
  sph_hankel_1l(unsigned int __n, std::complex<long double> __x)
  { return std::__detail::__sph_hankel_1<long double>(__n, __x); }

  /**
   * Return the complex spherical Hankel function of the first kind
   * @f$ h^{(1)}_n(x) @f$ of non-negative integral @f$ n @f$
   * and complex argument @f$ x @f$.
   *
   * The spherical Hankel function of the first kind is defined by
   * @f[
   *    h^{(1)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(1)}_{n+1/2}(x)
   *                 = j_n(x) + i n_n(x)
   * @f]
   * where @f$ j_n(x) @f$ and @f$ n_n(x) @f$ are the spherical Bessel
   * and Neumann functions respectively.
   *
   * @param __n The integral order >= 0
   * @param __x The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    sph_hankel_1(unsigned int __n, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sph_hankel_1<__type>(__n, __x);
    }

  // Spherical Hankel functions of the second kind.

  /**
   * Return the complex spherical Hankel function of the second kind
   * @f$ h^{(2)}_n(x) @f$ of non-negative integral @f$ n @f$
   * and <tt>std::complex<float></tt> argument @f$ x @f$.
   *
   * @see sph_hankel_2 for more details.
   */
  inline std::complex<float>
  sph_hankel_2f(unsigned int __n, std::complex<float> __x)
  { return std::__detail::__sph_hankel_2<float>(__n, __x); }

  /**
   * Return the complex spherical Hankel function of the second kind
   * @f$ h^{(2)}_n(x) @f$ of non-negative integral @f$ n @f$
   * and <tt>std::complex<long double></tt> argument @f$ x @f$.
   *
   * @see sph_hankel_2 for more details.
   */
  inline std::complex<long double>
  sph_hankel_2l(unsigned int __n, std::complex<long double> __x)
  { return std::__detail::__sph_hankel_2<long double>(__n, __x); }

  /**
   * Return the complex spherical Hankel function of the second kind
   * @f$ h^{(2)}_n(x) @f$ of nonnegative order @f$ n @f$
   * and complex argument @f$ x @f$.
   *
   * The spherical Hankel function of the second kind is defined by
   * @f[
   *    h^{(2)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(2)}_{n+1/2}(x)
   *                 = j_n(x) - i n_n(x)
   * @f]
   * where @f$ j_n(x) @f$ and @f$ n_n(x) @f$ are the spherical Bessel
   * and Neumann functions respectively.
   *
   * @param __n The integral order >= 0
   * @param __x The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    sph_hankel_2(unsigned int __n, std::complex<_Tp> __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__sph_hankel_2<__type>(__n, __x);
    }

  // Spherical harmonic functions

  /**
   * Return the complex spherical harmonic function of degree @f$ l @f$, order @f$ m @f$,
   * and @c float zenith angle @f$ \theta @f$, and azimuth angle @f$ \phi @f$.
   *
   * @see sph_harmonic for details.
   */
  inline std::complex<float>
  sph_harmonicf(unsigned int __l, int __m,
		float __theta, float __phi)
  { return std::__detail::__sph_harmonic<float>(__l, __m, __theta, __phi); }

  /**
   * Return the complex spherical harmonic function of degree @f$ l @f$, order @f$ m @f$,
   * and <tt>long double</tt> zenith angle @f$ \theta @f$,
   * and azimuth angle @f$ \phi @f$.
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
   * Return the complex spherical harmonic function of degree @f$ l @f$, order @f$ m @f$,
   * and real zenith angle @f$ \theta @f$, and azimuth angle @f$ \phi @f$.
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
    inline std::complex<__gnu_cxx::__promote_fp_t<_Ttheta, _Tphi>>
    sph_harmonic(unsigned int __l, int __m, _Ttheta __theta, _Tphi __phi)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Ttheta, _Tphi>;
      return std::__detail::__sph_harmonic<__type>(__l, __m, __theta, __phi);
    }

  // Polylogarithm functions

  /**
   * Return the real polylogarithm function of real thing @c s
   * and real argument @f$ w @f$.
   *
   * @see polylog for details.
   */
  inline float
  polylogf(float __s, float __w)
  { return std::__detail::__polylog<float>(__s, __w); }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @f$ w @f$.
   *
   * @see polylog for details.
   */
  inline long double
  polylogl(long double __s, long double __w)
  { return std::__detail::__polylog<long double>(__s, __w); }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @f$ w @f$.
   *
   * The polylogarithm function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __s 
   * @param __w 
   */
  template<typename _Tp, typename _Wp>
    inline __gnu_cxx::__promote_fp_t<_Tp, _Wp>
    polylog(_Tp __s, _Wp __w)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Wp>;
      return std::__detail::__polylog<__type>(__s, __w);
    }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @f$ w @f$.
   *
   * @see polylog for details.
   */
  inline std::complex<float>
  polylogf(float __s, std::complex<float> __w)
  { return std::__detail::__polylog<float>(__s, __w); }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @f$ w @f$.
   *
   * @see polylog for details.
   */
  inline std::complex<long double>
  polylogl(long double __s, std::complex<long double> __w)
  { return std::__detail::__polylog<long double>(__s, __w); }

  /**
   * Return the complex polylogarithm function of real thing @c s
   * and complex argument @f$ w @f$.
   *
   * The polylogarithm function is defined by
   * @f[
   *    
   * @f]
   *
   * @param __s 
   * @param __w 
   */
  template<typename _Tp, typename _Wp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp, _Wp>>
    polylog(_Tp __s, std::complex<_Tp> __w)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp, _Wp>;
      return std::__detail::__polylog<__type>(__s, __w);
    }

  // Dirichlet eta function

  /**
   * Return the Dirichlet eta function of real argument @f$ s @f$.
   *
   * @see dirichlet_eta for details.
   */
  inline float
  dirichlet_etaf(float __s)
  { return std::__detail::__dirichlet_eta<float>(__s); }

  /**
   * Return the Dirichlet eta function of real argument @f$ s @f$.
   *
   * @see dirichlet_eta for details.
   */
  inline long double
  dirichlet_etal(long double __s)
  { return std::__detail::__dirichlet_eta<long double>(__s); }

  /**
   * Return the Dirichlet eta function of real argument @f$ s @f$.
   *
   * The Dirichlet eta function is defined by
   * @f[
   *    \eta(s) = \sum_{k=1}^\infty \frac{(-1)^k}{k^s}
   *    = \left( 1 - 2^{1-s} \right) \zeta(s)
   * @f]
   * An important reflection formula is:
   * @f[
   *    \eta(-s) = 2 \frac{1-2^{-s-1}}{1-2^{-s}} \pi^{-s-1} 
   *              s \sin(\frac{\pi s}{2}) \Gamma(s) \eta(s+1)
   * @f]
   *
   * @param __s 
   */
  template<typename _Tp>
    inline _Tp
    dirichlet_eta(_Tp __s)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__dirichlet_eta<__type>(__s);
    }

  // Dirichlet beta function

  /**
   * Return the Dirichlet beta function of real argument @f$ s @f$.
   *
   * @see dirichlet_beta for details.
   */
  inline float
  dirichlet_betaf(float __s)
  { return std::__detail::__dirichlet_beta<float>(__s); }

  /**
   * Return the Dirichlet beta function of real argument @f$ s @f$.
   *
   * @see dirichlet_beta for details.
   */
  inline long double
  dirichlet_betal(long double __s)
  { return std::__detail::__dirichlet_beta<long double>(__s); }

  /**
   * Return the Dirichlet beta function of real argument @f$ s @f$.
   *
   * The Dirichlet beta function is defined by:
   * @f[
   *    \beta(s) = \sum_{k=0}^\infty \frac{(-1)^k}{(2k+1)^s}
   * @f]
   * An important reflection formula is:
   * @f[
   *    \beta(1-s) = \left( \frac{2}{\pi}\right)^s \sin(\frac{\pi s}{2})
   *               \Gamma(s) \beta(s)
   * @f]
   *
   * @param __s 
   */
  template<typename _Tp>
    inline _Tp
    dirichlet_beta(_Tp __s)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__dirichlet_beta<__type>(__s);
    }

  // Dirichlet lambda function

  /**
   * Return the Dirichlet lambda function of real argument @f$ s @f$.
   *
   * @see dirichlet_lambda for details.
   */
  inline float
  dirichlet_lambdaf(float __s)
  { return std::__detail::__dirichlet_lambda<float>(__s); }

  /**
   * Return the Dirichlet lambda function of real argument @f$ s @f$.
   *
   * @see dirichlet_lambda for details.
   */
  inline long double
  dirichlet_lambdal(long double __s)
  { return std::__detail::__dirichlet_lambda<long double>(__s); }

  /**
   * Return the Dirichlet lambda function of real argument @f$ s @f$.
   *
   * The Dirichlet lambda function is defined by
   * @f[
   *    \lambda(s) = \sum_{k=0}^\infty \frac{1}{(2k+1)^s}
   *    = \left( 1 - 2^{-s} \right) \zeta(s)
   * @f]
   *
   * @param __s 
   */
  template<typename _Tp>
    inline _Tp
    dirichlet_lambda(_Tp __s)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__dirichlet_lambda<__type>(__s);
    }

  // Clausen S functions

  /**
   * Return the Clausen sine function @f$ S_n(w) @f$ of order @f$ m @f$
   * and @c float argument @f$ w @f$.
   *
   * @see clausen_s for details.
   */
  inline float
  clausen_sf(unsigned int __m, float __w)
  { return std::__detail::__clausen_s<float>(__m, __w); }

  /**
   * Return the Clausen sine function @f$ S_n(w) @f$ of order @f$ m @f$
   * and <tt>long double</tt> argument @f$ w @f$.
   *
   * @see clausen_s for details.
   */
  inline long double
  clausen_sl(unsigned int __m, long double __w)
  { return std::__detail::__clausen_s<long double>(__m, __w); }

  /**
   * Return the Clausen sine function @f$ S_n(w) @f$ of order @f$ m @f$
   * and real argument @f$ w @f$.
   *
   * The Clausen sine function is defined by
   * @f[
   *    S_n(w) = \sum_{k=1}^\infty\frac{\sin(kx)}{k^n}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __m The unsigned integer order
   * @param __w The real argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    clausen_s(unsigned int __m, _Tp __w)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__clausen_s<__type>(__m, __w);
    }

  // Clausen C functions

  /**
   * Return the Clausen cosine function @f$ C_n(w) @f$ of order @f$ m @f$
   * and @c float argument @f$ w @f$.
   *
   * @see clausen_c for details.
   */
  inline float
  clausen_cf(unsigned int __m, float __w)
  { return std::__detail::__clausen_c<float>(__m, __w); }

  /**
   * Return the Clausen cosine function @f$ C_n(w) @f$ of order @f$ m @f$
   * and <tt>long double</tt> argument @f$ w @f$.
   *
   * @see clausen_c for details.
   */
  inline long double
  clausen_cl(unsigned int __m, long double __w)
  { return std::__detail::__clausen_c<long double>(__m, __w); }

  /**
   * Return the Clausen cosine function @f$ C_n(w) @f$ of order @f$ m @f$
   * and real argument @f$ w @f$.
   *
   * The Clausen cosine function is defined by
   * @f[
   *    C_n(w) = \sum_{k=1}^\infty\frac{\cos(kx)}{k^n}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __m The unsigned integer order
   * @param __w The real argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    clausen_c(unsigned int __m, _Tp __w)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__clausen_c<__type>(__m, __w);
    }

  // Clausen functions - real argument

  /**
   * Return the Clausen function @f$ Cl_n(w) @f$ of integer order @f$ m @f$
   * and @c float argument @f$ w @f$.
   *
   * @see clausen for details.
   */
  inline float
  clausenf(unsigned int __m, float __w)
  { return std::__detail::__clausen<float>(__m, __w); }

  /**
   * Return the Clausen function @f$ Cl_n(w) @f$ of integer order @f$ m @f$
   * and <tt>long double</tt> argument @f$ w @f$.
   *
   * @see clausen for details.
   */
  inline long double
  clausenl(unsigned int __m, long double __w)
  { return std::__detail::__clausen<long double>(__m, __w); }

  /**
   * Return the Clausen function @f$ Cl_n(w) @f$ of integer order @f$ m @f$
   * and real argument @f$ w @f$.
   *
   * The Clausen function is defined by
   * @f[
   *    Cl_n(w) = S_n(w) = \sum_{k=1}^\infty\frac{\sin(kx)}{k^n} \mbox{ for even } m
   * = C_n(w) = \sum_{k=1}^\infty\frac{\cos(kx)}{k^n} \mbox{ for odd } m
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __m The integral order
   * @param __w The complex argument
   */
  template<typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tp>
    clausen(unsigned int __m, _Tp __w)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__clausen<__type>(__m, __w);
    }

  // Clausen functions - complex argument

  /**
   * Return the Clausen function @f$ Cl_n(w) @f$ of integer order @f$ m @f$
   * and <tt>std::complex<float></tt> argument @f$ w @f$.
   *
   * @see clausen for details.
   */
  inline std::complex<float>
  clausenf(unsigned int __m, std::complex<float> __w)
  { return std::__detail::__clausen<float>(__m, __w); }

  /**
   * Return the Clausen function @f$ Cl_n(w) @f$ of integer order @f$ m @f$
   * and <tt>std::complex<long double></tt> argument @f$ w @f$.
   *
   * @see clausen for details.
   */
  inline std::complex<long double>
  clausenl(unsigned int __m, std::complex<long double> __w)
  { return std::__detail::__clausen<long double>(__m, __w); }

  /**
   * Return the Clausen function @f$ Cl_n(w) @f$ of integer order @f$ m @f$
   * and complex argument @f$ w @f$.
   *
   * The Clausen function is defined by
   * @f[
   *    Cl_n(w) = S_n(w) = \sum_{k=1}^\infty\frac{\sin(kx)}{k^n} \mbox{ for even } m
   * = C_n(w) = \sum_{k=1}^\infty\frac{\cos(kx)}{k^n} \mbox{ for odd } m
   * @f]
   *
   * @tparam _Tp The real type of the complex components
   * @param __m The integral order
   * @param __w The complex argument
   */
  template<typename _Tp>
    inline std::complex<__gnu_cxx::__promote_fp_t<_Tp>>
    clausen(unsigned int __m, std::complex<_Tp> __w)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__clausen<__type>(__m, __w);
    }

  // Exponential theta_1 functions.

  /**
   * Return the exponential theta-1 function @f$ \theta_1(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_1 for details.
   */
  inline float
  theta_1f(float __nu, float __x)
  { return std::__detail::__theta_1<float>(__nu, __x); }

  /**
   * Return the exponential theta-1 function @f$ \theta_1(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_1 for details.
   */
  inline long double
  theta_1l(long double __nu, long double __x)
  { return std::__detail::__theta_1<long double>(__nu, __x); }

  /**
   * Return the exponential theta-1 function @f$ \theta_1(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>
    theta_1(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__theta_1<__type>(__nu, __x);
    }

  // Exponential theta_2 functions.

  /**
   * Return the exponential theta-2 function @f$ \theta_2(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_2 for details.
   */
  inline float
  theta_2f(float __nu, float __x)
  { return std::__detail::__theta_2<float>(__nu, __x); }

  /**
   * Return the exponential theta-2 function @f$ \theta_2(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_2 for details.
   */
  inline long double
  theta_2l(long double __nu, long double __x)
  { return std::__detail::__theta_2<long double>(__nu, __x); }

  /**
   * Return the exponential theta-2 function @f$ \theta_2(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>
    theta_2(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__theta_2<__type>(__nu, __x);
    }

  // Exponential theta_3 functions.

  /**
   * Return the exponential theta-3 function @f$ \theta_3(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_3 for details.
   */
  inline float
  theta_3f(float __nu, float __x)
  { return std::__detail::__theta_3<float>(__nu, __x); }

  /**
   * Return the exponential theta-3 function @f$ \theta_3(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_3 for details.
   */
  inline long double
  theta_3l(long double __nu, long double __x)
  { return std::__detail::__theta_3<long double>(__nu, __x); }

  /**
   * Return the exponential theta-3 function @f$ \theta_3(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>
    theta_3(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__theta_3<__type>(__nu, __x);
    }

  // Exponential theta_4 functions.

  /**
   * Return the exponential theta-4 function @f$ \theta_4(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_4 for details.
   */
  inline float
  theta_4f(float __nu, float __x)
  { return std::__detail::__theta_4<float>(__nu, __x); }

  /**
   * Return the exponential theta-4 function @f$ \theta_4(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
   *
   * @see theta_4 for details.
   */
  inline long double
  theta_4l(long double __nu, long double __x)
  { return std::__detail::__theta_4<long double>(__nu, __x); }

  /**
   * Return the exponential theta-4 function @f$ \theta_4(\nu,x) @f$
   * of period @f$ nu @f$ and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>
    theta_4(_Tpnu __nu, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpnu, _Tp>;
      return std::__detail::__theta_4<__type>(__nu, __x);
    }

  // Elliptic nome function

  /**
   * Return the elliptic nome function @f$ q(k) @f$
   * of modulus @f$ k @f$.
   *
   * @see ellnome for details.
   */
  inline float
  ellnomef(float __k)
  { return std::__detail::__ellnome<float>(__k); }

  /**
   * Return the elliptic nome function @f$ q(k) @f$
   * of <tt>long double</tt> modulus @f$ k @f$.
   *
   * @see ellnome for details.
   */
  inline long double
  ellnomel(long double __k)
  { return std::__detail::__ellnome<long double>(__k); }

  /**
   * Return the elliptic nome function @f$ q(k) @f$ of modulus @f$ k @f$.
   *
   * The elliptic nome function is defined by
   * @f[
   *    q(k) = \exp \left(-\pi\frac{K(k)}{K(\sqrt{1-k^2})} \right)
   * @f]
   * where @f$ K(k) @f$ is the complete elliptic function of the first kind.
   *
   * @tparam _Tp The real type of the modulus
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   */
  template<typename _Tp>
    inline _Tp
    ellnome(_Tp __k)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tp>;
      return std::__detail::__ellnome<__type>(__k);
    }

  // Neville theta_s functions.

  /**
   * Return the Neville theta-s function @f$ \theta_s(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_s for details.
   */
  inline float
  theta_sf(float __k, float __x)
  { return std::__detail::__theta_s<float>(__k, __x); }

  /**
   * Return the Neville theta-s function @f$ \theta_s(k,x) @f$
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_s for details.
   */
  inline long double
  theta_sl(long double __k, long double __x)
  { return std::__detail::__theta_s<long double>(__k, __x); }

  /**
   * Return the Neville theta-s function @f$ \theta_s(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tpk, _Tp>
    theta_s(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpk, _Tp>;
      return std::__detail::__theta_s<__type>(__k, __x);
    }

  // Neville theta_c functions.

  /**
   * Return the Neville theta-c function @f$ \theta_c(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_c for details.
   */
  inline float
  theta_cf(float __k, float __x)
  { return std::__detail::__theta_c<float>(__k, __x); }

  /**
   * Return the Neville theta-c function @f$ \theta_c(k,x) @f$
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_c for details.
   */
  inline long double
  theta_cl(long double __k, long double __x)
  { return std::__detail::__theta_c<long double>(__k, __x); }

  /**
   * Return the Neville theta-c function @f$ \theta_c(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
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
    inline __gnu_cxx::__promote_fp_t<_Tpk, _Tp>
    theta_c(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpk, _Tp>;
      return std::__detail::__theta_c<__type>(__k, __x);
    }

  // Neville theta_d functions.

  /**
   * Return the Neville theta-d function @f$ \theta_d(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_d for details.
   */
  inline float
  theta_df(float __k, float __x)
  { return std::__detail::__theta_d<float>(__k, __x); }

  /**
   * Return the Neville theta-d function @f$ \theta_d(k,x) @f$
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_d for details.
   */
  inline long double
  theta_dl(long double __k, long double __x)
  { return std::__detail::__theta_d<long double>(__k, __x); }

  /**
   * Return the Neville theta-d function @f$ \theta_d(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * The Neville theta-d function is defined by
   * @f[
   *    \theta_d(k,x) = 
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __x The argument
   */
  template<typename _Tpk, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tpk, _Tp>
    theta_d(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpk, _Tp>;
      return std::__detail::__theta_d<__type>(__k, __x);
    }

  // Neville theta_n functions.

  /**
   * Return the Neville theta-n function @f$ \theta_n(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_n for details.
   */
  inline float
  theta_nf(float __k, float __x)
  { return std::__detail::__theta_n<float>(__k, __x); }

  /**
   * Return the Neville theta-n function @f$ \theta_n(k,x) @f$
   * of <tt>long double</tt> modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * @see theta_n for details.
   */
  inline long double
  theta_nl(long double __k, long double __x)
  { return std::__detail::__theta_n<long double>(__k, __x); }

  /**
   * Return the Neville theta-n function @f$ \theta_n(k,x) @f$
   * of modulus @f$ k @f$ and argument @f$ x @f$.
   *
   * The Neville theta-n function is defined by
   * @f[
   *    \theta_n(k,x) = 
   * @f]
   *
   * @param __k The modulus @f$ -1 <= k <= +1 @f$
   * @param __x The argument
   */
  template<typename _Tpk, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tpk, _Tp>
    theta_n(_Tpk __k, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tpk, _Tp>;
      return std::__detail::__theta_n<__type>(__k, __x);
    }

  // Owens T functions.

  /**
   * Return the Owens T function @f$ T(h,a) @f$
   * of shape factor @f$ h @f$ and integration limit @f$ a @f$.
   *
   * @see owens_t for details.
   */
  inline float
  owens_tf(float __h, float __a)
  { return std::__detail::__owens_t<float>(__h, __a); }

  /**
   * Return the Owens T function @f$ T(h,a) @f$ of <tt>long double</tt>
   * shape factor @f$ h @f$ and integration limit @f$ a @f$.
   *
   * @see owens_t for details.
   */
  inline long double
  owens_tl(long double __h, long double __a)
  { return std::__detail::__owens_t<long double>(__h, __a); }

  /**
   * Return the Owens T function @f$ T(h,a) @f$ of shape factor @f$ h @f$
   * and integration limit @f$ a @f$.
   *
   * The Owens T function is defined by
   * @f[
   *    T(h,a) = \frac{1}{2\pi}\int_0^a
   *           \frac{\exp\left[-\frac{1}{2}h^2(1+x^2)\right]}{1+x^2} dx
   * @f]
   *
   * @param __h The shape factor
   * @param __a The integration limit
   */
  template<typename _Tph, typename _Tpa>
    inline __gnu_cxx::__promote_fp_t<_Tph, _Tpa>
    owens_t(_Tph __h, _Tpa __a)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tph, _Tpa>;
      return std::__detail::__owens_t<__type>(__h, __a);
    }

  // Fermi-Dirac integrals.

  inline float
  fermi_diracf(float __s, float __x)
  { return std::__detail::__fermi_dirac<float>(__s, __x); }

  inline long double
  fermi_diracl(long double __s, long double __x)
  { return std::__detail::__fermi_dirac<long double>(__s, __x); }

  template<typename _Tps, typename _Tp>
    inline __gnu_cxx::__promote_fp_t<_Tps, _Tp>
    fermi_dirac(_Tps __s, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tps, _Tp>;
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
    inline __gnu_cxx::__promote_fp_t<_Tps, _Tp>
    bose_einstein(_Tps __s, _Tp __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tps, _Tp>;
      return std::__detail::__bose_einstein<__type>(__s, __x);
    }

  // Reperiodized sine function.

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see sin_pi for more details.
   */
  inline float
  sin_pif(float __x)
  { return std::__detail::__sin_pi<float>(__x); }

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see sin_pi for more details.
   */
  inline long double
  sin_pil(long double __x)
  { return std::__detail::__sin_pi<long double>(__x); }

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized sine function is defined by:
   * @f[
   * 	\sin_\pi(x) = \sin(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sin_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__sin_pi<__type>(__x);
    }

  // Reperiodized hyperbolic sine function.

  /**
   * Return the reperiodized hyperbolic sine function @f$ \sinh_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see sinh_pi for more details.
   */
  inline float
  sinh_pif(float __x)
  { return std::__detail::__sinh_pi<float>(__x); }

  /**
   * Return the reperiodized hyperbolic sine function @f$ \sinh_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see sinh_pi for more details.
   */
  inline long double
  sinh_pil(long double __x)
  { return std::__detail::__sinh_pi<long double>(__x); }

  /**
   * Return the reperiodized hyperbolic sine function @f$ \sinh_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized hyperbolic sine function is defined by:
   * @f[
   * 	\sinh_\pi(x) = \sinh(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sinh_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__sinh_pi<__type>(__x);
    }

  // Reperiodized cosine function.

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see cos_pi for more details.
   */
  inline float
  cos_pif(float __x)
  { return std::__detail::__cos_pi<float>(__x); }

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see cos_pi for more details.
   */
  inline long double
  cos_pil(long double __x)
  { return std::__detail::__cos_pi<long double>(__x); }

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized cosine function is defined by:
   * @f[
   * 	\cos_\pi(x) = \cos(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    cos_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__cos_pi<__type>(__x);
    }

  // Reperiodized hyperbolic cosine function.

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see cosh_pi for more details.
   */
  inline float
  cosh_pif(float __x)
  { return std::__detail::__cosh_pi<float>(__x); }

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see cosh_pi for more details.
   */
  inline long double
  cosh_pil(long double __x)
  { return std::__detail::__cosh_pi<long double>(__x); }

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized hyperbolic cosine function is defined by:
   * @f[
   * 	\cosh_\pi(x) = \cosh(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    cosh_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__cosh_pi<__type>(__x);
    }

  // Reperiodized tangent function.

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see tan_pi for more details.
   */
  inline float
  tan_pif(float __x)
  { return std::__detail::__tan_pi<float>(__x); }

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see tan_pi for more details.
   */
  inline long double
  tan_pil(long double __x)
  { return std::__detail::__tan_pi<long double>(__x); }

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized tangent function is defined by:
   * @f[
   * 	\tan_\pi(x) = \tan(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    tan_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__tan_pi<__type>(__x);
    }

  // Reperiodized hyperbolic tangent function.

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see tanh_pi for more details.
   */
  inline float
  tanh_pif(float __x)
  { return std::__detail::__tanh_pi<float>(__x); }

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see tanh_pi for more details.
   */
  inline long double
  tanh_pil(long double __x)
  { return std::__detail::__tanh_pi<long double>(__x); }

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized hyperbolic tangent function is defined by:
   * @f[
   * 	\tanh_\pi(x) = \tanh(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    tanh_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__tanh_pi<__type>(__x);
    }

  /**
   * Return both the sine and the cosine of a @c float argument.
   */
  inline __gnu_cxx::__sincos_t<float>
  sincosf(float __x)
  { return std::__detail::__sincos<float>(__x); }

  /**
   * Return both the sine and the cosine of a <tt> long double </tt> argument.
   *
   * @see sincos for details.
   */
  inline __gnu_cxx::__sincos_t<long double>
  sincosl(long double __x)
  { return std::__detail::__sincos<long double>(__x); }

  /**
   * Return both the sine and the cosine of a @c double argument.
   *
   * @see sincos for details.
   */
  inline __gnu_cxx::__sincos_t<double>
  sincos(double __x)
  { return std::__detail::__sincos<double>(__x); }

  /**
   * Return both the sine and the cosine of a reperiodized argument.
   * @f[
   *   sincos(x) = {\sin(x), \cos(x)}
   * @f]
   */
  template<typename _Tp>
    inline __gnu_cxx::__sincos_t<_Tp>
    sincos(_Tp __x)
    { return std::__detail::__sincos<_Tp>(__x); }

  /**
   * Return both the sine and the cosine of a reperiodized @c float argument.
   *
   * @see sincos_pi for details.
   */
  inline __gnu_cxx::__sincos_t<float>
  sincos_pif(float __x)
  { return std::__detail::__sincos_pi<float>(__x); }

  /**
   * Return both the sine and the cosine of a reperiodized
   * <tt> long double </tt> argument.
   *
   * @see sincos_pi for details.
   */
  inline __gnu_cxx::__sincos_t<long double>
  sincos_pil(long double __x)
  { return std::__detail::__sincos_pi<long double>(__x); }

  /**
   * Return both the sine and the cosine of a reperiodized real argument.
   * 
   * @f[
   *   sincos_\pi(x) = {\sin(\pi x), \cos(\pi x)}
   * @f]
   */
  template<typename _Tp>
    inline __gnu_cxx::__sincos_t<_Tp>
    sincos_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__sincos_pi<__type>(__x);
    }

#endif // __cplusplus >= 201103L

  /** @} */ // gnu_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#pragma GCC visibility pop

#endif // _GLIBCXX_BITS_SPECFUN_H
