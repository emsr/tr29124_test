
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

/** @file bits/specfun.h
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SPECFUN_H
#define SPECFUN_H 1

#include <emsr/sf_trig.h>
#include <emsr/sf_cardinal.h>
#include <emsr/sf_theta.h>
#include <emsr/sf_zeta.h>
#include <emsr/sf_prime.h>
#include <emsr/sf_ellint.h>
#include <emsr/sf_gamma.h>
#include <emsr/sf_beta.h>
#include <emsr/sf_bernoulli.h>
#include <emsr/sf_euler.h>
#include <emsr/sf_stirling.h>
#include <emsr/sf_owens_t.h>
#include <emsr/sf_distributions.h>
#include <emsr/sf_laguerre.h>
#include <emsr/sf_legendre.h>
#include <emsr/sf_hermite.h>
#include <emsr/sf_gegenbauer.h>
#include <emsr/sf_jacobi.h>
#include <emsr/sf_chebyshev.h>
#include <emsr/sf_expint.h>
#include <emsr/sf_hypint.h>
#include <emsr/sf_trigint.h>
#include <emsr/sf_coulomb.h>
#include <emsr/sf_dawson.h>
#include <emsr/sf_fresnel.h>
#include <emsr/sf_lerch.h>
#include <emsr/sf_mittag_leffler.h>
#include <emsr/sf_bessel.h>
#include <emsr/sf_mod_bessel.h>
#include <emsr/sf_hyperg.h>
#include <emsr/sf_polylog.h>
#include <emsr/sf_airy.h>
#include <emsr/sf_hankel.h>

  /**
   * @defgroup emsr_specfun Mathematical Special Functions
   * @ingroup numerics
   *
   * @section emsr_specfun_desc_desc Mathematical Special Functions
   *
   * A collection of advanced mathematical special functions,
   * defined by ISO/IEC IS 29124 and then added to ISO C++ 2017.
   *
   *
   * @subsection emsr_specfun_intro Introduction and History
   * The first significant library upgrade on the road to C++2011,
   * <a href="http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1836.pdf">
   * TR1</a>, included a set of 23 mathematical functions that significantly
   * extended the standard transcendental functions inherited from C and
   * declared in @<cmath@>.
   *
   * Although most components from TR1 were eventually adopted for C++11 these
   * math functions were left behind out of concern for implementability.
   * The math functions were published as a separate international standard
   * <a href="http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2010/n3060.pdf">
   * IS 29124 - Extensions to the C++ Library to Support Mathematical Special
   * Functions</a>.
   *
   * Follow-up proosals for new special functions have also been published:
   * <a href="http://open-std.org/JTC1/SC22/WG21/docs/papers/2013/n3494.pdf">
   * A proposal to add special mathematical functions according to
   * the ISO/IEC 80000-2:2009 standard, Vincent Reverdy</a>.
   *
   * <a href="http://open-std.org/JTC1/SC22/WG21/docs/papers/2004/n1668.pdf">
   * A Proposal to add Mathematical Functions for Statistics
   * to the C++ Standard Library, Paul A Bristow</a>.
   *
   * <a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/p0081r0.pdf">
   * A proposal to add sincos to the standard library, Paul Dreik</a>.

   * For C++17 these functions were incorporated into the main standard.
   *
   * @subsection emsr_specfun_contents Contents
   * The following functions (standardized in C++17) are implemented in namespace @c emsr:
   * - @ref emsr::assoc_laguerre    "assoc_laguerre - Associated Laguerre functions"
   * - @ref emsr::assoc_legendre    "assoc_legendre - Associated Legendre functions"
   * - @ref emsr::assoc_legendre_q  "assoc_legendre_q - Associated Legendre functions of the second kind"
   * - @ref emsr::beta              "beta - Beta functions"
   * - @ref emsr::comp_ellint_1     "comp_ellint_1 - Complete elliptic functions of the first kind"
   * - @ref emsr::comp_ellint_2     "comp_ellint_2 - Complete elliptic functions of the second kind"
   * - @ref emsr::comp_ellint_3     "comp_ellint_3 - Complete elliptic functions of the third kind"
   * - @ref emsr::cyl_bessel_i      "cyl_bessel_i - Regular modified cylindrical Bessel functions"
   * - @ref emsr::cyl_bessel_j      "cyl_bessel_j - Cylindrical Bessel functions of the first kind"
   * - @ref emsr::cyl_bessel_k      "cyl_bessel_k - Irregular modified cylindrical Bessel functions"
   * - @ref emsr::cyl_neumann       "cyl_neumann - Cylindrical Neumann functions or Cylindrical Bessel functions of the second kind"
   * - @ref emsr::ellint_1          "ellint_1 - Incomplete elliptic functions of the first kind"
   * - @ref emsr::ellint_2          "ellint_2 - Incomplete elliptic functions of the second kind"
   * - @ref emsr::ellint_3          "ellint_3 - Incomplete elliptic functions of the third kind"
   * - @ref emsr::expint            "expint - The exponential integral"
   * - @ref emsr::hermite           "hermite - Hermite polynomials"
   * - @ref emsr::laguerre          "laguerre - Laguerre functions"
   * - @ref emsr::legendre          "legendre - Legendre polynomials"
   * - @ref emsr::riemann_zeta      "riemann_zeta - The Riemann zeta function"
   * - @ref emsr::sph_bessel        "sph_bessel - Spherical Bessel functions"
   * - @ref emsr::sph_legendre      "sph_legendre - Spherical Legendre functions"
   * - @ref emsr::sph_neumann       "sph_neumann - Spherical Neumann functions"
   *
   * The hypergeometric functions were stricken from the TR29124 and C++17
   * versions of this math library because of implementation concerns.
   * However, since they were in the TR1 version and since they are popular
   * we kept them as an extension in namespace @c __gnu_cxx:
   * - @ref emsr::conf_hyperg       "conf_hyperg - - Confluent hypergeometric functions"
   * - @ref emsr::hyperg            "hyperg - Hypergeometric functions"
   *
   * - @ref emsr::airy_ai           "airy_ai - Airy functions of the first kind"
   * - @ref emsr::airy_bi           "airy_bi - Airy functions of the second kind"
   * - @ref emsr::bell              "bell - Bell numbers and polynomials"
   * - @ref emsr::bernoulli         "bernoulli - Bernoulli polynomials"
   * - @ref emsr::binomial          "binomial - Binomial coefficients"
   * - @ref emsr::bose_einstein     "bose_einstein - Bose-Einstein integrals"
   * - @ref emsr::chebyshev_t       "chebyshev_t - Chebyshev polynomials of the first kind"
   * - @ref emsr::chebyshev_u       "chebyshev_u - Chebyshev polynomials of the second kind"
   * - @ref emsr::chebyshev_v       "chebyshev_v - Chebyshev polynomials of the third kind"
   * - @ref emsr::chebyshev_w       "chebyshev_w - Chebyshev polynomials of the fourth kind"
   * - @ref emsr::clausen           "clausen - Clausen integrals"
   * - @ref emsr::clausen_cl        "clausen_cl - Clausen cosine integrals"
   * - @ref emsr::clausen_sl        "clausen_sl - Clausen sine integrals"
   * - @ref emsr::comp_ellint_d     "comp_ellint_d - Incomplete Legendre D elliptic integral"
   * - @ref emsr::conf_hyperg_lim   "conf_hyperg_lim - Confluent hypergeometric limit functions"
   * - @ref emsr::cos_pi            "cos_pi - Reperiodized cosine function."
   * - @ref emsr::cosh_pi           "cosh_pi - Reperiodized hyperbolic cosine function."
   * - @ref emsr::coshint           "coshint - Hyperbolic cosine integral"
   * - @ref emsr::cosint            "cosint - Cosine integral"
   * - @ref emsr::cyl_hankel_1      "cyl_hankel_1 - Cylindrical Hankel functions of the first kind"
   * - @ref emsr::cyl_hankel_2      "cyl_hankel_2 - Cylindrical Hankel functions of the second kind"
   * - @ref emsr::dawson            "dawson - Dawson integrals"
   * - @ref emsr::debye             "debye - Debye functions"
   * - @ref emsr::digamma           "digamma - Digamma or psi function"
   * - @ref emsr::dilog             "dilog - Dilogarithm functions"
   * - @ref emsr::dirichlet_beta    "dirichlet_beta - Dirichlet beta function"
   * - @ref emsr::dirichlet_eta     "dirichlet_eta - Dirichlet beta function"
   * - @ref emsr::dirichlet_lambda  "dirichlet_lambda - Dirichlet lambda function"
   * - @ref emsr::double_factorial  "double_factorial - Double factorials"
   * - @ref emsr::ellint_d          "ellint_d - Legendre D elliptic integrals"
   * - @ref emsr::ellint_rc 	    "ellint_rc - Carlson elliptic functions R_C"
   * - @ref emsr::ellint_rd 	    "ellint_rd - Carlson elliptic functions R_D"
   * - @ref emsr::ellint_rf 	    "ellint_rf - Carlson elliptic functions R_F"
   * - @ref emsr::ellint_rg 	    "ellint_rg - Carlson elliptic functions R_G"
   * - @ref emsr::ellint_rj 	    "ellint_rj - Carlson elliptic functions R_J"
   * - @ref emsr::ellnome  	    "ellnome - Elliptic nome"
   * - @ref emsr::euler  	    "euler - Euler numbers"
   * - @ref emsr::euler  	    "euler - Euler polynomials"
   * - @ref emsr::eulerian_1        "eulerian_1 - Eulerian numbers of the first kind"
   * - @ref emsr::eulerian_2        "eulerian_2 - Eulerian numbers of the second kind"
   * - @ref emsr::expint            "expint - Exponential integrals"
   * - @ref emsr::factorial         "factorial - Factorials"
   * - @ref emsr::falling_factorial "falling_factorial - Falling factorials"
   * - @ref emsr::fermi_dirac       "fermi_dirac - Fermi-Dirac integrals"
   * - @ref emsr::fresnel_c         "fresnel_c - Fresnel cosine integrals"
   * - @ref emsr::fresnel_s         "fresnel_s - Fresnel sine integrals"
   * - @ref emsr::gamma_p           "gamma_p - Regularized lower incomplete gamma functions"
   * - @ref emsr::gamma_q           "gamma_q - Regularized upper incomplete gamma functions"
   * - @ref emsr::gamma_reciprocal  "gamma_reciprocal - Reciprocal gamma function"
   * - @ref emsr::gegenbauer  	    "gegenbauer - Gegenbauer polynomials"
   * - @ref emsr::heuman_lambda     "heuman_lambda - Heuman lambda functions"
   * - @ref emsr::hurwitz_zeta      "hurwitz_zeta - Hurwitz zeta functions"
   * - @ref emsr::ibeta  	    "ibeta - Regularized incomplete beta functions"
   * - @ref emsr::jacobi  	    "jacobi - Jacobi polynomials"
   * - @ref emsr::jacobi_sn 	    "jacobi_sn - Jacobi sine amplitude functions"
   * - @ref emsr::jacobi_cn 	    "jacobi_cn - Jacobi cosine amplitude functions"
   * - @ref emsr::jacobi_dn 	    "jacobi_dn - Jacobi delta amplitude functions"
   * - @ref emsr::jacobi_theta_1    "theta_1 - Jacobi theta function 1"
   * - @ref emsr::jacobi_theta_2    "theta_2 - Jacobi theta function 2"
   * - @ref emsr::jacobi_theta_3    "theta_3 - Jacobi theta function 3"
   * - @ref emsr::jacobi_theta_4    "theta_4 - Jacobi theta function 4"
   * - @ref emsr::jacobi_zeta       "jacobi_zeta - Jacobi zeta functions"
   * - @ref emsr::lah               "lah - Lah numbers"
   * - @ref emsr::lbinomial         "lbinomial - Log binomial coefficients"
   * - @ref emsr::ldouble_factorial "ldouble_factorial - Log double factorials"
   * - @ref emsr::legendre_q        "legendre_q - Legendre functions of the second kind"
   * - @ref emsr::lerch_phi         "lerch_phi - The Lerch transcendent"
   * - @ref emsr::lfactorial        "lfactorial - Log factorials"
   * - @ref emsr::lfalling_factorial "lfalling_factorial - Log falling factorials"
   * - @ref emsr::lgamma "lgamma - Log gamma for complex arguments"
   * - @ref emsr::lrising_factorial "lrising_factorial - Log rising factorials"
   * - @ref emsr::mittag_leffler    "mittag_leffler - Mittag-Leffler functions"
   * - @ref emsr::owens_t           "owens_t - Owens T functions"
   * - @ref emsr::periodic_zeta     "periodic_zeta - Periodic zeta functions"
   * - @ref emsr::radpoly           "radpoly - Radial polynomials"
   * - @ref emsr::rising_factorial  "rising_factorial - Rising factorials"
   * - @ref emsr::sinhc             "sinhc - Hyperbolic sinus cardinal function"
   * - @ref emsr::sinhc_pi          "sinhc_pi - Reperiodized hyperbolic sinus cardinal function"
   * - @ref emsr::sinc              "sinc - Normalized sinus cardinal function"
   * - @ref emsr::sincos            "sincos - Sine + cosine function"
   * - @ref emsr::sincos_pi         "sincos_pi - Reperiodized sine + cosine function"
   * - @ref emsr::sin_pi            "sin_pi - Reperiodized sine function."
   * - @ref emsr::sinh_pi           "sinh_pi - Reperiodized hyperbolic sine function."
   * - @ref emsr::sinc_pi           "sinc_pi - Sinus cardinal function"
   * - @ref emsr::sinhint           "sinhint - Hyperbolic sine integral"
   * - @ref emsr::sinint            "sinint - Sine integral"
   * - @ref emsr::sph_bessel_i      "sph_bessel_i - Spherical regular modified Bessel functions"
   * - @ref emsr::sph_bessel_k      "sph_bessel_k - Spherical iregular modified Bessel functions"
   * - @ref emsr::sph_hankel_1      "sph_hankel_1 - Spherical Hankel functions of the first kind"
   * - @ref emsr::sph_hankel_2      "sph_hankel_2 - Spherical Hankel functions of the first kind"
   * - @ref emsr::sph_harmonic      "sph_harmonic - Spherical"
   * - @ref emsr::stirling_1        "stirling_1 - Stirling numbers of the first kind"
   * - @ref emsr::stirling_2        "stirling_2 - Stirling numbers of the second kind"
   * - @ref emsr::tan_pi            "tan_pi - Reperiodized tangent function."
   * - @ref emsr::tanh_pi           "tanh_pi - Reperiodized hyperbolic tangent function."
   * - @ref emsr::tgamma            "tgamma - Gamma for complex arguments"
   * - @ref emsr::tgamma            "tgamma - Upper incomplete gamma functions"
   * - @ref emsr::tgamma_lower      "tgamma_lower - Lower incomplete gamma functions"
   * - @ref emsr::theta_1 	    "theta_1 - Exponential theta function 1"
   * - @ref emsr::theta_2 	    "theta_2 - Exponential theta function 2"
   * - @ref emsr::theta_3 	    "theta_3 - Exponential theta function 3"
   * - @ref emsr::theta_4 	    "theta_4 - Exponential theta function 4"
   * - @ref emsr::tricomi_u         "tricomi_u - Tricomi confluent hypergeometric function"
   * - @ref emsr::zernike           "zernike - Zernike polynomials"
   *
   * <!-- @subsection emsr_specfun_general General Features -->
   *
   * @subsection emsr_specfun_promotion Argument Promotion
   * The arguments suppled to the non-suffixed functions will be promoted
   * according to the following rules:
   * 1. If any argument intended to be floating point is given an integral value
   * That integral value is promoted to double.
   * 2. All floating point arguments are promoted up to the largest floating
   *    point precision among them.
   *
   * @subsection emsr_specfun_NaN NaN Arguments
   * If any of the floating point arguments supplied to these functions is
   * invalid or NaN (std::numeric_limits<Tp>::quiet_NaN),
   * the value NaN is returned.
   *
   * @subsection emsr_specfun_impl Implementation
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
   * @subsection emsr_specfun_testing Testing
   *
   * These functions have been tested against equivalent implementations
   * from the <a href="http://www.gnu.org/software/gsl">
   * Gnu Scientific Library, GSL</a> and
   * <a href="http://www.boost.org/doc/libs/1_60_0/libs/math/doc/html/index.html>Boost</a>
   * and the ratio
   * @f[
   *   \frac{|f - f_{test}|}{|f_{test}|}
   * @f]
   * is generally found to be within 10<sup>-15</sup> for 64-bit double
   * on linux-x86_64 systems over most of the ranges of validity.
   * 
   * @todo Provide accuracy comparisons on a per-function basis for a small
   *       number of targets.
   *
   * @subsection emsr_specfun_bibliography General Bibliography
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
   * @see The PlanetMath website has a wealth of material.  In particular
   * https://planetmath.org/msc.html#33_Special_functions
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

  /** @} */ // emsr_specfun

#endif // SPECFUN_H
