// Copyright (C) 2015-2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
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

// specfun_testcase.h

//
// These are little PODs for special function inputs and
// expexted results for the testsuite.
//

#ifndef _GLIBCXX_SPECFUN_TESTCASE_H
#define _GLIBCXX_SPECFUN_TESTCASE_H

#include <complex>

// Associated Laguerre polynomials.
template<typename _Tp>
  struct testcase_assoc_laguerre
  {
    _Tp f0;
    unsigned int n;
    unsigned int m;
    _Tp x;
    _Tp f;
  };

// Associated Legendre functions.
template<typename _Tp>
  struct testcase_assoc_legendre
  {
    _Tp f0;
    unsigned int l;
    unsigned int m;
    _Tp x;
    _Tp f;
  };

// Beta function.
template<typename _Tp>
  struct testcase_beta
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp f;
  };

// Complete elliptic integrals of the first kind.
template<typename _Tp>
  struct testcase_comp_ellint_1
  {
    _Tp f0;
    _Tp k;
    _Tp f;
  };

// Complete elliptic integrals of the second kind.
template<typename _Tp>
  struct testcase_comp_ellint_2
  {
    _Tp f0;
    _Tp k;
    _Tp f;
  };

// Complete elliptic integrals of the third kind.
template<typename _Tp>
  struct testcase_comp_ellint_3
  {
    _Tp f0;
    _Tp k;
    _Tp nu;
    _Tp f;
  };

// Complete elliptic D integrals.
template<typename _Tp>
  struct testcase_comp_ellint_d
  {
    _Tp f0;
    _Tp k;
    _Tp f;
  };

// Confluent hypergeometric functions.
template<typename _Tp>
  struct testcase_conf_hyperg
  {
    _Tp f0;
    _Tp a;
    _Tp c;
    _Tp x;
    _Tp f;
  };

// Confluent hypergeometric functions.
template<typename _Tp>
  struct testcase_conf_hyperg_lim
  {
    _Tp f0;
    _Tp c;
    _Tp x;
    _Tp f;
  };

// Generic cylindrical Bessel functions.
template<typename _Tp>
  struct testcase_cyl_bessel
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Regular modified cylindrical Bessel functions.
template<typename _Tp>
  struct testcase_cyl_bessel_i
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Cylindrical Bessel functions (of the first kind).
template<typename _Tp>
  struct testcase_cyl_bessel_j
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Irregular modified cylindrical Bessel functions.
template<typename _Tp>
  struct testcase_cyl_bessel_k
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Cylindrical Neumann functions.
template<typename _Tp>
  struct testcase_cyl_neumann
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Elliptic integrals of the first kind.
template<typename _Tp>
  struct testcase_ellint_1
  {
    _Tp f0;
    _Tp k;
    _Tp phi;
    _Tp f;
  };

// Elliptic integrals of the second kind.
template<typename _Tp>
  struct testcase_ellint_2
  {
    _Tp f0;
    _Tp k;
    _Tp phi;
    _Tp f;
  };

// Elliptic integrals of the third kind.
template<typename _Tp>
  struct testcase_ellint_3
  {
    _Tp f0;
    _Tp k;
    _Tp nu;
    _Tp phi;
    _Tp f;
  };

// Elliptic D integrals.
template<typename _Tp>
  struct testcase_ellint_d
  {
    _Tp f0;
    _Tp k;
    _Tp phi;
    _Tp f;
  };

// Heuman lambda functions.
template<typename _Tp>
  struct testcase_heuman_lambda
  {
    _Tp f0;
    _Tp k;
    _Tp phi;
    _Tp f;
  };

// Jacobi zeta functions.
template<typename _Tp>
  struct testcase_jacobi_zeta
  {
    _Tp f0;
    _Tp k;
    _Tp phi;
    _Tp f;
  };

// Exponential integral.
template<typename _Tp>
  struct testcase_expint
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Hermite polynomials
template<typename _Tp>
  struct testcase_hermite
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Hypergeometric functions.
template<typename _Tp>
  struct testcase_hyperg
  {
    _Tp f0;
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp x;
    _Tp f;
  };

// Laguerre polynomials.
template<typename _Tp>
  struct testcase_laguerre
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Legendre polynomials.
template<typename _Tp>
  struct testcase_legendre
  {
    _Tp f0;
    unsigned int l;
    _Tp x;
    _Tp f;
  };

// Riemann zeta function.
template<typename _Tp>
  struct testcase_riemann_zeta
  {
    _Tp f0;
    _Tp s;
    _Tp f;
  };

// Hurwitz zeta function.
template<typename _Tp>
  struct testcase_hurwitz_zeta
  {
    _Tp f0;
    _Tp s;
    _Tp a;
    _Tp f;
  };

// Spherical Bessel functions.
template<typename _Tp>
  struct testcase_sph_bessel
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Regular modified spherical Bessel functions.
template<typename _Tp>
  struct testcase_sph_bessel_i
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Irregular modified spherical Bessel functions.
template<typename _Tp>
  struct testcase_sph_bessel_k
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Spherical Legendre functions.
template<typename _Tp>
  struct testcase_sph_legendre
  {
    _Tp f0;
    unsigned int l;
    unsigned int m;
    _Tp theta;
    _Tp f;
  };

// Spherical Neumann functions.
template<typename _Tp>
  struct testcase_sph_neumann
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Airy Ai functions.
template<typename _Tp>
  struct testcase_airy_ai
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Airy Bi functions.
template<typename _Tp>
  struct testcase_airy_bi
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Upper incomplete gamma functions.
template<typename _Tp>
  struct testcase_tgamma
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Lower incomplete gamma functions.
template<typename _Tp>
  struct testcase_tgamma_lower
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Regularized upper incomplete gamma functions.
template<typename _Tp>
  struct testcase_qgamma
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Regularized lower incomplete gamma functions.
template<typename _Tp>
  struct testcase_pgamma
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Dilogarithm functions.
template<typename _Tp>
  struct testcase_dilog
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Digamma functions.
template<typename _Tp>
  struct testcase_gamma
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rc
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rd
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_comp_ellint_rf
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rf
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_comp_ellint_rg
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rg
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rj
  {
    _Tp f0;
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp p;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_psi
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinint
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_cosint
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinhint
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_coshint
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_dawson
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_jacobi_sn
  {
    _Tp f0;
    _Tp k;
    _Tp u;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_jacobi_cn
  {
    _Tp f0;
    _Tp k;
    _Tp u;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_jacobi_dn
  {
    _Tp f0;
    _Tp k;
    _Tp u;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_expint_en
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_fresnel_c
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_fresnel_s
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinc
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinc_pi
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Log upper Pochhammer symbol.
template<typename _Tp>
  struct testcase_lpochhammer
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Upper Pochhammer symbols.
template<typename _Tp>
  struct testcase_pochhammer
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Log lower Pochhammer symbol.
template<typename _Tp>
  struct testcase_lpochhammer_lower
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Lower Pochhammer symbol.
template<typename _Tp>
  struct testcase_pochhammer_lower
  {
    _Tp f0;
    _Tp a;
    _Tp x;
    _Tp f;
  };

// Legendre functions of the second kind.
template<typename _Tp>
  struct testcase_legendre_q
  {
    _Tp f0;
    unsigned int l;
    _Tp x;
    _Tp f;
  };

// Factorial.
template<typename _Tp>
  struct testcase_factorial
  {
    _Tp f0;
    unsigned int n;
    _Tp f;
  };

// Log factorial.
template<typename _Tp>
  struct testcase_lfactorial
  {
    _Tp f0;
    unsigned int n;
    _Tp f;
  };

// Double factorial.
template<typename _Tp>
  struct testcase_double_factorial
  {
    _Tp f0;
    unsigned int n;
    _Tp f;
  };

// Log double factorial.
template<typename _Tp>
  struct testcase_ldouble_factorial
  {
    _Tp f0;
    unsigned int n;
    _Tp f;
  };

// Binomial coefficient.
template<typename _Tp>
  struct testcase_bincoef
  {
    _Tp f0;
    unsigned int n;
    unsigned int k;
    _Tp f;
  };

// Log binomial coefficient.
template<typename _Tp>
  struct testcase_lbincoef
  {
    _Tp f0;
    unsigned int n;
    unsigned int k;
    _Tp f;
  };

// Gegenbauer polynomials.
template<typename _Tp>
  struct testcase_gegenbauer
  {
    _Tp f0;
    unsigned int n;
    _Tp alpha;
    _Tp x;
    _Tp f;
  };

// Chebyshev polynomials of the first kind.
template<typename _Tp>
  struct testcase_chebyshev_t
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Chebyshev polynomials of the second kind.
template<typename _Tp>
  struct testcase_chebyshev_u
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Chebyshev polynomials of the third kind.
template<typename _Tp>
  struct testcase_chebyshev_v
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Chebyshev polynomials of the fourth kind.
template<typename _Tp>
  struct testcase_chebyshev_w
  {
    _Tp f0;
    unsigned int n;
    _Tp x;
    _Tp f;
  };

// Zernike polynomials.
template<typename _Tp>
  struct testcase_zernike
  {
    _Tp f0;
    unsigned int n;
    int m;
    _Tp rho;
    _Tp phi;
    _Tp f;
  };

// Radial polynomials.
template<typename _Tp>
  struct testcase_radpoly
  {
    _Tp f0;
    unsigned int n;
    int m;
    _Tp rho;
    _Tp f;
  };

// Incomplete beta function.
template<typename _Tp>
  struct testcase_ibeta
  {
    _Tp f0;
    _Tp a;
    _Tp b;
    _Tp x;
    _Tp f;
  };

// Complementary beta function.
template<typename _Tp>
  struct testcase_ibetac
  {
    _Tp f0;
    _Tp a;
    _Tp b;
    _Tp x;
    _Tp f;
  };

// Cylindrical Hankel functions.
template<typename _Tp>
  struct testcase_cyl_hankel_1
  {
    std::complex<_Tp> f0;
    _Tp nu;
    _Tp x;
    std::complex<_Tp> f;
  };

// Cylindrical Hankel functions.
template<typename _Tp>
  struct testcase_cyl_hankel_2
  {
    std::complex<_Tp> f0;
    _Tp nu;
    _Tp x;
    std::complex<_Tp> f;
  };

// Spherical Hankel functions.
template<typename _Tp>
  struct testcase_sph_hankel_1
  {
    std::complex<_Tp> f0;
    unsigned int n;
    _Tp x;
    std::complex<_Tp> f;
  };

// Spherical Hankel functions.
template<typename _Tp>
  struct testcase_sph_hankel_2
  {
    std::complex<_Tp> f0;
    unsigned int n;
    _Tp x;
    std::complex<_Tp> f;
  };

// Spherical Harmonic functions.
template<typename _Tp>
  struct testcase_sph_harmonic
  {
    std::complex<_Tp> f0;
    unsigned int l;
    int m;
    _Tp theta;
    _Tp phi;
    std::complex<_Tp> f;
  };

// Jacobi polynomials.
template<typename _Tp>
  struct testcase_jacobi
  {
    _Tp f0;
    unsigned int n;
    _Tp alpha;
    _Tp beta;
    _Tp x;
    _Tp f;
  };

// Polylogarithm functions.
template<typename _Tp>
  struct testcase_polylog
  {
    std::complex<_Tp> f0;
    _Tp s;
    std::complex<_Tp> w;
    std::complex<_Tp> f;
  };

// Clausen functions.
template<typename _Tp>
  struct testcase_clausen
  {
    std::complex<_Tp> f0;
    unsigned int m;
    std::complex<_Tp> w;
    std::complex<_Tp> f;
  };

// Dirichlet eta function.
template<typename _Tp>
  struct testcase_dirichlet_eta
  {
    _Tp f0;
    _Tp s;
    _Tp f;
  };

// Dirichlet beta function.
template<typename _Tp>
  struct testcase_dirichlet_beta
  {
    _Tp f0;
    _Tp s;
    _Tp f;
  };

// Dirichlet lambda function.
template<typename _Tp>
  struct testcase_dirichlet_lambda
  {
    _Tp f0;
    _Tp s;
    _Tp f;
  };

// Owens T functions.
template<typename _Tp>
  struct testcase_owens_t
  {
    _Tp f0;
    _Tp h;
    _Tp a;
    _Tp f;
  };

// Clausen Cl_2 function.
template<typename _Tp>
  struct testcase_clausen_c
  {
    _Tp f0;
    unsigned int m;
    _Tp w;
    _Tp f;
  };

// Exponential theta_1 functions.
template<typename _Tp>
  struct testcase_theta_1
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Exponential theta_2 functions.
template<typename _Tp>
  struct testcase_theta_2
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Exponential theta_3 functions.
template<typename _Tp>
  struct testcase_theta_3
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Exponential theta_4 functions.
template<typename _Tp>
  struct testcase_theta_4
  {
    _Tp f0;
    _Tp nu;
    _Tp x;
    _Tp f;
  };

// Elliptic nome.
template<typename _Tp>
  struct testcase_ellnome
  {
    _Tp f0;
    _Tp k;
    _Tp f;
  };

// Reperiodized sine function.
template<typename _Tp>
  struct testcase_sin_pi
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Reperiodized cosine function.
template<typename _Tp>
  struct testcase_cos_pi
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

// Reperiodized tangent function.
template<typename _Tp>
  struct testcase_tan_pi
  {
    _Tp f0;
    _Tp x;
    _Tp f;
  };

#endif // _GLIBCXX_SPECFUN_TESTCASE_H
