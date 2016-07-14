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
    unsigned int n;
    unsigned int m;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Associated Legendre functions.
template<typename _Tp>
  struct testcase_assoc_legendre
  {
    unsigned int l;
    unsigned int m;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Beta function.
template<typename _Tp>
  struct testcase_beta
  {
    _Tp x;
    _Tp y;
    _Tp f0;
    _Tp f;
  };

// Complete elliptic integrals of the first kind.
template<typename _Tp>
  struct testcase_comp_ellint_1
  {
    _Tp k;
    _Tp f0;
    _Tp f;
  };

// Complete elliptic integrals of the second kind.
template<typename _Tp>
  struct testcase_comp_ellint_2
  {
    _Tp k;
    _Tp f0;
    _Tp f;
  };

// Complete elliptic integrals of the third kind.
template<typename _Tp>
  struct testcase_comp_ellint_3
  {
    _Tp k;
    _Tp nu;
    _Tp f0;
    _Tp f;
  };

// Complete elliptic D integrals.
template<typename _Tp>
  struct testcase_comp_ellint_d
  {
    _Tp k;
    _Tp f0;
    _Tp f;
  };

// Confluent hypergeometric functions.
template<typename _Tp>
  struct testcase_conf_hyperg
  {
    _Tp a;
    _Tp c;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Confluent hypergeometric functions.
template<typename _Tp>
  struct testcase_conf_hyperg_lim
  {
    _Tp c;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Generic cylindrical Bessel functions.
template<typename _Tp>
  struct testcase_cyl_bessel
  {
    _Tp nu;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Cylindrical Bessel functions (of the first kind).
template<typename _Tp>
  struct testcase_cyl_bessel
  {
    _Tp nu;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Elliptic integrals of the first kind.
template<typename _Tp>
  struct testcase_ellint_1
  {
    _Tp k;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Elliptic integrals of the second kind.
template<typename _Tp>
  struct testcase_ellint_2
  {
    _Tp k;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Elliptic integrals of the third kind.
template<typename _Tp>
  struct testcase_ellint_3
  {
    _Tp k;
    _Tp nu;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Elliptic D integrals.
template<typename _Tp>
  struct testcase_ellint_d
  {
    _Tp k;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Heuman lambda functions.
template<typename _Tp>
  struct testcase_heuman_lambda
  {
    _Tp k;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Jacobi zeta functions.
template<typename _Tp>
  struct testcase_jacobi_zeta
  {
    _Tp k;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Exponential integral.
template<typename _Tp>
  struct testcase_expint
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Hermite polynomials
template<typename _Tp>
  struct testcase_hermite
  {
    unsigned int n;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Hypergeometric functions.
template<typename _Tp>
  struct testcase_hyperg
  {
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Laguerre polynomials.
template<typename _Tp>
  struct testcase_laguerre
  {
    unsigned int n;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Legendre polynomials.
template<typename _Tp>
  struct testcase_legendre
  {
    unsigned int l;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Riemann zeta function.
template<typename _Tp>
  struct testcase_riemann_zeta
  {
    _Tp s;
    _Tp f0;
    _Tp f;
  };

// Hurwitz zeta function.
template<typename _Tp>
  struct testcase_hurwitz_zeta
  {
    _Tp s;
    _Tp a;
    _Tp f0;
    _Tp f;
  };

// Spherical Bessel functions.
template<typename _Tp>
  struct testcase_sph_bessel
  {
    unsigned int n;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Spherical Legendre functions.
template<typename _Tp>
  struct testcase_sph_legendre
  {
    unsigned int l;
    unsigned int m;
    _Tp theta;
    _Tp f0;
    _Tp f;
  };

// Airy functions.
template<typename _Tp>
  struct testcase_airy
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Incomplete gamma functions.
template<typename _Tp>
  struct testcase_gamma_inc
  {
    _Tp a;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Regularized incomplete gamma functions.
template<typename _Tp>
  struct testcase_xgamma
  {
    _Tp a;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Dilogarithm functions.
template<typename _Tp>
  struct testcase_dilog
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Digamma functions.
template<typename _Tp>
  struct testcase_gamma
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rc
  {
    _Tp x;
    _Tp y;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rd
  {
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_comp_ellint_rf
  {
    _Tp x;
    _Tp y;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rf
  {
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_comp_ellint_rg
  {
    _Tp x;
    _Tp y;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rg
  {
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_ellint_rj
  {
    _Tp x;
    _Tp y;
    _Tp z;
    _Tp p;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_psi
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinint
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_cosint
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinhint
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_coshint
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_dawson
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Jacobi elliptic functions.
template<typename _Tp>
  struct testcase_jacobi_xn
  {
    _Tp k;
    _Tp u;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_expint_en
  {
    unsigned int n;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Fresnel functions
template<typename _Tp>
  struct testcase_fresnel
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinc
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

template<typename _Tp>
  struct testcase_sinc_pi
  {
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Log upper Pochhammer symbol.
template<typename _Tp>
  struct testcase_lpochhammer_u
  {
    _Tp a;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Upper Pochhammer symbols.
template<typename _Tp>
  struct testcase_pochhammer_u
  {
    _Tp a;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Log lower Pochhammer symbol.
template<typename _Tp>
  struct testcase_lpochhammer_l
  {
    _Tp a;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Lower Pochhammer symbol.
template<typename _Tp>
  struct testcase_pochhammer_l
  {
    _Tp a;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Factorial type functions.
template<typename _Tp>
  struct testcase_factorial
  {
    unsigned int n;
    _Tp f0;
    _Tp f;
  };

// Binomial coefficients.
template<typename _Tp>
  struct testcase_bincoef
  {
    unsigned int n;
    unsigned int k;
    _Tp f0;
    _Tp f;
  };

// Gegenbauer polynomials.
template<typename _Tp>
  struct testcase_gegenbauer
  {
    unsigned int n;
    _Tp alpha;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Chebyshev polynomials.
template<typename _Tp>
  struct testcase_chebyshev
  {
    unsigned int n;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Zernike polynomials.
template<typename _Tp>
  struct testcase_zernike
  {
    unsigned int n;
    int m;
    _Tp rho;
    _Tp phi;
    _Tp f0;
    _Tp f;
  };

// Radial polynomials.
template<typename _Tp>
  struct testcase_radpoly
  {
    unsigned int n;
    int m;
    _Tp rho;
    _Tp f0;
    _Tp f;
  };

// Incomplete beta functions.
template<typename _Tp>
  struct testcase_ibeta
  {
    _Tp a;
    _Tp b;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Cylindrical Hankel functions.
template<typename _Tp>
  struct testcase_cyl_hankel
  {
    _Tp nu;
    _Tp x;
    std::complex<_Tp> f0;
    std::complex<_Tp> f;
  };

// Spherical Hankel functions.
template<typename _Tp>
  struct testcase_sph_hankel
  {
    unsigned int n;
    _Tp x;
    std::complex<_Tp> f0;
    std::complex<_Tp> f;
  };

// Spherical Harmonic functions.
template<typename _Tp>
  struct testcase_sph_harmonic
  {
    unsigned int l;
    int m;
    _Tp theta;
    _Tp phi;
    std::complex<_Tp> f0;
    std::complex<_Tp> f;
  };

// Jacobi polynomials.
template<typename _Tp>
  struct testcase_jacobi
  {
    unsigned int n;
    _Tp alpha;
    _Tp beta;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Polylogarithm functions.
template<typename _Tp>
  struct testcase_polylog
  {
    _Tp s;
    std::complex<_Tp> w;
    std::complex<_Tp> f0;
    std::complex<_Tp> f;
  };

// Clausen functions.
template<typename _Tp>
  struct testcase_clausen
  {
    unsigned int m;
    std::complex<_Tp> w;
    std::complex<_Tp> f0;
    std::complex<_Tp> f;
  };

// Dirichlet eta function.
template<typename _Tp>
  struct testcase_dirichlet_eta
  {
    _Tp s;
    _Tp f0;
    _Tp f;
  };

// Dirichlet beta function.
template<typename _Tp>
  struct testcase_dirichlet_beta
  {
    _Tp s;
    _Tp f0;
    _Tp f;
  };

// Dirichlet lambda function.
template<typename _Tp>
  struct testcase_dirichlet_lambda
  {
    _Tp s;
    _Tp f0;
    _Tp f;
  };

// Owens T functions.
template<typename _Tp>
  struct testcase_owens_t
  {
    _Tp h;
    _Tp a;
    _Tp f0;
    _Tp f;
  };

// Clausen Cl_2 function.
template<typename _Tp>
  struct testcase_clausen_c
  {
    unsigned int m;
    _Tp w;
    _Tp f0;
    _Tp f;
  };

// Exponential theta_N functions.
template<typename _Tp>
  struct testcase_theta
  {
    _Tp nu;
    _Tp x;
    _Tp f0;
    _Tp f;
  };

// Elliptic nome.
template<typename _Tp>
  struct testcase_ellnome
  {
    _Tp k;
    _Tp f0;
    _Tp f;
  };

#endif // _GLIBCXX_SPECFUN_TESTCASE_H
