// Copyright (C) 2015-2018 Free Software Foundation, Inc.
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
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

// specfun_testcase.h

//
// These are little PODs for special function inputs and
// expexted results for the testsuite.
//

#ifndef SPECFUN_TESTCASE_H
#define SPECFUN_TESTCASE_H

#include <complex>

// Associated Laguerre polynomials.
template<typename _Tpa, typename Tp>
  struct testcase_assoc_laguerre
  {
    unsigned int n;
    _Tpa alpha;
    Tp x;
    Tp f0;
    Tp f;
  };

// Associated Legendre functions.
template<typename Tp>
  struct testcase_assoc_legendre
  {
    unsigned int l;
    unsigned int m;
    Tp x;
    Tp f0;
    Tp f;
  };

// Beta function.
template<typename Tp>
  struct testcase_beta
  {
    Tp x;
    Tp y;
    Tp f0;
    Tp f;
  };

// Complete elliptic integrals of the first kind.
template<typename Tp>
  struct testcase_comp_ellint_1
  {
    Tp k;
    Tp f0;
    Tp f;
  };

// Complete elliptic integrals of the second kind.
template<typename Tp>
  struct testcase_comp_ellint_2
  {
    Tp k;
    Tp f0;
    Tp f;
  };

// Complete elliptic integrals of the third kind.
template<typename Tp>
  struct testcase_comp_ellint_3
  {
    Tp k;
    Tp nu;
    Tp f0;
    Tp f;
  };

// Complete elliptic D integrals.
template<typename Tp>
  struct testcase_comp_ellint_d
  {
    Tp k;
    Tp f0;
    Tp f;
  };

// Confluent hypergeometric functions.
template<typename Tp>
  struct testcase_conf_hyperg
  {
    Tp a;
    Tp c;
    Tp x;
    Tp f0;
    Tp f;
  };

// Confluent hypergeometric functions.
template<typename Tp>
  struct testcase_conf_hyperg_lim
  {
    Tp c;
    Tp x;
    Tp f0;
    Tp f;
  };

// Generic cylindrical Bessel functions.
template<typename Tp>
  struct testcase_cyl_bessel
  {
    Tp nu;
    Tp x;
    Tp f0;
    Tp f;
  };

// Cylindrical Bessel functions (of the first kind).
template<typename Tp>
  struct testcase_cyl_bessel
  {
    Tp nu;
    Tp x;
    Tp f0;
    Tp f;
  };

// Elliptic integrals of the first kind.
template<typename Tp>
  struct testcase_ellint_1
  {
    Tp k;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Elliptic integrals of the second kind.
template<typename Tp>
  struct testcase_ellint_2
  {
    Tp k;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Elliptic integrals of the third kind.
template<typename Tp>
  struct testcase_ellint_3
  {
    Tp k;
    Tp nu;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Elliptic D integrals.
template<typename Tp>
  struct testcase_ellint_d
  {
    Tp k;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Heuman lambda functions.
template<typename Tp>
  struct testcase_heuman_lambda
  {
    Tp k;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Jacobi zeta functions.
template<typename Tp>
  struct testcase_jacobi_zeta
  {
    Tp k;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Exponential integral.
template<typename Tp>
  struct testcase_expint
  {
    Tp x;
    Tp f0;
    Tp f;
  };

// Hermite polynomials
template<typename Tp>
  struct testcase_hermite
  {
    unsigned int n;
    Tp x;
    Tp f0;
    Tp f;
  };

// Hypergeometric functions.
template<typename Tp>
  struct testcase_hyperg
  {
    Tp a;
    Tp b;
    Tp c;
    Tp x;
    Tp f0;
    Tp f;
  };

// Laguerre polynomials.
template<typename Tp>
  struct testcase_laguerre
  {
    unsigned int n;
    Tp x;
    Tp f0;
    Tp f;
  };

// Legendre polynomials.
template<typename Tp>
  struct testcase_legendre
  {
    unsigned int l;
    Tp x;
    Tp f0;
    Tp f;
  };

// Riemann zeta function.
template<typename Tp>
  struct testcase_riemann_zeta
  {
    Tp s;
    Tp f0;
    Tp f;
  };

// Hurwitz zeta function.
template<typename Tp>
  struct testcase_hurwitz_zeta
  {
    Tp s;
    Tp a;
    Tp f0;
    Tp f;
  };

// Spherical Bessel functions.
template<typename Tp>
  struct testcase_sph_bessel
  {
    unsigned int n;
    Tp x;
    Tp f0;
    Tp f;
  };

// Spherical Legendre functions.
template<typename Tp>
  struct testcase_sph_legendre
  {
    unsigned int l;
    unsigned int m;
    Tp theta;
    Tp f0;
    Tp f;
  };

// Airy functions.
template<typename Tp>
  struct testcase_airy
  {
    Tp x;
    Tp f0;
    Tp f;
  };

// Incomplete gamma functions.
template<typename Tp>
  struct testcase_gamma_inc
  {
    Tp a;
    Tp x;
    Tp f0;
    Tp f;
  };

// Regularized incomplete gamma functions.
template<typename Tp>
  struct testcase_xgamma
  {
    Tp a;
    Tp x;
    Tp f0;
    Tp f;
  };

// Dilogarithm functions.
template<typename Tp>
  struct testcase_dilog
  {
    Tp x;
    Tp f0;
    Tp f;
  };

// Gamma function.
template<typename Tp>
  struct testcase_gamma
  {
    Tp x;
    Tp f0;
    Tp f;
  };

// Reciprocal Gamma function.
template<typename Tp>
  struct testcase_gamma_reciprocal
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rc
  {
    Tp x;
    Tp y;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rd
  {
    Tp x;
    Tp y;
    Tp z;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_comp_ellint_rf
  {
    Tp x;
    Tp y;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rf
  {
    Tp x;
    Tp y;
    Tp z;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_comp_ellint_rg
  {
    Tp x;
    Tp y;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rg
  {
    Tp x;
    Tp y;
    Tp z;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rj
  {
    Tp x;
    Tp y;
    Tp z;
    Tp p;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_digamma
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_polygamma
  {
    unsigned int m;
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinint
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_cosint
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinhint
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_coshint
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_dawson
  {
    Tp x;
    Tp f0;
    Tp f;
  };

// Jacobi elliptic functions.
template<typename Tp>
  struct testcase_jacobi_xn
  {
    Tp k;
    Tp u;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_expint_en
  {
    unsigned int n;
    Tp x;
    Tp f0;
    Tp f;
  };

// Fresnel functions
template<typename Tp>
  struct testcase_fresnel
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinc
  {
    Tp x;
    Tp f0;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinc_pi
  {
    Tp x;
    Tp f0;
    Tp f;
  };

// Log rising factorials.
template<typename Tp>
  struct testcase_lrising_factorial
  {
    Tp a;
    Tp x;
    Tp f0;
    Tp f;
  };

// Rising factorials.
template<typename Tp>
  struct testcase_rising_factorial
  {
    Tp a;
    Tp x;
    Tp f0;
    Tp f;
  };

// Log falling factorials.
template<typename Tp>
  struct testcase_lfalling_factorial
  {
    Tp a;
    Tp x;
    Tp f0;
    Tp f;
  };

// Falling factorials.
template<typename Tp>
  struct testcase_falling_factorial
  {
    Tp a;
    Tp x;
    Tp f0;
    Tp f;
  };

// Factorial type functions.
template<typename Tp>
  struct testcase_factorial
  {
    unsigned int n;
    Tp f0;
    Tp f;
  };

// Binomial coefficients.
template<typename Tp>
  struct testcase_binomial
  {
    unsigned int n;
    unsigned int k;
    Tp f0;
    Tp f;
  };

// Gegenbauer polynomials.
template<typename Tp>
  struct testcase_gegenbauer
  {
    unsigned int n;
    Tp lambda;
    Tp x;
    Tp f0;
    Tp f;
  };

// Chebyshev polynomials.
template<typename Tp>
  struct testcase_chebyshev
  {
    unsigned int n;
    Tp x;
    Tp f0;
    Tp f;
  };

// Zernike polynomials.
template<typename Tp>
  struct testcase_zernike
  {
    unsigned int n;
    int m;
    Tp rho;
    Tp phi;
    Tp f0;
    Tp f;
  };

// Radial polynomials.
template<typename Tp>
  struct testcase_radpoly
  {
    unsigned int n;
    int m;
    Tp rho;
    Tp f0;
    Tp f;
  };

// Incomplete beta functions.
template<typename Tp>
  struct testcase_ibeta
  {
    Tp a;
    Tp b;
    Tp x;
    Tp f0;
    Tp f;
  };

// Cylindrical Hankel functions.
template<typename Tp>
  struct testcase_cyl_hankel
  {
    Tp nu;
    Tp x;
    std::complex<Tp> f0;
    std::complex<Tp> f;
  };

// Spherical Hankel functions.
template<typename Tp>
  struct testcase_sph_hankel
  {
    unsigned int n;
    Tp x;
    std::complex<Tp> f0;
    std::complex<Tp> f;
  };

// Spherical Harmonic functions.
template<typename Tp>
  struct testcase_sph_harmonic
  {
    unsigned int l;
    int m;
    Tp theta;
    Tp phi;
    std::complex<Tp> f0;
    std::complex<Tp> f;
  };

// Jacobi polynomials.
template<typename Tp>
  struct testcase_jacobi
  {
    unsigned int n;
    Tp alpha;
    Tp beta;
    Tp x;
    Tp f0;
    Tp f;
  };

// Polylogarithm functions.
template<typename Tp>
  struct testcase_polylog
  {
    Tp s;
    std::complex<Tp> w;
    std::complex<Tp> f0;
    std::complex<Tp> f;
  };

// Clausen functions.
template<typename Tp>
  struct testcase_clausen
  {
    unsigned int m;
    std::complex<Tp> w;
    std::complex<Tp> f0;
    std::complex<Tp> f;
  };

// Dirichlet eta function.
template<typename Tp>
  struct testcase_dirichlet_eta
  {
    Tp s;
    Tp f0;
    Tp f;
  };

// Dirichlet beta function.
template<typename Tp>
  struct testcase_dirichlet_beta
  {
    Tp s;
    Tp f0;
    Tp f;
  };

// Dirichlet lambda function.
template<typename Tp>
  struct testcase_dirichlet_lambda
  {
    Tp s;
    Tp f0;
    Tp f;
  };

// Owens T functions.
template<typename Tp>
  struct testcase_owens_t
  {
    Tp h;
    Tp a;
    Tp f0;
    Tp f;
  };

// Clausen Cl_m function.
template<typename Tp>
  struct testcase_clausen_cl
  {
    unsigned int m;
    Tp w;
    Tp f0;
    Tp f;
  };

// Exponential theta_N functions.
template<typename Tp>
  struct testcase_theta
  {
    Tp nu;
    Tp x;
    Tp f0;
    Tp f;
  };

// Elliptic nome.
template<typename Tp>
  struct testcase_ellnome
  {
    Tp k;
    Tp f0;
    Tp f;
  };

// Debye functions.
template<typename Tp>
  struct testcase_debye
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Bernoulli numbers.
template<typename Tp>
  struct testcase_bernoulli
  {
    Tp f0;
    unsigned int n;
    Tp f;
  };

// Euler numbers.
template<typename Tp>
  struct testcase_euler
  {
    unsigned int n;
    Tp f0;
    Tp f;
  };

// Eulerian numbers.
template<typename Tp>
  struct testcase_eulerian
  {
    unsigned int n;
    unsigned int m;
    Tp f0;
    Tp f;
  };

// Stirling numbers.
template<typename Tp>
  struct testcase_stirling
  {
    unsigned int n;
    unsigned int m;
    Tp f0;
    Tp f;
  };

#endif // SPECFUN_TESTCASE_H
