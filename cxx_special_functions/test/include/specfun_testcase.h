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
    Tp f0;
    unsigned int n;
    _Tpa alpha;
    Tp x;
    Tp f;
  };

// Associated Legendre functions.
template<typename Tp>
  struct testcase_assoc_legendre
  {
    Tp f0;
    unsigned int l;
    unsigned int m;
    Tp x;
    Tp f;
  };

// Beta function.
template<typename Tp>
  struct testcase_beta
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp f;
  };

// Complete elliptic integrals of the first kind.
template<typename Tp>
  struct testcase_comp_ellint_1
  {
    Tp f0;
    Tp k;
    Tp f;
  };

// Complete elliptic integrals of the second kind.
template<typename Tp>
  struct testcase_comp_ellint_2
  {
    Tp f0;
    Tp k;
    Tp f;
  };

// Complete elliptic integrals of the third kind.
template<typename Tp>
  struct testcase_comp_ellint_3
  {
    Tp f0;
    Tp k;
    Tp nu;
    Tp f;
  };

// Complete elliptic D integrals.
template<typename Tp>
  struct testcase_comp_ellint_d
  {
    Tp f0;
    Tp k;
    Tp f;
  };

// Confluent hypergeometric functions.
template<typename Tp>
  struct testcase_conf_hyperg
  {
    Tp f0;
    Tp a;
    Tp c;
    Tp x;
    Tp f;
  };

// Confluent hypergeometric functions.
template<typename Tp>
  struct testcase_conf_hyperg_lim
  {
    Tp f0;
    Tp c;
    Tp x;
    Tp f;
  };

// Regular Coulomb functions.
template<typename Tp>
  struct testcase_coulomb_f
  {
    Tp f0;
    Tp lambda;
    Tp eta;
    Tp rho;
    Tp f;
  };

// Regular Coulomb functions.
template<typename Tp>
  struct testcase_coulomb_g
  {
    Tp f0;
    Tp lambda;
    Tp eta;
    Tp rho;
    Tp f;
  };

// Generic cylindrical Bessel functions.
template<typename Tp>
  struct testcase_cyl_bessel
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Regular modified cylindrical Bessel functions.
template<typename Tp>
  struct testcase_cyl_bessel_i
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Scaled regular modified cylindrical Bessel functions.
template<typename Tp>
  struct testcase_cyl_bessel_i_scaled : public testcase_cyl_bessel_i<Tp>
  { };

// Cylindrical Bessel functions (of the first kind).
template<typename Tp>
  struct testcase_cyl_bessel_j
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Irregular modified cylindrical Bessel functions.
template<typename Tp>
  struct testcase_cyl_bessel_k
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Scaled irregular modified cylindrical Bessel functions.
template<typename Tp>
  struct testcase_cyl_bessel_k_scaled : public testcase_cyl_bessel_k<Tp>
  { };

// Cylindrical Neumann functions.
template<typename Tp>
  struct testcase_cyl_neumann
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Elliptic integrals of the first kind.
template<typename Tp>
  struct testcase_ellint_1
  {
    Tp f0;
    Tp k;
    Tp phi;
    Tp f;
  };

// Elliptic integrals of the second kind.
template<typename Tp>
  struct testcase_ellint_2
  {
    Tp f0;
    Tp k;
    Tp phi;
    Tp f;
  };

// Elliptic integrals of the third kind.
template<typename Tp>
  struct testcase_ellint_3
  {
    Tp f0;
    Tp k;
    Tp nu;
    Tp phi;
    Tp f;
  };

// Elliptic D integrals.
template<typename Tp>
  struct testcase_ellint_d
  {
    Tp f0;
    Tp k;
    Tp phi;
    Tp f;
  };

// Heuman lambda functions.
template<typename Tp>
  struct testcase_heuman_lambda
  {
    Tp f0;
    Tp k;
    Tp phi;
    Tp f;
  };

// Jacobi zeta functions.
template<typename Tp>
  struct testcase_jacobi_zeta
  {
    Tp f0;
    Tp k;
    Tp phi;
    Tp f;
  };

// Exponential integral.
template<typename Tp>
  struct testcase_expint
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Hermite polynomials
template<typename Tp>
  struct testcase_hermite
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Hypergeometric functions.
template<typename Tp>
  struct testcase_hyperg
  {
    Tp f0;
    Tp a;
    Tp b;
    Tp c;
    Tp x;
    Tp f;
  };

// Laguerre polynomials.
template<typename Tp>
  struct testcase_laguerre
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Legendre polynomials.
template<typename Tp>
  struct testcase_legendre
  {
    Tp f0;
    unsigned int l;
    Tp x;
    Tp f;
  };

// Riemann zeta function.
template<typename Tp>
  struct testcase_riemann_zeta
  {
    Tp f0;
    Tp s;
    Tp f;
  };

// Hurwitz zeta function.
template<typename Tp>
  struct testcase_hurwitz_zeta
  {
    Tp f0;
    Tp s;
    Tp a;
    Tp f;
  };

// Spherical Bessel functions.
template<typename Tp>
  struct testcase_sph_bessel
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Regular modified spherical Bessel functions.
template<typename Tp>
  struct testcase_sph_bessel_i
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Irregular modified spherical Bessel functions.
template<typename Tp>
  struct testcase_sph_bessel_k
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Spherical Legendre functions.
template<typename Tp>
  struct testcase_sph_legendre
  {
    Tp f0;
    unsigned int l;
    unsigned int m;
    Tp theta;
    Tp f;
  };

// Spherical Neumann functions.
template<typename Tp>
  struct testcase_sph_neumann
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Airy Ai functions.
template<typename Tp>
  struct testcase_airy_ai
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Airy Ai functions.
template<typename Tp>
  struct testcase_airy_ai_scaled
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Airy Bi functions.
template<typename Tp>
  struct testcase_airy_bi
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Airy Bi functions.
template<typename Tp>
  struct testcase_airy_bi_scaled
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Upper incomplete gamma functions.
template<typename Tp>
  struct testcase_tgamma
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Lower incomplete gamma functions.
template<typename Tp>
  struct testcase_tgamma_lower
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Regularized upper incomplete gamma functions.
template<typename Tp>
  struct testcase_gamma_q
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Regularized lower incomplete gamma functions.
template<typename Tp>
  struct testcase_gamma_p
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Dilogarithm functions.
template<typename Tp>
  struct testcase_dilog
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Gamma function.
template<typename Tp>
  struct testcase_gamma
  {
    Tp f0;
    Tp x;
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
    Tp f0;
    Tp x;
    Tp y;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rd
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp z;
    Tp f;
  };

template<typename Tp>
  struct testcase_comp_ellint_rf
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rf
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp z;
    Tp f;
  };

template<typename Tp>
  struct testcase_comp_ellint_rg
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rg
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp z;
    Tp f;
  };

template<typename Tp>
  struct testcase_ellint_rj
  {
    Tp f0;
    Tp x;
    Tp y;
    Tp z;
    Tp p;
    Tp f;
  };

template<typename Tp>
  struct testcase_digamma
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_polygamma
  {
    Tp f0;
    unsigned int m;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinint
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_cosint
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinhint
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_coshint
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_dawson
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_jacobi_sn
  {
    Tp f0;
    Tp k;
    Tp u;
    Tp f;
  };

template<typename Tp>
  struct testcase_jacobi_cn
  {
    Tp f0;
    Tp k;
    Tp u;
    Tp f;
  };

template<typename Tp>
  struct testcase_jacobi_dn
  {
    Tp f0;
    Tp k;
    Tp u;
    Tp f;
  };

template<typename Tp>
  struct testcase_expint_en
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_fresnel_c
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_fresnel_s
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinc
  {
    Tp f0;
    Tp x;
    Tp f;
  };

template<typename Tp>
  struct testcase_sinc_pi
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Log rising factorials.
template<typename Tp>
  struct testcase_lrising_factorial
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Rising factorials.
template<typename Tp>
  struct testcase_rising_factorial
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Log falling factorials.
template<typename Tp>
  struct testcase_lfalling_factorial
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Falling factorials.
template<typename Tp>
  struct testcase_falling_factorial
  {
    Tp f0;
    Tp a;
    Tp x;
    Tp f;
  };

// Legendre functions of the second kind.
template<typename Tp>
  struct testcase_legendre_q
  {
    Tp f0;
    unsigned int l;
    Tp x;
    Tp f;
  };

// Factorial.
template<typename Tp>
  struct testcase_factorial
  {
    Tp f0;
    unsigned int n;
    Tp f;
  };

// Log factorial.
template<typename Tp>
  struct testcase_lfactorial
  {
    Tp f0;
    unsigned int n;
    Tp f;
  };

// Double factorial.
template<typename Tp>
  struct testcase_double_factorial
  {
    Tp f0;
    unsigned int n;
    Tp f;
  };

// Log double factorial.
template<typename Tp>
  struct testcase_ldouble_factorial
  {
    Tp f0;
    unsigned int n;
    Tp f;
  };

// Binomial coefficient.
template<typename Tp>
  struct testcase_binomial
  {
    Tp f0;
    unsigned int n;
    unsigned int k;
    Tp f;
  };

// Log binomial coefficient.
template<typename Tp>
  struct testcase_lbinomial
  {
    Tp f0;
    unsigned int n;
    unsigned int k;
    Tp f;
  };

// Gegenbauer polynomials.
template<typename Tp>
  struct testcase_gegenbauer
  {
    Tp f0;
    unsigned int n;
    Tp lambda;
    Tp x;
    Tp f;
  };

// Chebyshev polynomials of the first kind.
template<typename Tp>
  struct testcase_chebyshev_t
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Chebyshev polynomials of the second kind.
template<typename Tp>
  struct testcase_chebyshev_u
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Chebyshev polynomials of the third kind.
template<typename Tp>
  struct testcase_chebyshev_v
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Chebyshev polynomials of the fourth kind.
template<typename Tp>
  struct testcase_chebyshev_w
  {
    Tp f0;
    unsigned int n;
    Tp x;
    Tp f;
  };

// Zernike polynomials.
template<typename Tp>
  struct testcase_zernike
  {
    Tp f0;
    unsigned int n;
    int m;
    Tp rho;
    Tp phi;
    Tp f;
  };

// Radial polynomials.
template<typename Tp>
  struct testcase_radpoly
  {
    Tp f0;
    unsigned int n;
    int m;
    Tp rho;
    Tp f;
  };

// Incomplete beta function.
template<typename Tp>
  struct testcase_ibeta
  {
    Tp f0;
    Tp a;
    Tp b;
    Tp x;
    Tp f;
  };

// Complementary beta function.
template<typename Tp>
  struct testcase_ibetac
  {
    Tp f0;
    Tp a;
    Tp b;
    Tp x;
    Tp f;
  };

// Cylindrical Hankel functions.
template<typename Tp>
  struct testcase_cyl_hankel_1
  {
    std::complex<Tp> f0;
    Tp nu;
    Tp x;
    std::complex<Tp> f;
  };

// Cylindrical Hankel functions.
template<typename Tp>
  struct testcase_cyl_hankel_2
  {
    std::complex<Tp> f0;
    Tp nu;
    Tp x;
    std::complex<Tp> f;
  };

// Spherical Hankel functions.
template<typename Tp>
  struct testcase_sph_hankel_1
  {
    std::complex<Tp> f0;
    unsigned int n;
    Tp x;
    std::complex<Tp> f;
  };

// Spherical Hankel functions.
template<typename Tp>
  struct testcase_sph_hankel_2
  {
    std::complex<Tp> f0;
    unsigned int n;
    Tp x;
    std::complex<Tp> f;
  };

// Spherical Harmonic functions.
template<typename Tp>
  struct testcase_sph_harmonic
  {
    std::complex<Tp> f0;
    unsigned int l;
    int m;
    Tp theta;
    Tp phi;
    std::complex<Tp> f;
  };

// Jacobi polynomials.
template<typename Tp>
  struct testcase_jacobi
  {
    Tp f0;
    unsigned int n;
    Tp alpha;
    Tp beta;
    Tp x;
    Tp f;
  };

// Polylogarithm functions.
template<typename Tp>
  struct testcase_polylog
  {
    std::complex<Tp> f0;
    Tp s;
    std::complex<Tp> w;
    std::complex<Tp> f;
  };

// Clausen functions.
template<typename Tp>
  struct testcase_clausen
  {
    std::complex<Tp> f0;
    unsigned int m;
    std::complex<Tp> w;
    std::complex<Tp> f;
  };

// Dirichlet eta function.
template<typename Tp>
  struct testcase_dirichlet_eta
  {
    Tp f0;
    Tp s;
    Tp f;
  };

// Dirichlet beta function.
template<typename Tp>
  struct testcase_dirichlet_beta
  {
    Tp f0;
    Tp s;
    Tp f;
  };

// Dirichlet lambda function.
template<typename Tp>
  struct testcase_dirichlet_lambda
  {
    Tp f0;
    Tp s;
    Tp f;
  };

// Owens T functions.
template<typename Tp>
  struct testcase_owens_t
  {
    Tp f0;
    Tp h;
    Tp a;
    Tp f;
  };

// Clausen Cl_m function.
template<typename Tp>
  struct testcase_clausen_cl
  {
    Tp f0;
    unsigned int m;
    Tp w;
    Tp f;
  };

// Exponential theta_1 functions.
template<typename Tp>
  struct testcase_theta_1
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Exponential theta_2 functions.
template<typename Tp>
  struct testcase_theta_2
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Exponential theta_3 functions.
template<typename Tp>
  struct testcase_theta_3
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Exponential theta_4 functions.
template<typename Tp>
  struct testcase_theta_4
  {
    Tp f0;
    Tp nu;
    Tp x;
    Tp f;
  };

// Elliptic nome.
template<typename Tp>
  struct testcase_ellnome
  {
    Tp f0;
    Tp k;
    Tp f;
  };

// Reperiodized sine function.
template<typename Tp>
  struct testcase_sin_pi
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Reperiodized cosine function.
template<typename Tp>
  struct testcase_cos_pi
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Reperiodized tangent function.
template<typename Tp>
  struct testcase_tan_pi
  {
    Tp f0;
    Tp x;
    Tp f;
  };

// Log gamma functions.
template<typename Tp>
  struct testcase_lgamma
  {
    Tp f0;
    Tp a;
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

//
template<typename Tp>
  struct testcase_boltzmann_p
  {
    Tp f0;
    Tp mu;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_boltzmann_pdf
  {
    Tp f0;
    Tp mu;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_laplace_p
  {
    Tp f0;
    Tp mu;
    Tp sigma;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_laplace_pdf
  {
    Tp f0;
    Tp mu;
    Tp sigma;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_maxwell_p
  {
    Tp f0;
    Tp mu;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_maxwell_pdf
  {
    Tp f0;
    Tp mu;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_normal_p
  {
    Tp f0;
    Tp mu;
    Tp sigma;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_normal_pdf
  {
    Tp f0;
    Tp mu;
    Tp sigma;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_rayleigh_p
  {
    Tp f0;
    Tp mu;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_rayleigh_pdf
  {
    Tp f0;
    Tp mu;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_lognormal_p
  {
    Tp f0;
    Tp mu;
    Tp sigma;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_lognormal_pdf
  {
    Tp f0;
    Tp mu;
    Tp sigma;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_logistic_p
  {
    Tp f0;
    Tp mu;
    Tp s;
    Tp x;
    Tp f;
  };

//
template<typename Tp>
  struct testcase_logistic_pdf
  {
    Tp f0;
    Tp mu;
    Tp s;
    Tp x;
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

// Euler numbers.
template<typename Tp>
  struct testcase_euler
  {
    Tp f0;
    unsigned int n;
    Tp f;
  };

// Eulerian numbers of the first kind.
template<typename Tp>
  struct testcase_eulerian_1
  {
    Tp f0;
    unsigned int n;
    unsigned int m;
    Tp f;
  };

// Eulerian numbers of the second kind.
template<typename Tp>
  struct testcase_eulerian_2
  {
    Tp f0;
    unsigned int n;
    unsigned int m;
    Tp f;
  };

// Stirling numbers of the first kind.
template<typename Tp>
  struct testcase_stirling_1
  {
    Tp f0;
    unsigned int n;
    unsigned int m;
    Tp f;
  };

// Stirling numbers of the second kind.
template<typename Tp>
  struct testcase_stirling_2
  {
    Tp f0;
    unsigned int n;
    unsigned int m;
    Tp f;
  };

#endif // SPECFUN_TESTCASE_H
