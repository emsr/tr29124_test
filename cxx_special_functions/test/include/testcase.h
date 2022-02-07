// 2007-02-04  Edward Smith-Rowland <3dw4rd@verizon.net>
//
// Copyright (C) 2007-2019 Free Software Foundation, Inc.
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

//  testcase.h

//
//  These are little PODs for special function inputs and
//  expected results for the testsuite.
//

//  5.2.1.1  Associated Laguerre polynomials.
template <typename Tp>
struct testcase_assoc_laguerre
{
  Tp f0;
  unsigned int n;
  unsigned int m;
  Tp x;
  Tp f;
};

//  5.2.1.2  Associated Legendre functions.
template <typename Tp>
struct testcase_assoc_legendre
{
  Tp f0;
  unsigned int l;
  unsigned int m;
  Tp x;
  Tp f;
};

//  5.2.1.3  Beta function.
template <typename Tp>
struct testcase_beta
{
  Tp f0;
  Tp x;
  Tp y;
  Tp f;
};

//  5.2.1.4  Complete elliptic integrals of the first kind.
template <typename Tp>
struct testcase_comp_ellint_1
{
  Tp f0;
  Tp k;
  Tp f;
};

//  5.2.1.5  Complete elliptic integrals of the second kind.
template <typename Tp>
struct testcase_comp_ellint_2
{
  Tp f0;
  Tp k;
  Tp f;
};

//  5.2.1.6  Complete elliptic integrals of the third kind.
template <typename Tp>
struct testcase_comp_ellint_3
{
  Tp f0;
  Tp k;
  Tp nu;
  Tp f;
};

//  5.2.1.7  Confluent hypergeometric functions.
template <typename Tp>
struct testcase_conf_hyperg
{
  Tp f0;
  Tp a;
  Tp c;
  Tp x;
  Tp f;
};

//  5.2.1.8  Regular modified cylindrical Bessel functions.
template <typename Tp>
struct testcase_cyl_bessel_i
{
  Tp f0;
  Tp nu;
  Tp x;
  Tp f;
};

//  5.2.1.9  Cylindrical Bessel functions (of the first kind).
template <typename Tp>
struct testcase_cyl_bessel_j
{
  Tp f0;
  Tp nu;
  Tp x;
  Tp f;
};

//  5.2.1.10  Irregular modified cylindrical Bessel functions.
template <typename Tp>
struct testcase_cyl_bessel_k
{
  Tp f0;
  Tp nu;
  Tp x;
  Tp f;
};

//  5.2.1.11  Cylindrical Neumann functions.
template <typename Tp>
struct testcase_cyl_neumann
{
  Tp f0;
  Tp nu;
  Tp x;
  Tp f;
};

//  5.2.1.12  Elliptic integrals of the first kind.
template <typename Tp>
struct testcase_ellint_1
{
  Tp f0;
  Tp k;
  Tp phi;
  Tp f;
};

//  5.2.1.13  Elliptic integrals of the second kind.
template <typename Tp>
struct testcase_ellint_2
{
  Tp f0;
  Tp k;
  Tp phi;
  Tp f;
};

//  5.2.1.14  Elliptic integrals of the third kind.
template <typename Tp>
struct testcase_ellint_3
{
  Tp f0;
  Tp k;
  Tp nu;
  Tp phi;
  Tp f;
};

//  5.2.1.15  Exponential integral.
template <typename Tp>
struct testcase_expint
{
  Tp f0;
  Tp x;
  Tp f;
};

//  5.2.1.16  Hermite polynomials
template <typename Tp>
struct testcase_hermite
{
  Tp f0;
  unsigned int n;
  Tp x;
  Tp f;
};

//  5.2.1.17  Hypergeometric functions.
template <typename Tp>
struct testcase_hyperg
{
  Tp f0;
  Tp a;
  Tp b;
  Tp c;
  Tp x;
  Tp f;
};

//  5.2.1.18  Laguerre polynomials.
template <typename Tp>
struct testcase_laguerre
{
  Tp f0;
  unsigned int n;
  Tp x;
  Tp f;
};

//  5.2.1.19  Legendre polynomials.
template <typename Tp>
struct testcase_legendre
{
  Tp f0;
  unsigned int l;
  Tp x;
  Tp f;
};

//  5.2.1.20  Riemann zeta function.
template <typename Tp>
struct testcase_riemann_zeta
{
  Tp f0;
  Tp x;
  Tp f;
};

//  5.2.1.21  Spherical Bessel functions.
template <typename Tp>
struct testcase_sph_bessel
{
  Tp f0;
  unsigned int n;
  Tp x;
  Tp f;
};

//  5.2.1.22  Spherical Legendre functions.
template <typename Tp>
struct testcase_sph_legendre
{
  Tp f0;
  unsigned int l;
  unsigned int m;
  Tp theta;
  Tp f;
};

//  5.2.1.23  Spherical Neumann functions.
template <typename Tp>
struct testcase_sph_neumann
{
  Tp f0;
  unsigned int n;
  Tp x;
  Tp f;
};
