// Special functions -*- C++ -*-

// Copyright (C) 2006-2017 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

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

/** @file bits/specfun_state.h
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SPECFUN_STATE_H
#define _GLIBCXX_BITS_SPECFUN_STATE_H 1

#pragma GCC system_header

#include <cmath>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A struct to store a cosine and a sine value.
   * A return for sincos-type functions.
   */
  template<typename _Tp>
    struct __quadrature_point_t
    {
      _Tp __zero;
      _Tp __weight;

      __quadrature_point_t() = default;

      __quadrature_point_t(_Tp __z, _Tp __w)
      : __zero(__z),
	__weight(__w)
      { }
    };

  /**
   * A struct to store a cosine and a sine value.
   * A return for sincos-type functions.
   */
  template<typename _Tp>
    struct __sincos_t
    {
      _Tp __sin_v;
      _Tp __cos_v;
    };

  // Slots for Jacobi elliptic function tuple.
  template<typename _Tp>
    struct __jacobi_t
    {
      /// Jacobi sine amplitude value.
      _Tp __sn_value;

      /// Jacobi cosine amplitude value.
      _Tp __cn_value;

      /// Jacobi delta amplitude value.
      _Tp __dn_value;

      _Tp __am() const
      { return std::asin(__sn_value); }

      _Tp __ns() const
      { return _Tp{1} / __sn_value; }

      _Tp __nc() const
      { return _Tp{1} / __cn_value; }

      _Tp __nd() const
      { return _Tp{1} / __dn_value; }

      _Tp __sc() const
      { return __sn_value / __cn_value; }

      _Tp __sd() const
      { return __sn_value / __dn_value; }

      _Tp __cd() const
      { return __cn_value / __dn_value; }

      _Tp __cs() const
      { return __cn_value / __sn_value; }

      _Tp __ds() const
      { return __dn_value / __sn_value; }

      _Tp __dc() const
      { return __dn_value / __cn_value; }
    };

  template<typename _Tx, typename _Tp>
    struct __airy_t
    {
      /// The argument of the Airy fuctions.
      _Tx __x_arg;

      /// The value of the Airy function Ai.
      _Tp __Ai_value;

      /// The derivative of the Airy function Ai.
      _Tp __Ai_deriv;

      /// The value of the Airy function Bi.
      _Tp __Bi_value;

      /// The derivative of the Airy function Bi.
      _Tp __Bi_deriv;

      /// Return the Wronskian of the Airy functions.
      _Tp __Wronskian() const
      { return __Ai_value * __Bi_deriv - __Bi_value * __Ai_deriv; }
    };

  /**
   * _Tp pretty much has to be complex.
   */
  template<typename _Tx, typename _Tp>
    struct __fock_airy_t
    {
      /// The argument of the Fock-type Airy fuctions.
      _Tx __x_arg;

      /// The value of the Fock-type Airy function w1.
      _Tp __w1_value;

      /// The derivative of the Fock-type Airy function w1.
      _Tp __w1_deriv;

      /// The value of the Fock-type Airy function w2.
      _Tp __w2_value;

      /// The derivative of the Fock-type Airy function w2.
      _Tp __w2_deriv;

      /// Return the Wronskian of the Fock-type Airy functions.
      _Tp __Wronskian() const
      { return __w1_value * __w2_deriv - __w2_value * __w1_deriv; }
    };

  /**
   * This struct captures the state of the cylindrical Bessel functions
   * at a given order and argument.
   */
  template<typename _Tnu, typename _Tx, typename _Tp>
    struct __cyl_bessel_t
    {
      /// The real order of the cylindrical Bessel functions.
      _Tnu __nu_arg;

      /// The argument of the cylindrical Bessel functions.
      _Tx __x_arg;

      /// The value of the Bessel function of the first kind.
      _Tp __J_value;

      /// The derivative of the Bessel function of the first kind.
      _Tp __J_deriv;

      /// The value of the Bessel function of the second kind.
      _Tp __N_value;

      /// The derivative of the Bessel function of the second kind.
      _Tp __N_deriv;

      /// Return the Wronskian of the cylindrical Bessel functions.
      _Tp __Wronskian() const
      { return __J_value * __N_deriv - __N_value * __J_deriv; }
    };

  /**
   * This struct captures the state of the Coulomb functions
   * at a given order and argument.
   */
  template<typename _Teta, typename _Trho, typename _Tp>
    struct __cyl_coulomb_t
    {
      /// The nonnegative order of the Coulomb functions.
      unsigned int __l;

      /// The real parameter of the Coulomb functions.
      _Teta __eta_arg;

      /// The argument of the Coulomb functions.
      _Trho __rho_arg;

      /// The value of the regular Coulomb function.
      _Tp __F_value;

      /// The derivative of the regular Coulomb function.
      _Tp __F_deriv;

      /// The value of the irregular Coulomb function.
      _Tp __G_value;

      /// The derivative of the irregular Coulomb function.
      _Tp __G_deriv;

      /// Return the Wronskian of the Coulomb functions.
      _Tp __Wronskian() const
      { return __F_value * __G_deriv - __G_value * __F_deriv; }
    };

  /**
   * This struct captures the state of the modified cylindrical Bessel functions
   * at a given order and argument.
   */
  template<typename _Tnu, typename _Tx, typename _Tp>
    struct __cyl_mod_bessel_t
    {
      /// The real order of the modified cylindrical Bessel functions.
      _Tnu __nu_arg;

      /// The argument of the modified cylindrical Bessel functions.
      _Tx __x_arg;

      /// The value of the modified cylindrical Bessel function
      /// of the first kind.
      _Tp __I_value;

      /// The derivative of the modified cylindrical Bessel function
      /// of the first kind.
      _Tp __I_deriv;

      /// The value of the modified cylindrical Bessel function
      /// of the second kind.
      _Tp __K_value;

      /// The derivative of the modified cylindrical Bessel function
      /// of the second kind.
      _Tp __K_deriv;

      /// Return the Wronskian of the modified cylindrical Bessel functions.
      _Tp __Wronskian() const
      { return __I_value * __K_deriv - __K_value * __I_deriv; }
    };

  /**
   * _Tp pretty much has to be complex.
   */
  template<typename _Tnu, typename _Tx, typename _Tp>
    struct __cyl_hankel_t
    {
      /// The real order of the cylindrical Hankel functions.
      _Tnu __nu_arg;

      /// The argument of the modified Hankel functions.
      _Tx __x_arg;

      /// The value of the cylindrical Hankel function of the first kind.
      _Tp __H1_value;

      /// The derivative of the cylindrical Hankel function of the first kind.
      _Tp __H1_deriv;

      /// The value of the cylindrical Hankel function of the second kind.
      _Tp __H2_value;

      /// The derivative of the cylindrical Hankel function of the second kind.
      _Tp __H2_deriv;

      /// Return the Wronskian of the cylindrical Hankel functions.
      _Tp __Wronskian() const
      { return __H1_value * __H2_deriv - __H2_value * __H1_deriv; }
    };

  template<typename _Tn, typename _Tx, typename _Tp>
    struct __sph_bessel_t
    {
      /// The integral order of the spherical Bessel functions.
      _Tn __n_arg;

      /// The argument of the spherical Bessel functions.
      _Tx __x_arg;

      /// The value of the spherical Bessel function of the first kind.
      _Tp __j_value;

      /// The derivative of the spherical Bessel function of the first kind.
      _Tp __j_deriv;

      /// The value of the spherical Bessel function of the second kind.
      _Tp __n_value;

      /// The derivative of the spherical Bessel function of the second kind.
      _Tp __n_deriv;

      /// Return the Wronskian of the spherical Bessel functions.
      _Tp __Wronskian() const
      { return __j_value * __n_deriv - __n_value * __j_deriv; }
    };

  template<typename _Tn, typename _Tx, typename _Tp>
    struct __sph_mod_bessel_t
    {
      /// The argument of the modified spherical Bessel functions.
      _Tx __x_arg;

      /// The integral order of the modified spherical Bessel functions.
      _Tn n_arg;

      /// The value of the modified spherical Bessel function
      /// of the first kind.
      _Tp __i_value;

      /// The derivative of the modified spherical Bessel function
      /// of the first kind.
      _Tp __i_deriv;

      /// The value of the modified spherical Bessel function
      /// of the second kind.
      _Tp __k_value;

      /// The derivative of the modified spherical Bessel function
      /// of the second kind.
      _Tp __k_deriv;

      /// Return the Wronskian of the modified cylindrical Bessel functions.
      _Tp __Wronskian() const
      { return __i_value * __k_deriv - __k_value * __i_deriv; }
    };

  /**
   * _Tp pretty much has to be complex.
   */
  template<typename _Tn, typename _Tx, typename _Tp>
    struct __sph_hankel_t
    {
      /// The integral order of the spherical Hankel functions.
      _Tn __n_arg;

      /// The argument of the spherical Hankel functions.
      _Tx __x_arg;

      /// The velue of the spherical Hankel function of the first kind.
      _Tp __h1_value;

      /// The derivative of the spherical Hankel function of the first kind.
      _Tp __h1_deriv;

      /// The velue of the spherical Hankel function of the second kind.
      _Tp __h2_value;

      /// The derivative of the spherical Hankel function of the second kind.
      _Tp __h2_deriv;

      /// Return the Wronskian of the cylindrical Hankel functions.
      _Tp __Wronskian() const
      { return __h1_value * __h2_deriv - __h2_value * __h1_deriv; }
    };

  template<typename _Tp>
    struct __pqgamma_t
    {
      /// 
      _Tp __pgamma_value;

      /// 
      _Tp __qgamma_value;
    };

  /**
   * The log of the absolute value of the gamma function
   * The sign of the exponentiated log(gamma) is stored in sign.
   */
  template<typename _Tp>
    struct __lgamma_t
    {
      /// The value log gamma function.
      _Tp __lgamma_value;

      /// The sign of the exponent of the log gamma value.
      int __lgamma_sign;
    };

  /**
   * The sign of the exponentiated log(gamma) is appied to the tgamma value.
   */
  template<typename _Tp>
    struct __gamma_inc_t
    {
      /// The value of the total gamma function.
      _Tp __tgamma_value;
      /// The value of the log of the incomplete gamma function
      _Tp __lgamma_value;
    };

  /**
   * @brief A structure for the gamma functions required by the Temme series
   * 	    expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   * @f[
   *   \Gamma_1 = \frac{1}{2\mu}
   * 	 \left[\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * and
   * @f[
   *   \Gamma_2 = \frac{1}{2}
   *     \left[\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   * is the nearest integer to @f$ \nu @f$.
   * The values of @f$ \Gamma(1 + \mu) @f$ and @f$ \Gamma(1 - \mu) @f$
   * are returned as well.
   *
   * The accuracy requirements on this are high for @f$ |\mu| < 0 @f$.
   */
  template<typename _Tp>
    struct __gamma_temme_t
    {
      /// The input parameter of the gamma functions
      _Tp __mu_arg;

      /// The output function @f$ 1/\Gamma(1 + \mu) @f$
      _Tp __gamma_plus_value;

      /// The output function @f$ 1/\Gamma(1 - \mu) @f$
      _Tp __gamma_minus_value;

      /// The output function @f$ \Gamma_1(\mu) @f$
      _Tp __gamma_1_value;

      /// The output function @f$ \Gamma_2(\mu) @f$
      _Tp __gamma_2_value;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SPECFUN_STATE_H
