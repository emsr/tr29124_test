// Special functions -*- C++ -*-

// Copyright (C) 2006-2018 Free Software Foundation, Inc.
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
   * \brief Enumeration gor differing types of Gauss quadrature.
   * The gauss_quad_type is used to determine the boundary condition
   * modifications applied to orthogonal polynomials for quadrature
   * rules.
   */
  enum gauss_quad_type
  {
    Gauss,             ///< Gauss quadrature
    Gauss_Lobatto,     ///< Gauss-Lobatto quadrature
    Gauss_Radau_lower, ///< Gauss-Radau quadrature including the node -1
    Gauss_Radau_upper  ///< Gauss-Radau quadrature including the node +1
  };

  /**
   * A structure to store quadrature rules.
   */
  template<typename _Tp>
    struct __quadrature_point_t
    {
      _Tp __point;
      _Tp __weight;

      __quadrature_point_t() = default;

      __quadrature_point_t(_Tp __pt, _Tp __wt)
      : __point(__pt),
	__weight(__wt)
      { }
    };

  /**
   * A struct to store the state of a Hermite polynomial.
   */
  template<typename _Tp>
    struct __hermite_t
    {
      std::size_t __n;
      _Tp __x;
      _Tp __H_n;
      _Tp __H_nm1;
      _Tp __H_nm2;

      _Tp
      deriv() const
      { return _Tp(2 * __n) * __H_nm1; }

      _Tp
      deriv2() const
      { return _Tp(4 * __n * (__n - 1)) * __H_nm2; }
    };

  /**
   * A struct to store the state of a probabilists Hermite polynomial.
   */
  template<typename _Tp>
    struct __hermite_he_t
    {
      std::size_t __n;
      _Tp __x;
      _Tp __He_n;
      _Tp __He_nm1;
      _Tp __He_nm2;

      _Tp
      deriv() const
      { return _Tp(__n) * __He_nm1; }

      _Tp
      deriv2() const
      { return _Tp(__n * (__n - 1)) * __He_nm2; }
    };

  /**
   * A struct to store the state of a Legendre polynomial.
   */
  template<typename _Tp>
    struct __legendre_p_t
    {
      std::size_t __l;
      _Tp __z;
      _Tp __P_l;
      _Tp __P_lm1;
      _Tp __P_lm2;

      // @todo endpoints?
      _Tp
      deriv() const
      { return __l * (__z * __P_l - __P_lm1) / (__z * __z - _Tp{1}); }
    };

  /**
   * A struct to store the state of a Laguerre polynomial.
   */
  template<typename _Tpa, typename _Tp>
    struct __laguerre_t
    {
      std::size_t __n;
      _Tpa __alpha1;
      _Tp __x;
      _Tp __L_n;
      _Tp __L_nm1;
      _Tp __L_nm2;

      _Tp
      deriv() const
      { return (_Tp(__n) * __L_nm1 - _Tp(__n + __alpha1) * __L_nm2) / __x; }
    };

  /**
   * A struct to store the state of a Jacobi polynomial.
   */
  template<typename _Tp>
    struct __jacobi_t
    {
      std::size_t __n;
      _Tp __alpha1;
      _Tp __beta1;
      _Tp __x;
      _Tp __P_n;
      _Tp __P_nm1;
      _Tp __P_nm2;

      _Tp
      deriv() const
      {
	auto __apbp2k = __alpha1 + __beta1 + _Tp(2 * __n);
	return (__n * (__alpha1 - __beta1 - __apbp2k * __x) * __P_nm1
		   + _Tp{2} * (__n + __alpha1) * (__n + __beta1) * __P_nm2)
		/ (__apbp2k * (_Tp{1} - __x * __x));
      }
    };

  /**
   * A struct to store the state of a Gegenbauer polynomial.
   */
  template<typename _Tp>
    struct __gegenbauer_t
    {
      std::size_t __n;
      _Tp __alpha1;
      _Tp __x;
      _Tp __C_n;
      _Tp __C_nm1;
      _Tp __C_nm2;

      _Tp
      deriv() const
      {
	auto __apbp2k = _Tp{2} * __alpha1 + _Tp(2 * __n);
	return (__n * (-__apbp2k * __x) * __C_nm1
		   + _Tp{2} * (__n + __alpha1) * (__n + __alpha1) * __C_nm2)
		/ (__apbp2k * (_Tp{1} - __x * __x));
      }
    };

  /**
   * A struct to store the state of a Chebyshev polynomial of the first kind.
   */
  template<typename _Tp>
    struct __chebyshev_t_t
    {
      std::size_t __n;
      _Tp __x;
      _Tp __T_n;
      _Tp __T_nm1;
      _Tp __T_nm2;

      _Tp
      deriv() const
      { return _Tp(__n) * (__T_nm1 - __x * __T_n) / (_Tp{1} - __x * __x); }

      _Tp
      deriv2() const
      {
	const auto __xx = __x * __x;
	const auto __num = _Tp{1} - __xx;
	return _Tp(__n)
	     * (__x * __T_nm1 + (_Tp(__n - 1) * __xx - _Tp(__n)) * __T_n)
	     / __num / __num;
      }
    };

  /**
   * A struct to store the state of a Chebyshev polynomial of the second kind.
   */
  template<typename _Tp>
    struct __chebyshev_u_t
    {
      std::size_t __n;
      _Tp __x;
      _Tp __U_n;
      _Tp __U_nm1;
      _Tp __U_nm2;

      _Tp
      deriv() const
      {
	return (_Tp(__n + 1) * __U_nm1 - _Tp(__n) * __x * __U_n)
		/ (_Tp{1} - __x * __x);
      }
    };

  /**
   * A struct to store the state of a Chebyshev polynomial of the third kind.
   */
  template<typename _Tp>
    struct __chebyshev_v_t
    {
      std::size_t __n;
      _Tp __x;
      _Tp __V_n;
      _Tp __V_nm1;
      _Tp __V_nm2;

      _Tp
      deriv() const
      {
	auto __apbp2k = _Tp(2 * __n);
	return (__n * (_Tp{1} - __apbp2k * __x) * __V_nm1
		   + _Tp(2 * (__n + 0.5L) * (__n + -0.5L)) * __V_nm2)
		/ (__apbp2k * (_Tp{1} - __x * __x));
      }
    };

  /**
   * A struct to store the state of a Chebyshev polynomial of the fourth kind.
   */
  template<typename _Tp>
    struct __chebyshev_w_t
    {
      std::size_t __n;
      _Tp __x;
      _Tp __W_n;
      _Tp __W_nm1;
      _Tp __W_nm2;

      _Tp
      deriv() const
      {
	auto __apbp2k = _Tp(2 * __n);
	return (__n * (_Tp{-1} - __apbp2k * __x) * __W_nm1
		   + _Tp(2 * (__n - 0.5L) * (__n + 0.5L)) * __W_nm2)
		/ (__apbp2k * (_Tp{1} - __x * __x));
      }
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

  /**
   * Slots for Jacobi elliptic function tuple.
   */
  template<typename _Tp>
    struct __jacobi_ellint_t
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

      _Tp __sn_deriv() const
      { return __cn_value * __dn_value; }

      _Tp __cn_deriv() const
      { return -__sn_value * __dn_value; }
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

      /// Return the Wronskian of this Airy function state.
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

      /// Return the Wronskian of this Fock-type Airy function state.
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

      /// Return the Wronskian of this cylindrical Bessel function state.
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

      /// Return the Wronskian of this Coulomb function state.
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

      /// Return the Wronskian of this modified cylindrical Bessel function
      /// state.
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

      /// Return the Wronskian of this cylindrical Hankel function state.
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

      /// Return the Wronskian of this spherical Bessel function state.
      _Tp __Wronskian() const
      { return __j_value * __n_deriv - __n_value * __j_deriv; }
    };

  template<typename _Tn, typename _Tx, typename _Tp>
    struct __sph_mod_bessel_t
    {
      /// The argument of the modified spherical Bessel functions.
      _Tx __x_arg;

      /// The integral order of the modified spherical Bessel functions.
      _Tn __n_arg;

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

      /// Return the Wronskian of this modified cylindrical Bessel function
      /// state.
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

      /// Return the Wronskian of this cylindrical Hankel function state.
      _Tp __Wronskian() const
      { return __h1_value * __h2_deriv - __h2_value * __h1_deriv; }
    };

  template<typename _Tp>
    struct __gappa_pq_t
    {
      /// 
      _Tp __gappa_p_value;

      /// 
      _Tp __gappa_q_value;
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