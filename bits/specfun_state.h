// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
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

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // Slots for Jacobi elliptic function tuple.
  enum
  {
    _GLIBCXX_JACOBI_SN,
    _GLIBCXX_JACOBI_CN,
    _GLIBCXX_JACOBI_DN
  };

  template<typename _Tp>
    struct __airy_t
    {
      /// The argument of the Airy fuctions.
      _Tp x_arg;
      /// The value of the Airy function Ai.
      _Tp Ai_value;
      /// The derivative of the Airy function Ai.
      _Tp Ai_deriv;
      /// The value of the Airy function Bi.
      _Tp Bi_value;
      /// The derivative of the Airy function Bi.
      _Tp Bi_deriv;
    };

  template<typename _Tp>
    struct __cyl_bessel_t
    {
      /// The argument of the cylindrical Bessel functions.
      _Tp x_arg;
      /// The real order of the cylindrical Bessel functions.
      _Tp nu_arg;
      /// The value of the Bessel function of the first kind.
      _Tp J_value;
      /// The derivative of the Bessel function of the first kind.
      _Tp J_deriv;
      /// The value of the Bessel function of the second kind.
      _Tp N_value;
      /// The derivative of the Bessel function of the second kind.
      _Tp N_deriv;
    };

  template<typename _Tp>
    struct __cyl_mod_bessel_t
    {
      /// The argument of the modified cylindrical Bessel functions.
      _Tp x_arg;
      /// The real order of the modified cylindrical Bessel functions.
      _Tp nu_arg;
      /// The value of the modified cylindrical Bessel function
      /// of the first kind.
      _Tp I_value;
      /// The derivative of the modified cylindrical Bessel function
      /// of the first kind.
      _Tp I_deriv;
      /// The value of the modified cylindrical Bessel function
      /// of the second kind.
      _Tp K_value;
      /// The derivative of the modified cylindrical Bessel function
      /// of the second kind.
      _Tp K_deriv;
    };

  template<typename _Tp>
    struct __cyl_hankel_t
    {
      /// The argument of the modified Hankel functions.
      _Tp x_arg;
      /// The real order of the cylindrical Hankel functions.
      _Tp nu_arg;
      /// The value of the cylindrical Hankel function of the first kind.
      _Tp H1_value;
      /// The derivative of the cylindrical Hankel function of the first kind.
      _Tp H1_deriv;
      /// The value of the cylindrical Hankel function of the second kind.
      _Tp H2_value;
      /// The derivative of the cylindrical Hankel function of the second kind.
      _Tp H2_deriv;
    };

  template<typename _Tp>
    struct __sph_bessel_t
    {
      /// The argument of the spherical Bessel functions.
      _Tp x_arg;
      /// The integral order of the spherical Bessel functions.
      unsigned int n_arg;
      /// The value of the spherical Bessel function of the first kind.
      _Tp j_value;
      /// The derivative of the spherical Bessel function of the first kind.
      _Tp j_deriv;
      /// The value of the spherical Bessel function of the second kind.
      _Tp n_value;
      /// The derivative of the spherical Bessel function of the second kind.
      _Tp n_deriv;
    };

  template<typename _Tp>
    struct __sph_mod_bessel_t
    {
      /// The argument of the modified spherical Bessel functions.
      _Tp x_arg;
      /// The integral order of the modified spherical Bessel functions.
      unsigned int n_arg;
      /// The value of the modified spherical Bessel function
      /// of the first kind.
      _Tp i_value;
      /// The derivative of the modified spherical Bessel function
      /// of the first kind.
      _Tp i_deriv;
      /// The value of the modified spherical Bessel function
      /// of the second kind.
      _Tp k_value;
      /// The derivative of the modified spherical Bessel function
      /// of the second kind.
      _Tp k_deriv;
    };

  template<typename _Tp>
    struct __sph_hankel_t
    {
      /// The argument of the spherical Hankel functions.
      _Tp x_arg;
      /// The integral order of the spherical Hankel functions.
      unsigned int n_arg;
      /// The velue of the spherical Hankel function of the first kind.
      _Tp h1_value;
      /// The derivative of the spherical Hankel function of the first kind.
      _Tp h1_deriv;
      /// The velue of the spherical Hankel function of the second kind.
      _Tp h2_value;
      /// The derivative of the spherical Hankel function of the second kind.
      _Tp h2_deriv;
    };

  template<typename _Tp>
    struct __pqgamma_t
    {
      /// 
      _Tp pgamma_value;
      /// 
      _Tp qgamma_value;
    };

  /**
   * The log of the absolute value of the gamma function
   * The sign of the exponentiated log(gamma) is stored in sign.
   */
  template<typename _Tp>
    struct __lgamma_t
    {
      /// The value log gamma function.
      _Tp lgamma_value;
      /// The sign of the exponent of the log gamma value.
      int lgamma_sign;
    };

  /**
   * The sign of the exponentiated log(gamma) is appied to the tgamma value.
   */
  template<typename _Tp>
    struct __gamma_inc_t
    {
      /// The value of the total gamma function.
      _Tp tgamma_value;
      /// The value of the log of the incomplete gamma function
      _Tp lgamma_value;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SPECFUN_STATE_H
