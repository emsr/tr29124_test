// Special functions -*- C++ -*-

// Copyright (C) 2006-2019 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.

/** @file tr1/sf_bessel.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

//
// ISO C++ 14882 TR1: 5.2  Special functions
//

// Written by Edward Smith-Rowland.
//
// References:
//   (1) Handbook of Mathematical Functions,
//       ed. Milton Abramowitz and Irene A. Stegun,
//       Dover Publications,
//       Section 9 pp. 355-434, Section 10 pp. 435-478
//   (2) Numerical Recipes in C,
//       

#ifndef _TR1_BESSEL_FUNCTION_TCC
#define _TR1_BESSEL_FUNCTION_TCC 1

#include <emsr/fp_type_util.h>

namespace std
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // [5.2] Special functions

  /**
   * @addtogroup tr1_math_spec_func Mathematical Special Functions
   * A collection of advanced mathematical special functions.
   * @{
   */

  //
  // Implementation-space details.
  //
  namespace __detail
  {

    ///
    ///
    ///
    template <typename _Tp>
    void
    bessel_jn(const _Tp nu, const _Tp x,
                _Tp & J_nu, _Tp & N_nu, _Tp & Jp_nu, _Tp & Np_nu)
    {

      if (std::isnan(nu) || std::isnan(x))
        {
          J_nu = std::numeric_limits<_Tp>::quiet_NaN();
          N_nu = std::numeric_limits<_Tp>::quiet_NaN();
          Jp_nu = std::numeric_limits<_Tp>::quiet_NaN();
          Np_nu = std::numeric_limits<_Tp>::quiet_NaN();
          return;
        }

      if (x == _Tp(0))
        {
          if (nu == _Tp(0))
            {
              J_nu = _Tp(1);
              N_nu = -std::numeric_limits<_Tp>::infinity();
              Jp_nu = _Tp(0);
              Np_nu = std::numeric_limits<_Tp>::infinity();
            }
          else
            {
              J_nu = _Tp(0);
              N_nu = -std::numeric_limits<_Tp>::infinity();
              //Jp_nu = ???
              //Np_nu = ???
            }
          return;
        }

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp fp_min = _Tp(10) * std::numeric_limits<_Tp>::min();
      const int max_iter = 15000;
      const _Tp x_min = _Tp(2);
      const _Tp PI = _Tp(3.1415926535897932384626433832795029L);

      if (x < _Tp(0) || nu < _Tp(0))
        throw_runtime_error(N("Bad arguments in bessel_jn."));

      const int nl = (x < x_min
                    ? static_cast<int>(nu + _Tp(0.5L))
                    : std::max(0, static_cast<int>(nu - x + _Tp(1.5L))));

      const _Tp mu = nu - nl;
      const _Tp mu2 = mu * mu;
      const _Tp xi = _Tp(1) / x;
      const _Tp xi2 = _Tp(2) * xi;
      _Tp w = xi2 / PI;
      int isign = 1;
      _Tp h = nu * xi;
      if (h < fp_min)
        h = fp_min;
      _Tp b = xi2 * nu;
      _Tp d = _Tp(0);
      _Tp c = h;
      int i;
      for (i = 1; i <= max_iter; ++i)
        {
          b += xi2;
          d = b - d;
          if (std::abs(d) < fp_min)
            d = fp_min;
          c = b - _Tp(1) / c;
          if (std::abs(c) < fp_min)
            c = fp_min;
          d = _Tp(1) / d;
          const _Tp del = c * d;
          h = del * h;
          if (d < _Tp(0))
            isign = -isign;
          if (std::abs(del - _Tp(1)) < eps)
            break;
        }
      if (i > max_iter)
        throw_runtime_error(N("Argument x too large in bessel_jn; "
                                  "try asymptotic expansion."));
      _Tp J_nul = isign * fp_min;
      _Tp Jp_nul = h * J_nul;
      _Tp J_nul1 = J_nul;
      _Tp Jp_nu1 = Jp_nul;
      _Tp fact = nu * xi;
      for ( int l = nl; l >= 1; --l )
        {
          const _Tp J_nutemp = fact * J_nul + Jp_nul;
          fact -= xi;
          Jp_nul = fact * J_nutemp - J_nul;
          J_nul = J_nutemp;
        }
      if (J_nul == _Tp(0))
        J_nul = eps;
      _Tp f= Jp_nul / J_nul;
      _Tp N_mu, N_nu1, N_mup, J_mu;
      if (x < x_min)
        {
          const _Tp x2 = x / _Tp(2);
          const _Tp pimu = PI * mu;
          _Tp fact = (std::abs(pimu) < eps ? _Tp(1)
                      : pimu / std::sin(pimu));
          _Tp d = -std::log(x2);
          _Tp e = mu * d;
          _Tp fact2 = (std::abs(e) < eps
                       ? _Tp(1) : std::sinh(e) / e);
          _Tp gam1, gam2, gampl, gammi;
          gamma_temme(mu, gampl, gammi, gam1, gam2);
          _Tp ff = (_Tp(2) / PI) * fact
                   * (gam1 * std::cosh(e) + gam2 * fact2 * d);
          e = std::exp(e);
          _Tp p = e / (PI * gampl);
          _Tp q = _Tp(1) / (e * PI * gammi);
          const _Tp pimu2 = pimu / _Tp(2);
          _Tp fact3 = (std::abs(pimu2) < eps
                       ? _Tp(1) : std::sin(pimu2) / pimu2 );
          _Tp r = PI * pimu2 * fact3 * fact3;
          _Tp c = _Tp(1);
          d = -x2 * x2;
          _Tp sum = ff + r * q;
          _Tp sum1 = p;
          for (i = 1; i <= max_iter; ++i)
            {
              ff = (i * ff + p + q) / (i * i - mu2);
              c *= d / _Tp(i);
              p /= _Tp(i) - mu;
              q /= _Tp(i) + mu;
              const _Tp del = c * (ff + r * q);
              sum += del; 
              const _Tp del1 = c * p - i * del;
              sum1 += del1;
              if ( std::abs(del) < eps * (_Tp(1) + std::abs(sum)) )
                break;
            }
          if ( i > max_iter )
            throw_runtime_error(N("Bessel y series failed to converge "
                                      "in bessel_jn."));
          N_mu = -sum;
          N_nu1 = -sum1 * xi2;
          N_mup = mu * xi * N_mu - N_nu1;
          J_mu = w / (N_mup - f * N_mu);
        }
      else
        {
          _Tp a = _Tp(0.25L) - mu2;
          _Tp q = _Tp(1);
          _Tp p = -xi / _Tp(2);
          _Tp br = _Tp(2) * x;
          _Tp bi = _Tp(2);
          _Tp fact = a * xi / (p * p + q * q);
          _Tp cr = br + q * fact;
          _Tp ci = bi + p * fact;
          _Tp den = br * br + bi * bi;
          _Tp dr = br / den;
          _Tp di = -bi / den;
          _Tp dlr = cr * dr - ci * di;
          _Tp dli = cr * di + ci * dr;
          _Tp temp = p * dlr - q * dli;
          q = p * dli + q * dlr;
          p = temp;
          int i;
          for (i = 2; i <= max_iter; ++i)
            {
              a += _Tp(2 * (i - 1));
              bi += _Tp(2);
              dr = a * dr + br;
              di = a * di + bi;
              if (std::abs(dr) + std::abs(di) < fp_min)
                dr = fp_min;
              fact = a / (cr * cr + ci * ci);
              cr = br + cr * fact;
              ci = bi - ci * fact;
              if (std::abs(cr) + std::abs(ci) < fp_min)
                cr = fp_min;
              den = dr * dr + di * di;
              dr /= den;
              di /= -den;
              dlr = cr * dr - ci * di;
              dli = cr * di + ci * dr;
              temp = p * dlr - q * dli;
              q = p * dli + q * dlr;
              p = temp;
              if (std::abs(dlr - _Tp(1)) + std::abs(dli) < eps)
                break;
          }
          if (i > max_iter)
            throw_runtime_error(N("Lentz's method failed "
                                      "in bessel_jn."));
          const _Tp gam = (p - f) / q;
          J_mu = std::sqrt(w / ((p - f) * gam + q));
          J_mu = ::copysign( J_mu, J_nul );
          N_mu = J_mu * gam;
          N_mup = N_mu * (p + q / gam);
          N_nu1 = mu * xi * N_mu - N_mup;
      }
      fact = J_mu / J_nul;
      J_nu = J_nul1 * fact;
      Jp_nu = Jp_nu1 * fact;
      for (i = 1; i <= nl; ++i)
        {
          const _Tp N_nutemp = (mu + i) * xi2 * N_nu1 - N_mu;
          N_mu = N_nu1;
          N_nu1 = N_nutemp;
        }
      N_nu = N_mu;
      Np_nu = nu * xi * N_mu - N_nu1;

      return;
    }


    ///
    ///
    ///
    template <typename _Tp>
    void
    sph_bessel_jn(const int n, const _Tp x,
                    _Tp & j_n, _Tp & n_n, _Tp & jp_n, _Tp & np_n)
    {

      if ( n < 0 || x < _Tp(0) )
        throw_runtime_error(N("Bad arguments in sph_bessel."));

      const _Tp nu = _Tp(n) + _Tp(0.5L);

      _Tp J_nu, Jp_nu, N_nu, Np_nu;
      bessel_jn(nu, x, J_nu, N_nu, Jp_nu, Np_nu);

      const _Tp SQRTPIO2 = _Tp(1.2533141373155002512078826424055226L);
      const _Tp factor = SQRTPIO2 / std::sqrt(x);

      j_n = factor * J_nu;
      n_n = factor * N_nu;
      jp_n = factor * Jp_nu - j_n / (_Tp(2) * x);
      np_n = factor * Np_nu - n_n / (_Tp(2) * x);

      return;
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _TR1_BESSEL_FUNCTION_TCC

