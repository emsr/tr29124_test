
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

/** @file new_bessel_function.tcc
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

#ifndef NEW_BESSEL_FUNCTION_TCC
#define NEW_BESSEL_FUNCTION_TCC 1

#include <emsr/fp_type_util.h>

namespace std
{
namespace detail
{

    ///
    ///
    ///
    template <typename Tp>
    void
    bessel_jn(const Tp nu, const Tp x,
                Tp & J_nu, Tp & N_nu, Tp & Jp_nu, Tp & Np_nu)
    {

      if (std::isnan(nu) || std::isnan(x))
        {
          J_nu = std::numeric_limits<Tp>::quiet_NaN();
          N_nu = std::numeric_limits<Tp>::quiet_NaN();
          Jp_nu = std::numeric_limits<Tp>::quiet_NaN();
          Np_nu = std::numeric_limits<Tp>::quiet_NaN();
          return;
        }

      if (x == Tp(0))
        {
          if (nu == Tp(0))
            {
              J_nu = Tp(1);
              N_nu = -std::numeric_limits<Tp>::infinity();
              Jp_nu = Tp(0);
              Np_nu = std::numeric_limits<Tp>::infinity();
            }
          else
            {
              J_nu = Tp(0);
              N_nu = -std::numeric_limits<Tp>::infinity();
              //Jp_nu = ???
              //Np_nu = ???
            }
          return;
        }

      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp fp_min = Tp(10) * std::numeric_limits<Tp>::min();
      const int max_iter = 15000;
      const Tp x_min = Tp(2);
      const Tp PI = Tp(3.1415926535897932384626433832795029L);

      if (x < Tp(0) || nu < Tp(0))
        throw_runtime_error(N("Bad arguments in bessel_jn."));

      const int nl = (x < x_min
                    ? static_cast<int>(nu + Tp(0.5L))
                    : std::max(0, static_cast<int>(nu - x + Tp(1.5L))));

      const Tp mu = nu - nl;
      const Tp mu2 = mu * mu;
      const Tp xi = Tp(1) / x;
      const Tp xi2 = Tp(2) * xi;
      Tp w = xi2 / PI;
      int isign = 1;
      Tp h = nu * xi;
      if (h < fp_min)
        h = fp_min;
      Tp b = xi2 * nu;
      Tp d = Tp(0);
      Tp c = h;
      int i;
      for (i = 1; i <= max_iter; ++i)
        {
          b += xi2;
          d = b - d;
          if (std::abs(d) < fp_min)
            d = fp_min;
          c = b - Tp(1) / c;
          if (std::abs(c) < fp_min)
            c = fp_min;
          d = Tp(1) / d;
          const Tp del = c * d;
          h = del * h;
          if (d < Tp(0))
            isign = -isign;
          if (std::abs(del - Tp(1)) < eps)
            break;
        }
      if (i > max_iter)
        throw_runtime_error(N("Argument x too large in bessel_jn; "
                                  "try asymptotic expansion."));
      Tp J_nul = isign * fp_min;
      Tp Jp_nul = h * J_nul;
      Tp J_nul1 = J_nul;
      Tp Jp_nu1 = Jp_nul;
      Tp fact = nu * xi;
      for ( int l = nl; l >= 1; --l )
        {
          const Tp J_nutemp = fact * J_nul + Jp_nul;
          fact -= xi;
          Jp_nul = fact * J_nutemp - J_nul;
          J_nul = J_nutemp;
        }
      if (J_nul == Tp(0))
        J_nul = eps;
      Tp f= Jp_nul / J_nul;
      Tp N_mu, N_nu1, N_mup, J_mu;
      if (x < x_min)
        {
          const Tp x2 = x / Tp(2);
          const Tp pimu = PI * mu;
          Tp fact = (std::abs(pimu) < eps ? Tp(1)
                      : pimu / std::sin(pimu));
          Tp d = -std::log(x2);
          Tp e = mu * d;
          Tp fact2 = (std::abs(e) < eps
                       ? Tp(1) : std::sinh(e) / e);
          Tp gam1, gam2, gampl, gammi;
          gamma_temme(mu, gampl, gammi, gam1, gam2);
          Tp ff = (Tp(2) / PI) * fact
                   * (gam1 * std::cosh(e) + gam2 * fact2 * d);
          e = std::exp(e);
          Tp p = e / (PI * gampl);
          Tp q = Tp(1) / (e * PI * gammi);
          const Tp pimu2 = pimu / Tp(2);
          Tp fact3 = (std::abs(pimu2) < eps
                       ? Tp(1) : std::sin(pimu2) / pimu2 );
          Tp r = PI * pimu2 * fact3 * fact3;
          Tp c = Tp(1);
          d = -x2 * x2;
          Tp sum = ff + r * q;
          Tp sum1 = p;
          for (i = 1; i <= max_iter; ++i)
            {
              ff = (i * ff + p + q) / (i * i - mu2);
              c *= d / Tp(i);
              p /= Tp(i) - mu;
              q /= Tp(i) + mu;
              const Tp del = c * (ff + r * q);
              sum += del; 
              const Tp del1 = c * p - i * del;
              sum1 += del1;
              if ( std::abs(del) < eps * (Tp(1) + std::abs(sum)) )
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
          Tp a = Tp(0.25L) - mu2;
          Tp q = Tp(1);
          Tp p = -xi / Tp(2);
          Tp br = Tp(2) * x;
          Tp bi = Tp(2);
          Tp fact = a * xi / (p * p + q * q);
          Tp cr = br + q * fact;
          Tp ci = bi + p * fact;
          Tp den = br * br + bi * bi;
          Tp dr = br / den;
          Tp di = -bi / den;
          Tp dlr = cr * dr - ci * di;
          Tp dli = cr * di + ci * dr;
          Tp temp = p * dlr - q * dli;
          q = p * dli + q * dlr;
          p = temp;
          int i;
          for (i = 2; i <= max_iter; ++i)
            {
              a += Tp(2 * (i - 1));
              bi += Tp(2);
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
              if (std::abs(dlr - Tp(1)) + std::abs(dli) < eps)
                break;
          }
          if (i > max_iter)
            throw_runtime_error(N("Lentz's method failed "
                                      "in bessel_jn."));
          const Tp gam = (p - f) / q;
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
          const Tp N_nutemp = (mu + i) * xi2 * N_nu1 - N_mu;
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
    template <typename Tp>
    void
    sph_bessel_jn(const int n, const Tp x,
                    Tp & j_n, Tp & n_n, Tp & jp_n, Tp & np_n)
    {

      if ( n < 0 || x < Tp(0) )
        throw_runtime_error(N("Bad arguments in sph_bessel."));

      const Tp nu = Tp(n) + Tp(0.5L);

      Tp J_nu, Jp_nu, N_nu, Np_nu;
      bessel_jn(nu, x, J_nu, N_nu, Jp_nu, Np_nu);

      const Tp SQRTPIO2 = Tp(1.2533141373155002512078826424055226L);
      const Tp factor = SQRTPIO2 / std::sqrt(x);

      j_n = factor * J_nu;
      n_n = factor * N_nu;
      jp_n = factor * Jp_nu - j_n / (Tp(2) * x);
      np_n = factor * Np_nu - n_n / (Tp(2) * x);

      return;
    }

} // detail
} // emsr

#endif // NEW_BESSEL_FUNCTION_TCC

