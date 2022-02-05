// Special functions -*- C++ -*-

// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file good_modified_bessel_func.tcc
 */

//
// ISO C++ 14882 TR1: 5.2  Special functions
//

// Written by Edward Smith-Rowland based on numerous mathematics books.

#ifndef GOOD_MODIFIED_BESSEL_FUNC_TCC
#define GOOD_MODIFIED_BESSEL_FUNC_TCC 1

#include <emsr/fp_type_util.h>

namespace std
{
namespace detail
{

    ///
    ///  @brief  Compute the modified Bessel functions @f$ I_\nu(x) @f$ and
    ///          @f$ K_\nu(x) @f$ and their first derivatives
    ///          @f$ I'_\nu(x) @f$ and @f$ K'_\nu(x) @f$ respectively.
    ///          These four functions are computed together for numerical
    ///          stability.
    ///
    template <typename Tp>
    void
    bessel_ik(const Tp nu, const Tp x,
                Tp & I_nu, Tp & K_nu, Tp & Ip_nu, Tp & Kp_nu)
    {

      if (x < Tp(0) || nu < Tp(0))
        throw_runtime_error(N("Bad arguments in "
                                  "bessel_ik."));

      if (std::isnan(nu) || std::isnan(x))
        {
          I_nu = std::numeric_limits<Tp>::quiet_NaN();
          K_nu = std::numeric_limits<Tp>::quiet_NaN();
          Ip_nu = std::numeric_limits<Tp>::quiet_NaN();
          Kp_nu = std::numeric_limits<Tp>::quiet_NaN();
          return;
        }

      if (x == Tp(0))
        {
          if (nu == Tp(0))
            {
              I_nu = Tp(1);
              K_nu = std::numeric_limits<Tp>::infinity();
              Ip_nu = Tp(0);
              Kp_nu = -std::numeric_limits<Tp>::infinity();
            }
          else
            {
              I_nu = Tp(0);
              K_nu = std::numeric_limits<Tp>::infinity();
              //Ip_nu = ???
              //Kp_nu = ???
            }
          return;
        }

      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp fp_min = Tp(10) * std::numeric_limits<Tp>::epsilon();
      const int max_iter = 15000;
      const Tp x_min = Tp(2);

      const int nl = static_cast<int>(nu + Tp(0.5L));

      const Tp mu = nu - nl;
      const Tp mu2 = mu * mu;
      const Tp xi = Tp(1) / x;
      const Tp xi2 = Tp(2) * xi;
      Tp h = nu * xi;
      if ( h < fp_min )
        h = fp_min;
      Tp b = xi2 * nu;
      Tp d = Tp(0);
      Tp c = h;
      int i;
      for ( i = 1; i <= max_iter; ++i )
        {
          b += xi2;
          d = Tp(1) / (b + d);
          c = b + Tp(1) / c;
          const Tp del = c * d;
          h = del * h;
          if (std::abs(del - Tp(1)) < eps)
            break;
        }
      if (i > max_iter)
        throw_runtime_error(N("Argument x too large in bessel_ik; "
                                  "try asymptotic expansion."));
      Tp I_nul = fp_min;
      Tp Ip_nul = h * I_nul;
      Tp I_nul1 = I_nul;
      Tp Ip_nu1 = Ip_nul;
      Tp fact = nu * xi;
      for (int l = nl; l >= 1; --l)
        {
          const Tp I_nutemp = fact * I_nul + Ip_nul;
          fact -= xi;
          Ip_nul = fact * I_nutemp + I_nul;
          I_nul = I_nutemp;
        }
      Tp f = Ip_nul / I_nul;
      Tp K_mu, K_nu1;
      if (x < x_min)
        {
          const Tp x2 = x / Tp(2);
          const Tp pimu = PI * mu;
          const Tp fact = (std::abs(pimu) < eps
                            ? Tp(1) : pimu / std::sin(pimu));
          Tp d = -std::log(x2);
          Tp e = mu * d;
          const Tp fact2 = (std::abs(e) < eps
                            ? Tp(1) : sinh(e) / e);
          Tp gam1, gam2, gampl, gammi;
          gamma_temme(mu, gampl, gammi, gam1, gam2);
          Tp ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
          Tp sum = ff;
          e = std::exp(e);
          Tp p = e / (Tp(2) * gampl);
          Tp q = Tp(1) / (Tp(2) * e * gammi);
          Tp c = Tp(1);
          d = x2 * x2;
          Tp sum1 = p;
          int i;
          for (i = 1; i <= max_iter; ++i)
            {
              ff = (i * ff + p + q) / (i * i - mu2);
              c *= d / i;
              p /= i - mu;
              q /= i + mu;
              const Tp del = c * ff;
              sum += del; 
              const Tp del1 = c * (p - i * ff);
              sum1 += del1;
              if (std::abs(del) < eps * std::abs(sum))
                break;
            }
          if (i > max_iter)
            throw_runtime_error(N("Bessel k series failed to converge "
                                      "in bessel_ik."));
          K_mu = sum;
          K_nu1 = sum1 * xi2;
        }
      else
        {
          Tp b = Tp(2) * (Tp(1) + x);
          Tp d = Tp(1) / b;
          Tp delh = d;
          Tp h = delh;
          Tp q1 = Tp(0);
          Tp q2 = Tp(1);
          Tp a1 = 0.25L - mu2;
          Tp q = c = a1;
          Tp a = -a1;
          Tp s = Tp(1) + q * delh;
          int i;
          for (i = 2; i <= max_iter; ++i)
            {
              a -= 2 * (i - 1);
              c = -a * c / i;
              const Tp qnew = (q1 - b * q2) / a;
              q1 = q2;
              q2 = qnew;
              q += c * qnew;
              b += Tp(2);
              d = Tp(1) / (b + a * d);
              delh = (b * d - Tp(1)) * delh;
              h += delh;
              const Tp dels = q * delh;
              s += dels;
              if ( std::abs(dels / s) < eps )
                break;
            }
          if (i > max_iter)
            throw_runtime_error(N("Steed's method failed "
                                      "in bessel_ik."));
          h = a1 * h;
          K_mu = std::sqrt(PI / (Tp(2) * x)) * std::exp(-x) / s;
          K_nu1 = K_mu * (mu + x + Tp(0.5L) - h) * xi;
        }

      Tp K_mup = mu * xi * K_mu - K_nu1;
      Tp I_mu = xi / (f * K_mu - K_mup);
      I_nu = I_mu * I_nul1 / I_nul;
      Ip_nu = I_mu * Ip_nu1 / I_nul;
      for ( i = 1; i <= nl; ++i )
        {
          const Tp K_nutemp = (mu + i) * xi2 * K_nu1 + K_mu;
          K_mu = K_nu1;
          K_nu1 = K_nutemp;
        }
      K_nu = K_mu;
      Kp_nu = nu * xi * K_mu - K_nu1;

      return;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ I_{\nu}(x) \f$.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i(const Tp nu, const Tp x)
    {
      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_i.");

      if (std::isnan(nu) || std::isnan(x))
        return std::numeric_limits<Tp>::quiet_NaN();

      Tp I_nu, K_nu, Ip_nu, Kp_nu;
      bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);

      return I_nu;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu}(x) \f$.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_k(const Tp nu, const Tp x)
    {
      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_k.");

      if (std::isnan(nu) || std::isnan(x))
        return std::numeric_limits<Tp>::quiet_NaN();

      Tp I_nu, K_nu, Ip_nu, Kp_nu;
      bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);

      return K_nu;
    }


    ///
    ///  @brief  Compute the spherical modified Bessel functions
    ///          @f$ i_n(x) @f$ and @f$ k_n(x) @f$ and their first
    ///          derivatives @f$ i'_n(x) @f$ and @f$ k'_n(x) @f$
    ///          respectively.
    ///
    template <typename Tp>
    void
    sph_bessel_ik(const unsigned int n, const Tp x,
                    Tp & i_n, Tp & k_n, Tp & ip_n, Tp & kp_n)
    {

      if (x < Tp(0))
        throw_runtime_error(N("Bad arguments "
                                  "in sph_bessel."));

      const Tp nu = Tp(n) + Tp(0.5L);

      Tp I_nu, Ip_nu, K_nu, Kp_nu;
      bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);

      const Tp factor = Tp(SQRTPIO2) / std::sqrt(x);

      i_n = factor * I_nu;
      k_n = factor * K_nu;
      ip_n = factor * Ip_nu - i_n / (Tp(2) * x);
      kp_n = factor * Kp_nu - k_n / (Tp(2) * x);

      return;
    }


    ///
    ///  @brief  Compute the Airy functions
    ///          @f$ Ai(x) @f$ and @f$ Bi(x) @f$ and their first
    ///          derivatives @f$ Ai'(x) @f$ and @f$ Bi(x) @f$
    ///          respectively.
    ///
    template <typename Tp>
    void
    airy(const Tp x,
           Tp & Ai, Tp & Bi, Tp & Aip, Tp & Bip)
    {
      const Tp absx = std::abs(x);
      const Tp rootx = std::sqrt(absx);
      const Tp z = Tp(2) * absx * rootx / Tp(3);
      if (x > Tp(0))
        {
          Tp I_nu, Ip_nu, K_nu, Kp_nu;

          bessel_jn(Tp(1)/Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu);
          Ai = rootx * K_nu / (SQRT3 * PI);
          Bi = rootx * (K_nu / PI + Tp(2) * I_nu / SQRT3);

          bessel_jn(Tp(2)/Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu);
          Aip = -x * K_nu / (SQRT3 * PI);
          Bip = x * (K_nu / PI + Tp(2) * I_nu / SQRT3);
        }
      else if (x < Tp(0))
        {
          Tp J_nu, Jp_nu, N_nu, Np_nu;

          bessel_jn(Tp(1)/Tp(3), z, J_nu, N_nu, Jp_nu, Np_nu);
          Ai = rootx * (J_nu - N_nu / SQRT3) / Tp(2);
          Bi = -rootx * (N_nu + J_nu / SQRT3) / Tp(2);

          bessel_jn(Tp(2)/Tp(3), z, J_nu, N_nu, Jp_nu, Np_nu);
          Aip = absx * (N_nu / SQRT3 + J_nu) / Tp(2);
          Bip = absx * (J_nu / SQRT3 - N_nu) / Tp(2);
        }
      else
        {
          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
          // The number is Ai(0) = 3^{-2/3}/\Gamma(2/3).
          Ai = Tp(0.35502805388781723926L);
          Bi = Ai * SQRT3;

          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
          // The number is Ai'(0) = -3^{-1/3}/\Gamma(1/3).
          Aip = -Tp(0.25881940379280679840L);
          Bip = -Aip * SQRT3;
        }

      return;
    }

} // namespace detail
} // namespace emsr

#endif // GOOD_MODIFIED_BESSEL_FUNC_TCC
