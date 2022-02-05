
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

/** @file old_modified_bessel_func.tcc
 */

//
// ISO C++ 14882 TR1: 5.2  Special functions
//

// Written by Edward Smith-Rowland based on numerous mathematics books.

#ifndef OLD_MODIFIED_BESSEL_FUNC_TCC
#define OLD_MODIFIED_BESSEL_FUNC_TCC 1

#include <emsr/fp_type_util.h>

namespace emsr
{
namespace detail
{

    ///
    ///
    ///
    template <typename Tp>
    void
    cyl_bessel_k_scaled_temme(const Tp nu, const Tp x,
                                Tp & K_nu, Tp & K_nup1, Tp & Kp_nu)
    {
      const Tp x2    = x / Tp(2);
      const Tp ln_x2 = std::log(x2);
      const Tp x2_nu = std::exp(nu * ln_x2);
      const Tp pi_nu   = _TR1_M_PI * nu;
      const Tp sigma   = -nu * ln_x2;
      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp sinrat  = (std::abs(pi_nu) < eps ? Tp(1) : pi_nu / std::sin(pi_nu));
      const Tp sinhrat = (std::abs(sigma) < eps ? Tp(1) : ::sinh(sigma) / sigma);
      const Tp ex = std::exp(x);

      Tp g_1pnu, g_1mnu, g1, g2;
      gamma_temme(nu, g_1pnu, g_1mnu, g1, g2);

      Tp fk = sinrat * (::cosh(sigma) * g1
                           - sinhrat * ln_x2 * g2);
      Tp pk = 0.5 / x2_nu * g_1pnu;
      Tp qk = 0.5 * x2_nu * g_1mnu;
      Tp hk = pk;
      Tp ck = Tp(1);
      Tp sum0 = fk;
      Tp sum1 = hk;

      const unsigned int max_iter = 15000;
      unsigned int k = 0;
      while(k < max_iter)
        {
          ++k;
          fk  = (k * fk + pk + qk) / (k * k - nu * nu);
          ck *= x2 * x2 / k;
          pk /= (k - nu);
          qk /= (k + nu);
          hk  = -k * fk + pk;
          const Tp del0 = ck * fk;
          const Tp del1 = ck * hk;
          sum0 += del0;
          sum1 += del1;
          if (std::abs(del0) < 0.5L * std::abs(sum0) * eps)
            break;
        }
  
      K_nu   = sum0 * ex;
      K_nup1 = sum1 * (Tp(2) / x) * ex;
      Kp_nu  = - K_nup1 + (nu / x) * K_nu;

      return;
    }


    ///
    ///  @brief  Return the A series used in the asymptotic expansions
    ///          of the Bessel functions.
    ///
    template<typename Tp>
    Tp
    cyl_mod_bessel_a_series(const Tp nu, const Tp x,
                              const bool alternate)
    {

      Tp sign = Tp(1);
      if (alternate)
        sign = Tp(-1);

      const Tp fnu2 = Tp(4) * nu * nu;
      const Tp xx = Tp(8) * x;

      Tp term = Tp(1);
      Tp sum = Tp(1);
      for (unsigned int i = 1; i < 20; ++i)
        {
          const Tp i2 = Tp(2 * i);
          const Tp i4 = Tp(4 * i);
          const Tp b = i4 - Tp(1);
          const Tp numer = sign * (fnu2 - b * b);
          const Tp denom = i2 * xx;
          term *= numer / denom;
          if (std::abs(term) < std::numeric_limits<Tp>::epsilon())
            {
              break;
            }
          sum += term;
        }

      return sum;
    }


    ///
    ///  @brief  Return the asymptotic solution of the Bessel function.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i_asymp(const Tp n, const Tp x)
    {

      Tp sum = cyl_mod_bessel_a_series(n, x, true);

      return std::exp(x) * sum / std::sqrt(2 * _TR1_M_PI * x);
    }


    ///
    /// x >> nu*nu+1
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i_scaled_asymp(const Tp nu, const Tp x)
    {
      const Tp mu = Tp(4) * nu * nu;
      const Tp mum1 = mu - Tp(1);
      const Tp mum9 = mu - Tp(9);
      const Tp pre = Tp(1) / std::sqrt(Tp(2) * _TR1_M_PI * x);
      const Tp r = mu / x;
      const Tp I = pre * (Tp(1) - mum1 / (Tp(8) * x)
                    + mum1*mum9 / (Tp(128) * x * x));

      return I;
    }


    //
    //  Debye functions Abramowitz & Stegun, 9.3.9-10
    //
    template<typename Tp>
    Tp 
    debye_u1(const std::vector<Tp> & tpow)
    {
      return (Tp(3UL) * tpow[1]
            - Tp(5UL) * tpow[3]) / Tp(24UL);
    }

    template<typename Tp>
    Tp 
    debye_u2(const std::vector<Tp> & tpow)
    {
      return (Tp(81UL) * tpow[2]
           - Tp(462UL) * tpow[4]
           + Tp(385UL) * tpow[6]) / Tp(1152UL);
    }

    template<typename Tp>
    Tp
    debye_u3(const std::vector<Tp> & tpow)
    {
      return (Tp(30375UL) * tpow[3]
           - Tp(369603UL) * tpow[5]
           + Tp(765765UL) * tpow[7]
           - Tp(425425UL) * tpow[9]) / Tp(414720UL);
    }

    template<typename Tp>
    Tp
    debye_u4(const std::vector<Tp> & tpow)
    {
      return (Tp(4465125UL) * tpow[4]
           - Tp(94121676UL) * tpow[6]
          + Tp(349922430UL) * tpow[8]
          - Tp(446185740UL) * tpow[10]
          + Tp(185910725UL) * tpow[12]) / Tp(39813120UL);
    }

    template<typename Tp>
    Tp
    debye_u5(const std::vector<Tp> & tpow)
    {
      return (Tp(1519035525UL) * tpow[5]
           - Tp(49286948607UL) * tpow[7]
          + Tp(284499769554UL) * tpow[9]
          - Tp(614135872350UL) * tpow[11]
          + Tp(566098157625UL) * tpow[13]
          - Tp(188699385875UL) * tpow[15]) / Tp(6688604160UL);
    }


#if ! _GLIBCXX_USE_C99_MATH_TR1
    template<typename Tp>
    Tp
    xhypot(const Tp x, const Tp y)
    {
      const Tp ax = std::abs(x);
      const Tp ay = std::abs(y);
      if (ay > ax)
        {
          const Tp rho = ax / ay;
          return ay * std::sqrt(Tp(1) + rho * rho);
        }
      else
        {
          const Tp rho = ay / ax;
          return ax * std::sqrt(Tp(1) + rho * rho);
        }
    }
#endif


    ///
    ///  @brief  Return the uniform asymptotic approximation for the
    ///          modofied Bessel function @f$ K_\nu(x) @f$.
    ///
    ///  nu -> Inf; uniform in x > 0.  (See Abramowitz & Stegun, 9.7.7)
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i_scaled_asymp_unif(const Tp nu, const Tp x)
    {
      const Tp z = x / nu;
#if _GLIBCXX_USE_C99_MATH_TR1
      const Tp hypot = std::tr1::hypot(Tp(1), z);
#else
      const Tp hypot = xhypot(Tp(1), z);
#endif
      const Tp pre = Tp(1) / std::sqrt(Tp(2) * _TR1_M_PI * nu * hypot);
      const Tp eta = hypot + std::log(z / (Tp(1) + hypot));
      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp root3_eps = std::exp(std::log(eps) / Tp(3));
      const Tp arg = ( z < Tp(1) / root3_eps
                          ? nu * (-z + eta)
                          : -Tp(0.5L) * nu / z * (Tp(1) - Tp(1) / (Tp(12) * z * z)) );
      const Tp ex = std::exp(arg);
      const Tp t = Tp(1) / hypot;
      std::vector<Tp> tpow;
      tpow.push_back(Tp(1));
      for (int i = 1; i < 16; ++i)
        tpow.push_back(t * tpow.back());
      Tp sum = Tp(1);
                + debye_u1(tpow) / (nu)
                + debye_u2(tpow) / (nu*nu)
                + debye_u3(tpow) / (nu*nu*nu)
                + debye_u4(tpow) / (nu*nu*nu*nu)
                + debye_u5(tpow) / (nu*nu*nu*nu*nu);
      const Tp I = pre * ex * sum;
      return I;
      //On exponent failure: Tp I = Tp(0);
      //return I;
    }


    ///
    ///  @brief  Return the uniform asymptotic approximation for the
    ///          modofied Bessel function @f$ K_\nu(x) @f$.
    ///
    ///  nu -> Inf; uniform in x > 0.  (See Abramowitz & Stegun, 9.7.8)
    ///
    template<typename Tp>
    Tp
    gsl_sf_bessel_k_scaled_asymp_unif(const Tp nu, const Tp x)
    {
      const Tp z = x / nu;
#if _GLIBCXX_USE_C99_MATH_TR1
      const Tp hypot = std::tr1::hypot(Tp(1), z);
#else
      const Tp hypot = xhypot(Tp(1), z);
#endif
      const Tp pre = std::sqrt(_TR1_M_PI / (Tp(2) * nu * hypot));
      const Tp eta = hypot + std::log(z / (Tp(1) + hypot));
      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp root3_eps = std::exp(std::log(eps) / Tp(3));
      const Tp arg = ( z < Tp(1) / root3_eps
                          ? nu * (z - eta)
                          : Tp(0.5L) * nu / z * (Tp(1) + Tp(1) / (Tp(12) * z * z)) );
      const Tp ex = std::exp(arg);
      const Tp t = Tp(1) / hypot;
      std::vector<Tp> tpow;
      tpow.push_back(Tp(1));
      for(int i = 1; i < 16; ++i)
        tpow.push_back(t * tpow.back());
      Tp fact = Tp(1);
      Tp sum = Tp(1)
                - debye_u1(tpow) / (nu)
                + debye_u2(tpow) / (nu*nu)
                - debye_u3(tpow) / (nu*nu*nu)
                + debye_u4(tpow) / (nu*nu*nu*nu)
                - debye_u5(tpow) / (nu*nu*nu*nu*nu);
      const Tp K = pre * ex * sum;
      return K;
    }


    ///
    ///  @brief  Evaluate the continued fraction CF1 for @f$ I_{nu+1}/I_nu @f$
    ///          using Gautschi (Euler) equivalent series.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i_cont_frac_1_series(const Tp nu, const Tp x)
    {
      const int maxk = 20000;
      Tp tk = Tp(1);
      Tp sum = Tp(1);
      Tp rhok = Tp(0);

      for (int k = 1; k < maxk; ++k)
        {
          Tp ak = (x / (nu + k))
                   * x / (Tp(4) * (nu + Tp(k + 1)));
          rhok = -ak * (Tp(1) + rhok)
                 / (Tp(1) + ak * (Tp(1) + rhok));
          tk *= rhok;
          sum += tk;
          if (std::abs(tk / sum) < std::numeric_limits<Tp>::epsilon())
            break;
        }

      Tp ratio = x / (Tp(2) * (nu + Tp(1))) * sum;

      //if (k == maxk)
      //  throw_logic;

      return ratio;
    }


    /// Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
    /// to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
    ///
    /// This is unstable for small x; x > 2 is a good cutoff.
    /// Also requires |nu| < 1/2.
    template<typename Tp>
    void
    cyl_bessel_k_scaled_steed_temme_cont_frac_2(const Tp nu,
                                                  const Tp x,
                                                  Tp & K_nu,
                                                  Tp & K_nup1,
                                                  Tp & Kp_nu)
    {
      const int max_iter = 10000;

      Tp bi = Tp(2) * (Tp(1) + x);
      Tp di = Tp(1) / bi;
      Tp delhi = di;
      Tp hi = di;

      Tp qi   = Tp(0);
      Tp qip1 = Tp(1);

      Tp ai = -(Tp(0.25L) - nu * nu);
      Tp a1 = ai;
      Tp ci = -ai;
      Tp Qi = -ai;

      Tp s = Tp(1) + Qi * delhi;

      for (int i = 2; i <= max_iter; ++i)
        {
          ai -= Tp(2 * (i - 1));
          ci  = -ai * ci / i;
          const Tp tmp = (qi - bi * qip1) / ai;
          qi = qip1;
          qip1 = tmp;
          Qi += ci * qip1;
          bi += Tp(2);
          di = Tp(1) / (bi + ai * di);
          delhi = (bi * di - Tp(1)) * delhi;
          hi += delhi;
          const Tp dels = Qi * delhi;
          s += dels;
          if (std::abs(dels / s) < std::numeric_limits<Tp>::epsilon())
            break;
        }

      hi *= -a1;

      K_nu = std::sqrt(_TR1_M_PI/(Tp(2) * x)) / s;
      K_nup1 = K_nu * (nu + x + Tp(0.5L) - hi) / x;
      Kp_nu = - K_nup1 + nu / x * K_nu;
      //if(i == max_iter)
      //  GSL_ERROR ("error", GSL_EMAXITER);

      return;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ I_{\nu} \f$.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i_scaled(const Tp nu, const Tp x)
    {

      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_i_scaled.");

      if (std::isnan(x))
        return -std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp(0))
        {
          if (nu == Tp(0))
            return Tp(1);
          else
            return Tp(0);
        }

      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp root3_eps = std::exp(std::log(eps) / Tp(3));

      if (x * x < Tp(10) * (nu + Tp(1)))
        {
          const Tp ex = exp(-x);
          const Tp b = cyl_bessel_ij_series(nu, x, Tp(1), 100);
          return b * ex;
        }
      else if (Tp(0.5L) / (nu * nu + x * x) < root3_eps)
        {
      //  return cyl_bessel_i_asymp(nu, x);
          return cyl_bessel_i_scaled_asymp_unif(nu, x);
        }
      else
        {
          int Nmu = static_cast<int>(nu + Tp(0.5L));
          Tp mu = nu - Nmu;  // -1/2 <= mu <= 1/2
          Tp K_mu, K_mup1, Kp_mu;
          Tp K_nu, K_nup1, K_num1;
          Tp I_nu_ratio;

          //  Obtain K_mu, K_mup1
          if (x < Tp(2))
            {
              cyl_bessel_k_scaled_temme(mu, x, K_mu, K_mup1, Kp_mu);
            }
          else
            {
              cyl_bessel_k_scaled_steed_temme_cont_frac_2(mu, x, K_mu, K_mup1, Kp_mu);
            }

          //  Recurse forward to obtain K_num1, K_nu
          K_nu   = K_mu;
          K_nup1 = K_mup1;

          for(int n = 0; n < Nmu; ++n)
            {
              K_num1 = K_nu;
              K_nu   = K_nup1;
              K_nup1 = Tp(2) * (mu + Tp(n + 1)) * K_nu / x + K_num1;
            }

          //  Calculate I_{nu+1}/I_nu
          I_nu_ratio = cyl_bessel_i_cont_frac_1_series(nu, x);

          // solve for I_nu
          Tp I_nu = Tp(1) / (x * (K_nup1 + I_nu_ratio * K_nu));
          I_nu = eps * (Tp(0.5L) * Nmu + Tp(2)) * std::abs(I_nu);
          return I_nu;
       }
    }


    ///
    ///
    ///
    template<typename Tp>
    Tp
    cyl_bessel_i(Tp nu, Tp x)
    {
      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<Tp>::quiet_NaN();

      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_i.");

      Tp I_scaled = cyl_bessel_i_scaled(nu, x);
      return I_scaled * std::exp(x);
    }


    ///
    ///  @brief  Return the asymptotic solution of the modified
    ///          Bessel function \f$ K_{\nu} \f$.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_k_asymp(const Tp n, const Tp x)
    {

      Tp sum = cyl_mod_bessel_a_series(n, x, false);

      return std::sqrt(_TR1_M_PI_2 / x) * std::exp(-x) * sum;
    }


    ///
    ///  x >> nu*nu + 1
    ///
    template<typename Tp>
    Tp
    cyl_bessel_k_scaled_asymp(const Tp nu, const Tp x)
    {
      const Tp mu = Tp(4) * nu * nu;
      const Tp mum1 = mu - Tp(1);
      const Tp mum9 = mu - Tp(9);
      const Tp pre = std::sqrt(_TR1_M_PI / (Tp(2) * x));
      const Tp r = nu / x;
      const Tp K = pre * (Tp(1) + mum1 / (Tp(8) * x) + mum1 * mum9 / (Tp(128) * x * x));

      return K;
    }


    ///
    ///  @brief  Return the scaled modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu} \f$.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_k_scaled(const Tp nu, const Tp x)
    {

      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<Tp>::quiet_NaN();

      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_k_scaled.");

      int Nmu = static_cast<int>(nu + Tp(0.5L));
      Tp mu = nu - Nmu;  // -1/2 <= mu <= 1/2
      Tp K_mu, K_mup1, Kp_mu;

      if (x < Tp(2))
        cyl_bessel_k_scaled_temme(mu, x, K_mu, K_mup1, Kp_mu);
      else
        cyl_bessel_k_scaled_steed_temme_cont_frac_2(mu, x, K_mu, K_mup1, Kp_mu);

      //  Recurse forward to obtain K_num1, K_nu.

      Tp K_nu   = K_mu;
      Tp K_nup1 = K_mup1;

      for (int n = 0; n < Nmu; ++n)
        {
          const Tp K_num1 = K_nu;
          K_nu   = K_nup1;
          K_nup1 = Tp(2) * (mu + n + 1) / x * K_nu + K_num1;
        }

      return K_nu;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu} \f$.
    ///
    template<typename Tp>
    Tp
    cyl_bessel_k(const Tp nu, const Tp x)
    {
      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<Tp>::quiet_NaN();

      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_k.");

      Tp K_scaled = cyl_bessel_k_scaled(nu, x);
      return K_scaled * std::exp(-x);
    }

} // namespace detail
} // namespace emsr

#endif // OLD_MODIFIED_BESSEL_FUNC_TCC
