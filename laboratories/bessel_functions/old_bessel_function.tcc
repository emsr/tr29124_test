
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

/** @file old_bessel_function.tcc
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

#ifndef OLD_BESSEL_FUNCTION_TCC
#define OLD_BESSEL_FUNCTION_TCC 1

#include <emsr/fp_type_util.h>

namespace emsr
{
namespace detail
{

    ///
    ///  @brief  Evaluate the Neumann functions of order @f$ N_\nu(x) @f$
    ///          and @f$ N_{\nu + 1}(x) @f$ by Temme's method, see
    ///          Temme, Journal of Computational Physics, vol 21, 343 (1976)
    ///
    ///  This method works best for @f$ |\nu| < 1/2 @f$ and @f$ x <= 2 @f$.
    ///
    template <typename Tp>
    void
    cyl_neumann_temme(const Tp nu, const Tp x,
                        Tp & Ynu, Tp & Ynup1)
    {
      const Tp x2 = x / Tp(2);
      const Tp ln_x2 = std::log(x2);
      const Tp x2_nu = std::exp(nu * ln_x2);
      const Tp pi_nu   = _TR1_M_PI * nu;
      const Tp alpha   = pi_nu / Tp(2);
      const Tp sigma   = -nu * ln_x2;
      const Tp eps = std::numeric_limits<Tp>::epsilon();
      const Tp sinrat  = (std::abs(pi_nu) < eps
                            ? Tp(1) : pi_nu / std::sin(pi_nu));
      const Tp sinhrat = (std::abs(sigma) < eps
                            ? Tp(1) : std::sinh(sigma) / sigma);
      const Tp sinhalf = (std::abs(alpha) < eps
                            ? Tp(1) : std::sin(alpha) / alpha);
      const Tp sin_sqr = nu * _TR1_M_PI * _TR1_M_PI
                          * sinhalf * sinhalf / Tp(2);

      Tp g_1pnu, g_1mnu, g1, g2;
      gamma_temme(nu, g_1pnu, g_1mnu, g1, g2);

      Tp fk = (Tp(2) / _TR1_M_PI) * sinrat
               * (std::cosh(sigma) * g1 - sinhrat * ln_x2 * g2);
      Tp pk = Tp(1) / _TR1_M_PI / x2_nu * g_1pnu;
      Tp qk = Tp(1) / _TR1_M_PI * x2_nu * g_1mnu;
      Tp hk = pk;
      Tp ck = Tp(1);

      Tp sum0 = fk + sin_sqr * qk;
      Tp sum1 = pk;

      const unsigned int max_iter = 15000;
      unsigned int k = 0;
      while (k < max_iter)
        {
          ++k;
          fk  = (k * fk + pk + qk) / (k * k - nu * nu);
          ck *= -x2 * x2 / k;
          pk /= (k - nu);
          qk /= (k + nu);
          const Tp gk  = fk + sin_sqr * qk;
          hk  = -k * gk + pk; 
          const Tp del0 = ck * gk;
          const Tp del1 = ck * hk;
          sum0 += del0;
          sum1 += del1;
          if (std::abs(del0) < 0.5L * (Tp(1) + std::abs(sum0)) * eps)
            break;
        }

      Ynu   = -sum0;
      Ynup1 = -sum1 * Tp(2) / x;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Bessel function
    ///         of order nu: \f$J_{\nu}\f$.
    ///
    template <typename Tp>
    Tp
    cyl_bessel_j_asymp(const Tp nu, const Tp x)
    {
      const Tp coef = std::sqrt(Tp(2)/(Tp(_TR1_M_PI) * x));
      const Tp mu   = Tp(4) * nu * nu;
      const Tp mum1 = mu - Tp(1);
      const Tp mum9 = mu - Tp(9);
      const Tp mum25 = mu - Tp(25);
      const Tp P = Tp(1) - mum1 * mum9 / (Tp(128) * x * x);
      const Tp Q = mum1/(Tp(8) * x)
                    * (Tp(1) - mum9 * mum25 / (Tp(384) * x * x));

      const Tp arg = x - (nu + Tp(0.5L)) * Tp(_TR1_M_PI_2);
      const Tp c = std::cos(arg);
      const Tp s = std::sin(arg);

      const Tp Jnu = coef * (c * P - s * Q);

      return Jnu;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Bessel function
    ///         amplitude and phase functions of order nu: \f$ M_{\nu} \f$
    ///         and \f$ \theta_{\nu} \f$
    ///
    template <typename Tp>
    void
    cyl_bessel_asymp_amp_phase(const Tp nu, const Tp x,
                                 Tp & amp, Tp & phase)
    {
      const Tp r  = Tp(2) * nu / x;
      const Tp rr = r * r;
      const Tp xx = x * x;

      const Tp amp_term1 = (rr - Tp(1) / xx) / Tp(8);
      const Tp amp_term2 = amp_term1
                            * (rr - Tp(9) / xx) * Tp(3) / Tp(16);
      amp = (Tp(2) / Tp(_TR1_M_PI))
            * (Tp(1) + amp_term1 + amp_term2);
      amp = std::sqrt(amp) / std::sqrt(x);

      const Tp ph_term1 = x * (rr - Tp(1) / xx) / Tp(8);
      const Tp ph_term2 = ph_term1
                           * (rr - Tp(25) / xx) / Tp(48);
      phase = -Tp(_TR1_M_PI) / Tp(4) + ph_term1 + ph_term2;

      return;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Neumann
    ///         function of order \f$ \nu \f$: \f$ N_{\nu} \f$.
    ///
    template <typename Tp>
    Tp
    cyl_neumann_asymp(const Tp nu, const Tp x)
    {
      Tp amp, phase;
      cyl_bessel_asymp_amp_phase(nu, x, amp, phase);
      const Tp sin_term = std::sin(x) * std::cos(phase)
                           + std::cos(x) * std::sin(phase);
      const Tp N_nu = amp * sin_term;

      return N_nu;
    }


    ///
    ///  @brief This routine returns the cylindrical Bessel function
    ///         of order \f$ \nu \f$: \f$ J_{\nu} \f$ by series expansion.
    ///
    ///  See Abramowitz & Stegun, 9.1.10
    ///      Abramowitz & Stegun, 9.6.7
    ///
    template <typename Tp>
    Tp
    cyl_bessel_ij_series(const Tp nu, const Tp x, const Tp sgn,
                           const unsigned int max_iter)
    {

      const Tp x2 = x / Tp(2);
      Tp fact = nu * std::log(x2);
#if _GLIBCXX_USE_C99_MATH_TR1
      fact -= std::tr1::lgamma(nu + Tp(1));
#else
      fact -= log_gamma(nu + Tp(1));
#endif
      fact = std::exp(fact);
      const Tp xx4 = sgn * x2 * x2;
      Tp Jn = Tp(1);
      Tp term = Tp(1);

      for (unsigned int i = 1; i < max_iter; ++i)
        {
          term *= xx4 / (Tp(i) * (nu + Tp(i)));
          Jn += term;
          if (std::abs(term / Jn) < std::numeric_limits<Tp>::epsilon())
            break;
        }

      return fact * Jn;
    }


    ///
    ///  @brief Find the logarithmic derivative of the Bessel function
    ///         of order @f$ \nu @f$ by a continued fraction method.
    ///  @f[
    ///   f_\nu(x) = \frac{J'_\nu(x)}{J_\nu(x)}
    ///            = \frac{\nu}{x} - \frac{J_{\nu+1}(x)}{J_\nu(x)}
    ///  @f]
    ///
    ///  See Numerical Recipes in C, 2nd Ed.,
    ///      W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
    ///      Cambridge University Press, 1992, pp 242-243
    ///
    template <typename Tp>
    void
    bessel_j_cont_frac_1(const Tp nu, const Tp x,
                           Tp & ratio, Tp & sgn)
    {
      const Tp big = std::numeric_limits<Tp>::max();
      const int maxiter = 10000;
      const Tp eps = std::numeric_limits<Tp>::epsilon();
      int n = 1;
      Tp Anm2 = Tp(1);
      Tp Bnm2 = Tp(0);
      Tp Anm1 = Tp(0);
      Tp Bnm1 = Tp(1);
      Tp a1 = x / (Tp(2) * (nu + Tp(1)));
      Tp An = Anm1 + a1 * Anm2;
      Tp Bn = Bnm1 + a1 * Bnm2;
      Tp an;
      Tp dn = a1;
      ratio = An / Bn;
      sgn  = Tp(1);

      while (n < maxiter)
        {
          ++n;
          Anm2 = Anm1;
          Bnm2 = Bnm1;
          Anm1 = An;
          Bnm1 = Bn;
          an = -x * x / (Tp(4)
               * (nu + Tp(n - 1)) * (nu + Tp(n)));
          An = Anm1 + an * Anm2;
          Bn = Bnm1 + an * Bnm2;

          if (std::abs(An) > big || std::abs(Bn) > big)
            {
              An /= big;
              Bn /= big;
              Anm1 /= big;
              Bnm1 /= big;
              Anm2 /= big;
              Bnm2 /= big;
            }

          const Tp old_ratio = ratio;
          ratio = An / Bn;
          const Tp del = old_ratio / ratio;

          dn = Tp(1) / (Tp(2) * (nu + n) / x - dn);
          if (dn < Tp(0))
            sgn = -sgn;

          if (std::abs(del - Tp(1)) < Tp(2) * eps)
            break;
      }

      return;
    }


    ///
    ///  @brief  Compute the logarithmic derivative of the Hankel function
    ///          of the first kind by Steed's continued fraction method.
    ///  @f[
    ///   p_\nu(x) + iq_\nu(x)
    ///            = \frac{H^{(1)}'_\nu(x)}{H^{(1)}_\nu(x)}
    ///            = \frac{J'_\nu(x) + iN'\nu(x)}{J_\nu(x) + iN\nu(x)}
    ///            = \frac{\nu}{x} - \frac{J_{\nu+1}(x)}{J_\nu(x)}
    ///  @f]
    ///
    template <typename Tp>
    void
    bessel_jn_steed_cont_frac_2(const Tp nu, const Tp x,
                                  Tp & p, Tp & q)
    {
      const int max_iter = 10000;
      const Tp _small = std::sqrt(std::numeric_limits<Tp>::min());
      const Tp eps = std::numeric_limits<Tp>::epsilon();

      Tp x_inv = Tp(1) / x;
      Tp a = Tp(0.25L) - nu * nu;
      p = -x_inv / Tp(2);
      q = Tp(1);
      Tp br = Tp(2) * x;
      Tp bi = Tp(2);
      Tp fact = a * x_inv / (p * p + q * q);
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
      for (int i = 2; i <= max_iter; ++i)
        {
          a  += Tp(2 * (i - 1));
          bi += Tp(2);
          dr = a * dr + br;
          di = a * di + bi;
          if (std::abs(dr) + std::abs(di) < _small)
            dr = _small;
          fact = a / (cr * cr + ci * ci);
          cr = br + cr * fact;
          ci = bi - ci * fact;
          if (std::abs(cr) + std::abs(ci) < _small)
            cr = _small;
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

      return;
    }


    ///
    ///  @brief  Return the A series used in the asymptotic expansions
    ///          of the Bessel functions.
    ///
    template <typename Tp>
    std::pair<Tp,Tp>
    cyl_bessel_a_series(const Tp nu, const Tp x)
    {

      const Tp fnu2 = Tp(4) * nu * nu;
      const Tp xx = Tp(64) * x * x;

      Tp eterm = Tp(1);
      Tp esum = Tp(1);
      Tp oterm = (fnu2 - Tp(1)) / (Tp(8) * x);
      Tp osum = oterm;
      for (unsigned int i = 1; i < 10; ++i)
        {
          const Tp i2 = Tp(2 * i);
          const Tp i4 = Tp(4 * i);
          const Tp a = i4 - Tp(3);
          const Tp b = i4 - Tp(1);
          const Tp c = i4 + Tp(1);
          const Tp enumer = -(fnu2 - a * a) * (fnu2 - b * b);
          const Tp edenom = (i2 - Tp(1)) * (i2) * xx;
          eterm *= enumer / edenom;
          const Tp onumer = -(fnu2 - b * b) * (fnu2 - c * c);
          const Tp odenom = (i2) * (i2 + Tp(1)) * xx;
          oterm *= onumer / odenom;
          if (std::abs(eterm) < std::numeric_limits<Tp>::epsilon()
           && std::abs(oterm) < std::numeric_limits<Tp>::epsilon())
            {
              break;
            }
          esum += eterm;
          osum += oterm;
        }

      return std::make_pair(esum, osum);
    }


    ///
    ///  @brief This routine returns the cylindrical Bessel function
    ///         of order \f$ \nu \f$: \f$ J_{\nu} \f$
    ///         by asymptotic expansion.  Olver p ~ 130.
    ///
    template <typename Tp>
    Tp
    cyl_bessel_j_olver_asymp(const Tp nu, const Tp x)
    {

      std::pair<Tp,Tp> sum = cyl_bessel_a_series(nu, x);

      const Tp arg = x - (nu + Tp(0.5L)) * Tp(_TR1_M_PI_2);
      Tp Jn = std::sqrt(Tp(2) / (Tp(_TR1_M_PI) * x))
               * (std::cos(arg) * sum.first
                - std::sin(arg) * sum.second);

      return Jn;
    }


    ///
    ///  Return the cylindrical Bessel function of order \f$ \nu \f$:
    ///  \f$ J_{\nu} \f$.
    ///
    template <typename Tp>
    Tp
    cyl_bessel_j(const Tp nu, const Tp x)
    {

      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_j.");

      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp(0))
        {
          if (nu == Tp(0))
            return Tp(1);
          else
            return Tp(0);
        }

      if (x * x < Tp(10) * (nu + Tp(1)))
        return cyl_bessel_ij_series(nu, x, Tp(-1), 100);
      else if (nu > Tp(50))
        return cyl_bessel_j_olver_asymp(nu, x);
      else if (x > Tp(1000))
        return cyl_bessel_j_asymp(nu, x);
      else
        {
          // -1/2 <= mu <= 1/2
          int Nmu = static_cast<int>(nu + 0.5L);
          Tp mu = nu - Tp(Nmu);

          //  Get J ratio.
          Tp Jnup1_Jnu, sgn_Jnu;
          bessel_j_cont_frac_1(nu, x, Jnup1_Jnu, sgn_Jnu);

          if (x < Tp(2))
            {
              // Determine Y_mu, Y_mup1 directly and recurse forward to nu.
              // Then use the CF1 information to solve for J_nu and J_nup1.
              Tp Y_mu, Y_mup1;
              cyl_neumann_temme(mu, x, Y_mu, Y_mup1);

              Tp Ynm1 = Y_mu;
              Tp Yn   = Y_mup1;
              Tp Ynp1 = Tp(0);
              for (int n = 1; n < Nmu; ++n)
                {
                  Ynp1 = Tp(2) * (mu + n) * Yn / x - Ynm1;
                  Ynm1 = Yn;
                  Yn = Ynp1;
                }

              const Tp result = Tp(2) / (_TR1_M_PI * x)
                                 / (Jnup1_Jnu * Yn - Ynp1);
              return result;
            }
          else
            {
              //  Recurse backward from nu to mu, determining the J ratio
              //  at mu. Use this together with a Steed method CF2 to
              //  determine the actual J_mu, and thus obtain the normalization.
              Tp Jmu;
              Tp Jmup1_Jmu;
              Tp sgn_Jmu;
              Tp Jmuprime_Jmu;
              Tp P, Q;
              bessel_jn_steed_cont_frac_2(mu, x, P, Q);
              Tp gamma;
              const Tp sqrt_min = std::sqrt(std::numeric_limits<Tp>::min());
 
              Tp Jnp1 = sgn_Jnu * sqrt_min * Jnup1_Jnu;
              Tp Jn   = sgn_Jnu * sqrt_min;
              Tp Jnm1;
              for (int n = Nmu; n > 0; --n)
                {
                  Jnm1 = Tp(2) * (mu + n) * Jn / x - Jnp1;
                  Jnp1 = Jn;
                  Jn = Jnm1;
                }
              Jmup1_Jmu = Jnp1 / Jn;
              sgn_Jmu = (Jn < Tp(0) ? Tp(-1) : Tp(1));
              Jmuprime_Jmu = mu / x - Jmup1_Jmu;

              gamma = (P - Jmuprime_Jmu) / Q;
              Jmu = sgn_Jmu
                  * std::sqrt(2 / (_TR1_M_PI * x) / (Q + gamma * (P - Jmuprime_Jmu)));

              const Tp result = Jmu * (sgn_Jnu * sqrt_min) / Jn;

              return result;
            }
        }
    }


    ///
    ///  @brief This routine returns the cylindrical Neumann function
    ///         of order \f$ \nu \f$: \f$ N_{\nu} \f$
    ///         by asymptotic expansion.  Olver p ~ 130.
    ///
    template <typename Tp>
    Tp
    cyl_neumann_n_olver_asymp(const Tp nu, const Tp x)
    {

      std::pair<Tp,Tp> sum = cyl_bessel_a_series(nu, x);

      const Tp arg = x - (nu + Tp(0.5L)) * Tp(_TR1_M_PI_2);
      const Tp N_nu = std::sqrt(Tp(2) / (Tp(_TR1_M_PI) * x))
                       * (std::sin(arg) * sum.first
                        + std::cos(arg) * sum.second);

      return N_nu;
    }


    ///
    ///  Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
    ///
    template <typename Tp>
    void
    bessel_jn_mu_restricted(const Tp mu, const Tp x,
                              Tp & Jmu, Tp & Jmup1,
                              Tp & Ymu, Tp & Ymup1)
    {

      if (x < Tp(0) || std::abs(mu) > Tp(0.5L))
        {
          Jmu = Tp(0);
          Jmup1 = Tp(0);
          Ymu = Tp(0);
          Ymup1 = Tp(0);
        }
      else if (x == Tp(0))
        {
          if (mu == Tp(0))
            Jmu = Tp(1);
          else
            Jmu = Tp(0);
          Jmup1 = Tp(0);
          Ymu = Tp(0);
          Ymup1 = Tp(0);
        }
      else {

        if (x < Tp(2))
          {
            //  Use Taylor series for J and the Temme series for Y.
            //  The Taylor series for J requires nu > 0, so we shift
            //  up one and use the recursion relation to get Jmu, in
            //  case mu < 0.
            Tp Jmup2;
            Jmup1 = cyl_bessel_ij_series(mu + Tp(1), x, Tp(-1), 100);
            Jmup2 = cyl_bessel_ij_series(mu + Tp(2), x, Tp(-1), 100);
            Jmu  = Tp(2) * (mu + Tp(1)) * Jmup1 / x - Jmup2;
            cyl_neumann_temme(mu, x, Ymu, Ymup1);
          }
        else if (x < Tp(1000))
          {
            Tp J_ratio, J_sgn;
            bessel_j_cont_frac_1(mu, x, J_ratio, J_sgn);
            Tp P, Q;
            bessel_jn_steed_cont_frac_2(mu, x, P, Q);
            Tp Jprime_J_ratio = mu / x - J_ratio;
            Tp gamma = (P - Jprime_J_ratio) / Q;
            Jmu = J_sgn * std::sqrt(Tp(2) / (_TR1_M_PI * x)
                  / (Q + gamma * (P - Jprime_J_ratio)));
            Jmup1 = J_ratio * Jmu;
            Ymu = gamma * Jmu;
            Ymup1 = Ymu * (mu / x - P - Q / gamma);
            return;
          }
        else
          {
            //  Use asymptotics for large argument.
            Jmu = cyl_bessel_j_asymp(mu, x );
            Jmup1 = cyl_bessel_j_asymp(mu + Tp(1), x);
            Ymu = cyl_neumann_asymp(mu, x);
            Ymup1 = cyl_neumann_asymp(mu + Tp(1), x);
            return;
          }
      }
    }


    ///
    ///  Return the cylindrical Neumann function of order nu: \f$ N_{\nu}(x) \f$.
    ///
    template <typename Tp>
    Tp
    cyl_neumann_n(const Tp nu, const Tp x)
    {

      if (std::isnan(nu) || std::isnan(x))
        return - std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp(0))
        return std::numeric_limits<Tp>::infinity();

      if (x < Tp(0))
        throw std::domain_error("Bad argument in cyl_neumann_n.");

      if (nu > Tp(50))
        {
          return cyl_neumann_n_olver_asymp(nu, x);
        }
      else
        {
          int Nmu = static_cast<int>(nu + Tp(0.5L));
          // -1/2 <= mu <= 1/2
          Tp mu = nu - Nmu;

          Tp Y_mu, Y_mup1;

          if (x < Tp(2))
            {
              //  Determine Ymu, Ymup1 directly.
              cyl_neumann_temme(mu, x, Y_mu, Y_mup1);
            }
          else
            {
              //  Determine Ymu, Ymup1 and Jmu, Jmup1.
              Tp J_mu, J_mup1;
              bessel_jn_mu_restricted(mu, x, J_mu, J_mup1, Y_mu, Y_mup1);
            }

          //  Forward recursion to get Ynu, Ynup1.
          Tp Ynm1 = Y_mu;
          Tp Yn   = Y_mup1;
          for (int n = 1; n <= Nmu; ++n)
            {
              Tp Ynp1 = Tp(2) * (mu + n) * Yn / x - Ynm1;
              Ynm1 = Yn;
              Yn   = Ynp1;
            }

          Tp result  = Ynm1; // Y_nu

          return result;
        }
    }


    ///
    ///  Return the spherical Bessel function.
    ///
    template <typename Tp>
    Tp
    sph_bessel(const unsigned int n, const Tp x)
    {

      if (std::isnan(x))
        return std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp(0))
        return std::numeric_limits<Tp>::infinity();  //  TODO: This is wrong?

      if (x < Tp(0))
        throw std::domain_error("Bad argument in sph_bessel.");

      return std::sqrt(Tp(_TR1_M_PI_2) / x)
           * cyl_bessel_j(Tp(n) + Tp(0.5L), x);
    }


    ///
    ///  Return the spherical Neumann function.
    ///
    template <typename Tp>
    Tp
    sph_neumann(const unsigned int n, const Tp x)
    {

      if (std::isnan(x))
        return -std::numeric_limits<Tp>::quiet_NaN();

      if (x == Tp(0))
        return std::numeric_limits<Tp>::infinity();  //  TODO: This is wrong?

      if (x < Tp(0))
        throw std::domain_error("Bad argument in sph_neumann.");

      return std::sqrt(Tp(_TR1_M_PI_2) / x)
           * cyl_neumann_n(Tp(n) + Tp(0.5L), x);
    }

} // namespace detail
} // namespace emsr

#endif // OLD_BESSEL_FUNCTION_TCC
