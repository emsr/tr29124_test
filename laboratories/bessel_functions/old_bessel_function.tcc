// Special functions -*- C++ -*-

// Copyright (C) 2006-2019 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
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
    ///  @brief  Evaluate the Neumann functions of order @f$ N_\nu(x) @f$
    ///          and @f$ N_{\nu + 1}(x) @f$ by Temme's method, see
    ///          Temme, Journal of Computational Physics, vol 21, 343 (1976)
    ///
    ///  This method works best for @f$ |\nu| < 1/2 @f$ and @f$ x <= 2 @f$.
    ///
    template <typename _Tp>
    void
    cyl_neumann_temme(const _Tp nu, const _Tp x,
                        _Tp & Ynu, _Tp & Ynup1)
    {
      const _Tp x2 = x / _Tp(2);
      const _Tp ln_x2 = std::log(x2);
      const _Tp x2_nu = std::exp(nu * ln_x2);
      const _Tp pi_nu   = _TR1_M_PI * nu;
      const _Tp alpha   = pi_nu / _Tp(2);
      const _Tp sigma   = -nu * ln_x2;
      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp sinrat  = (std::abs(pi_nu) < eps
                            ? _Tp(1) : pi_nu / std::sin(pi_nu));
      const _Tp sinhrat = (std::abs(sigma) < eps
                            ? _Tp(1) : std::sinh(sigma) / sigma);
      const _Tp sinhalf = (std::abs(alpha) < eps
                            ? _Tp(1) : std::sin(alpha) / alpha);
      const _Tp sin_sqr = nu * _TR1_M_PI * _TR1_M_PI
                          * sinhalf * sinhalf / _Tp(2);

      _Tp g_1pnu, g_1mnu, g1, g2;
      gamma_temme(nu, g_1pnu, g_1mnu, g1, g2);

      _Tp fk = (_Tp(2) / _TR1_M_PI) * sinrat
               * (std::cosh(sigma) * g1 - sinhrat * ln_x2 * g2);
      _Tp pk = _Tp(1) / _TR1_M_PI / x2_nu * g_1pnu;
      _Tp qk = _Tp(1) / _TR1_M_PI * x2_nu * g_1mnu;
      _Tp hk = pk;
      _Tp ck = _Tp(1);

      _Tp sum0 = fk + sin_sqr * qk;
      _Tp sum1 = pk;

      const unsigned int max_iter = 15000;
      unsigned int k = 0;
      while (k < max_iter)
        {
          ++k;
          fk  = (k * fk + pk + qk) / (k * k - nu * nu);
          ck *= -x2 * x2 / k;
          pk /= (k - nu);
          qk /= (k + nu);
          const _Tp gk  = fk + sin_sqr * qk;
          hk  = -k * gk + pk; 
          const _Tp del0 = ck * gk;
          const _Tp del1 = ck * hk;
          sum0 += del0;
          sum1 += del1;
          if (std::abs(del0) < 0.5L * (_Tp(1) + std::abs(sum0)) * eps)
            break;
        }

      Ynu   = -sum0;
      Ynup1 = -sum1 * _Tp(2) / x;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Bessel function
    ///         of order nu: \f$J_{\nu}\f$.
    ///
    template <typename _Tp>
    _Tp
    cyl_bessel_j_asymp(const _Tp nu, const _Tp x)
    {
      const _Tp coef = std::sqrt(_Tp(2)/(_Tp(_TR1_M_PI) * x));
      const _Tp mu   = _Tp(4) * nu * nu;
      const _Tp mum1 = mu - _Tp(1);
      const _Tp mum9 = mu - _Tp(9);
      const _Tp mum25 = mu - _Tp(25);
      const _Tp P = _Tp(1) - mum1 * mum9 / (_Tp(128) * x * x);
      const _Tp Q = mum1/(_Tp(8) * x)
                    * (_Tp(1) - mum9 * mum25 / (_Tp(384) * x * x));

      const _Tp arg = x - (nu + _Tp(0.5L)) * _Tp(_TR1_M_PI_2);
      const _Tp c = std::cos(arg);
      const _Tp s = std::sin(arg);

      const _Tp Jnu = coef * (c * P - s * Q);

      return Jnu;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Bessel function
    ///         amplitude and phase functions of order nu: \f$ M_{\nu} \f$
    ///         and \f$ \theta_{\nu} \f$
    ///
    template <typename _Tp>
    void
    cyl_bessel_asymp_amp_phase(const _Tp nu, const _Tp x,
                                 _Tp & amp, _Tp & phase)
    {
      const _Tp r  = _Tp(2) * nu / x;
      const _Tp rr = r * r;
      const _Tp xx = x * x;

      const _Tp amp_term1 = (rr - _Tp(1) / xx) / _Tp(8);
      const _Tp amp_term2 = amp_term1
                            * (rr - _Tp(9) / xx) * _Tp(3) / _Tp(16);
      amp = (_Tp(2) / _Tp(_TR1_M_PI))
            * (_Tp(1) + amp_term1 + amp_term2);
      amp = std::sqrt(amp) / std::sqrt(x);

      const _Tp ph_term1 = x * (rr - _Tp(1) / xx) / _Tp(8);
      const _Tp ph_term2 = ph_term1
                           * (rr - _Tp(25) / xx) / _Tp(48);
      phase = -_Tp(_TR1_M_PI) / _Tp(4) + ph_term1 + ph_term2;

      return;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Neumann
    ///         function of order \f$ \nu \f$: \f$ N_{\nu} \f$.
    ///
    template <typename _Tp>
    _Tp
    cyl_neumann_asymp(const _Tp nu, const _Tp x)
    {
      _Tp amp, phase;
      cyl_bessel_asymp_amp_phase(nu, x, amp, phase);
      const _Tp sin_term = std::sin(x) * std::cos(phase)
                           + std::cos(x) * std::sin(phase);
      const _Tp N_nu = amp * sin_term;

      return N_nu;
    }


    ///
    ///  @brief This routine returns the cylindrical Bessel function
    ///         of order \f$ \nu \f$: \f$ J_{\nu} \f$ by series expansion.
    ///
    ///  See Abramowitz & Stegun, 9.1.10
    ///      Abramowitz & Stegun, 9.6.7
    ///
    template <typename _Tp>
    _Tp
    cyl_bessel_ij_series(const _Tp nu, const _Tp x, const _Tp sgn,
                           const unsigned int max_iter)
    {

      const _Tp x2 = x / _Tp(2);
      _Tp fact = nu * std::log(x2);
#if _GLIBCXX_USE_C99_MATH_TR1
      fact -= std::tr1::lgamma(nu + _Tp(1));
#else
      fact -= log_gamma(nu + _Tp(1));
#endif
      fact = std::exp(fact);
      const _Tp xx4 = sgn * x2 * x2;
      _Tp Jn = _Tp(1);
      _Tp term = _Tp(1);

      for (unsigned int i = 1; i < max_iter; ++i)
        {
          term *= xx4 / (_Tp(i) * (nu + _Tp(i)));
          Jn += term;
          if (std::abs(term / Jn) < std::numeric_limits<_Tp>::epsilon())
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
    template <typename _Tp>
    void
    bessel_j_cont_frac_1(const _Tp nu, const _Tp x,
                           _Tp & ratio, _Tp & sgn)
    {
      const _Tp big = std::numeric_limits<_Tp>::max();
      const int maxiter = 10000;
      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      int n = 1;
      _Tp Anm2 = _Tp(1);
      _Tp Bnm2 = _Tp(0);
      _Tp Anm1 = _Tp(0);
      _Tp Bnm1 = _Tp(1);
      _Tp a1 = x / (_Tp(2) * (nu + _Tp(1)));
      _Tp An = Anm1 + a1 * Anm2;
      _Tp Bn = Bnm1 + a1 * Bnm2;
      _Tp an;
      _Tp dn = a1;
      ratio = An / Bn;
      sgn  = _Tp(1);

      while (n < maxiter)
        {
          ++n;
          Anm2 = Anm1;
          Bnm2 = Bnm1;
          Anm1 = An;
          Bnm1 = Bn;
          an = -x * x / (_Tp(4)
               * (nu + _Tp(n - 1)) * (nu + _Tp(n)));
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

          const _Tp old_ratio = ratio;
          ratio = An / Bn;
          const _Tp del = old_ratio / ratio;

          dn = _Tp(1) / (_Tp(2) * (nu + n) / x - dn);
          if (dn < _Tp(0))
            sgn = -sgn;

          if (std::abs(del - _Tp(1)) < _Tp(2) * eps)
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
    template <typename _Tp>
    void
    bessel_jn_steed_cont_frac_2(const _Tp nu, const _Tp x,
                                  _Tp & p, _Tp & q)
    {
      const int max_iter = 10000;
      const _Tp _small = std::sqrt(std::numeric_limits<_Tp>::min());
      const _Tp eps = std::numeric_limits<_Tp>::epsilon();

      _Tp x_inv = _Tp(1) / x;
      _Tp a = _Tp(0.25L) - nu * nu;
      p = -x_inv / _Tp(2);
      q = _Tp(1);
      _Tp br = _Tp(2) * x;
      _Tp bi = _Tp(2);
      _Tp fact = a * x_inv / (p * p + q * q);
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
      for (int i = 2; i <= max_iter; ++i)
        {
          a  += _Tp(2 * (i - 1));
          bi += _Tp(2);
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
          if (std::abs(dlr - _Tp(1)) + std::abs(dli) < eps)
            break;
        }

      return;
    }


    ///
    ///  @brief  Return the A series used in the asymptotic expansions
    ///          of the Bessel functions.
    ///
    template <typename _Tp>
    std::pair<_Tp,_Tp>
    cyl_bessel_a_series(const _Tp nu, const _Tp x)
    {

      const _Tp fnu2 = _Tp(4) * nu * nu;
      const _Tp xx = _Tp(64) * x * x;

      _Tp eterm = _Tp(1);
      _Tp esum = _Tp(1);
      _Tp oterm = (fnu2 - _Tp(1)) / (_Tp(8) * x);
      _Tp osum = oterm;
      for (unsigned int i = 1; i < 10; ++i)
        {
          const _Tp i2 = _Tp(2 * i);
          const _Tp i4 = _Tp(4 * i);
          const _Tp a = i4 - _Tp(3);
          const _Tp b = i4 - _Tp(1);
          const _Tp c = i4 + _Tp(1);
          const _Tp enumer = -(fnu2 - a * a) * (fnu2 - b * b);
          const _Tp edenom = (i2 - _Tp(1)) * (i2) * xx;
          eterm *= enumer / edenom;
          const _Tp onumer = -(fnu2 - b * b) * (fnu2 - c * c);
          const _Tp odenom = (i2) * (i2 + _Tp(1)) * xx;
          oterm *= onumer / odenom;
          if (std::abs(eterm) < std::numeric_limits<_Tp>::epsilon()
           && std::abs(oterm) < std::numeric_limits<_Tp>::epsilon())
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
    template <typename _Tp>
    _Tp
    cyl_bessel_j_olver_asymp(const _Tp nu, const _Tp x)
    {

      std::pair<_Tp,_Tp> sum = cyl_bessel_a_series(nu, x);

      const _Tp arg = x - (nu + _Tp(0.5L)) * _Tp(_TR1_M_PI_2);
      _Tp Jn = std::sqrt(_Tp(2) / (_Tp(_TR1_M_PI) * x))
               * (std::cos(arg) * sum.first
                - std::sin(arg) * sum.second);

      return Jn;
    }


    ///
    ///  Return the cylindrical Bessel function of order \f$ \nu \f$:
    ///  \f$ J_{\nu} \f$.
    ///
    template <typename _Tp>
    _Tp
    cyl_bessel_j(const _Tp nu, const _Tp x)
    {

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_j.");

      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (x == _Tp(0))
        {
          if (nu == _Tp(0))
            return _Tp(1);
          else
            return _Tp(0);
        }

      if (x * x < _Tp(10) * (nu + _Tp(1)))
        return cyl_bessel_ij_series(nu, x, _Tp(-1), 100);
      else if (nu > _Tp(50))
        return cyl_bessel_j_olver_asymp(nu, x);
      else if (x > _Tp(1000))
        return cyl_bessel_j_asymp(nu, x);
      else
        {
          // -1/2 <= mu <= 1/2
          int Nmu = static_cast<int>(nu + 0.5L);
          _Tp mu = nu - _Tp(Nmu);

          //  Get J ratio.
          _Tp Jnup1_Jnu, sgn_Jnu;
          bessel_j_cont_frac_1(nu, x, Jnup1_Jnu, sgn_Jnu);

          if (x < _Tp(2))
            {
              // Determine Y_mu, Y_mup1 directly and recurse forward to nu.
              // Then use the CF1 information to solve for J_nu and J_nup1.
              _Tp Y_mu, Y_mup1;
              cyl_neumann_temme(mu, x, Y_mu, Y_mup1);

              _Tp Ynm1 = Y_mu;
              _Tp Yn   = Y_mup1;
              _Tp Ynp1 = _Tp(0);
              for (int n = 1; n < Nmu; ++n)
                {
                  Ynp1 = _Tp(2) * (mu + n) * Yn / x - Ynm1;
                  Ynm1 = Yn;
                  Yn = Ynp1;
                }

              const _Tp result = _Tp(2) / (_TR1_M_PI * x)
                                 / (Jnup1_Jnu * Yn - Ynp1);
              return result;
            }
          else
            {
              //  Recurse backward from nu to mu, determining the J ratio
              //  at mu. Use this together with a Steed method CF2 to
              //  determine the actual J_mu, and thus obtain the normalization.
              _Tp Jmu;
              _Tp Jmup1_Jmu;
              _Tp sgn_Jmu;
              _Tp Jmuprime_Jmu;
              _Tp P, Q;
              bessel_jn_steed_cont_frac_2(mu, x, P, Q);
              _Tp gamma;
              const _Tp sqrt_min = std::sqrt(std::numeric_limits<_Tp>::min());
 
              _Tp Jnp1 = sgn_Jnu * sqrt_min * Jnup1_Jnu;
              _Tp Jn   = sgn_Jnu * sqrt_min;
              _Tp Jnm1;
              for (int n = Nmu; n > 0; --n)
                {
                  Jnm1 = _Tp(2) * (mu + n) * Jn / x - Jnp1;
                  Jnp1 = Jn;
                  Jn = Jnm1;
                }
              Jmup1_Jmu = Jnp1 / Jn;
              sgn_Jmu = (Jn < _Tp(0) ? _Tp(-1) : _Tp(1));
              Jmuprime_Jmu = mu / x - Jmup1_Jmu;

              gamma = (P - Jmuprime_Jmu) / Q;
              Jmu = sgn_Jmu
                  * std::sqrt(2 / (_TR1_M_PI * x) / (Q + gamma * (P - Jmuprime_Jmu)));

              const _Tp result = Jmu * (sgn_Jnu * sqrt_min) / Jn;

              return result;
            }
        }
    }


    ///
    ///  @brief This routine returns the cylindrical Neumann function
    ///         of order \f$ \nu \f$: \f$ N_{\nu} \f$
    ///         by asymptotic expansion.  Olver p ~ 130.
    ///
    template <typename _Tp>
    _Tp
    cyl_neumann_n_olver_asymp(const _Tp nu, const _Tp x)
    {

      std::pair<_Tp,_Tp> sum = cyl_bessel_a_series(nu, x);

      const _Tp arg = x - (nu + _Tp(0.5L)) * _Tp(_TR1_M_PI_2);
      const _Tp N_nu = std::sqrt(_Tp(2) / (_Tp(_TR1_M_PI) * x))
                       * (std::sin(arg) * sum.first
                        + std::cos(arg) * sum.second);

      return N_nu;
    }


    ///
    ///  Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
    ///
    template <typename _Tp>
    void
    bessel_jn_mu_restricted(const _Tp mu, const _Tp x,
                              _Tp & Jmu, _Tp & Jmup1,
                              _Tp & Ymu, _Tp & Ymup1)
    {

      if (x < _Tp(0) || std::abs(mu) > _Tp(0.5L))
        {
          Jmu = _Tp(0);
          Jmup1 = _Tp(0);
          Ymu = _Tp(0);
          Ymup1 = _Tp(0);
        }
      else if (x == _Tp(0))
        {
          if (mu == _Tp(0))
            Jmu = _Tp(1);
          else
            Jmu = _Tp(0);
          Jmup1 = _Tp(0);
          Ymu = _Tp(0);
          Ymup1 = _Tp(0);
        }
      else {

        if (x < _Tp(2))
          {
            //  Use Taylor series for J and the Temme series for Y.
            //  The Taylor series for J requires nu > 0, so we shift
            //  up one and use the recursion relation to get Jmu, in
            //  case mu < 0.
            _Tp Jmup2;
            Jmup1 = cyl_bessel_ij_series(mu + _Tp(1), x, _Tp(-1), 100);
            Jmup2 = cyl_bessel_ij_series(mu + _Tp(2), x, _Tp(-1), 100);
            Jmu  = _Tp(2) * (mu + _Tp(1)) * Jmup1 / x - Jmup2;
            cyl_neumann_temme(mu, x, Ymu, Ymup1);
          }
        else if (x < _Tp(1000))
          {
            _Tp J_ratio, J_sgn;
            bessel_j_cont_frac_1(mu, x, J_ratio, J_sgn);
            _Tp P, Q;
            bessel_jn_steed_cont_frac_2(mu, x, P, Q);
            _Tp Jprime_J_ratio = mu / x - J_ratio;
            _Tp gamma = (P - Jprime_J_ratio) / Q;
            Jmu = J_sgn * std::sqrt(_Tp(2) / (_TR1_M_PI * x)
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
            Jmup1 = cyl_bessel_j_asymp(mu + _Tp(1), x);
            Ymu = cyl_neumann_asymp(mu, x);
            Ymup1 = cyl_neumann_asymp(mu + _Tp(1), x);
            return;
          }
      }
    }


    ///
    ///  Return the cylindrical Neumann function of order nu: \f$ N_{\nu}(x) \f$.
    ///
    template <typename _Tp>
    _Tp
    cyl_neumann_n(const _Tp nu, const _Tp x)
    {

      if (std::isnan(nu) || std::isnan(x))
        return - std::numeric_limits<_Tp>::quiet_NaN();

      if (x == _Tp(0))
        return std::numeric_limits<_Tp>::infinity();

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_neumann_n.");

      if (nu > _Tp(50))
        {
          return cyl_neumann_n_olver_asymp(nu, x);
        }
      else
        {
          int Nmu = static_cast<int>(nu + _Tp(0.5L));
          // -1/2 <= mu <= 1/2
          _Tp mu = nu - Nmu;

          _Tp Y_mu, Y_mup1;

          if (x < _Tp(2))
            {
              //  Determine Ymu, Ymup1 directly.
              cyl_neumann_temme(mu, x, Y_mu, Y_mup1);
            }
          else
            {
              //  Determine Ymu, Ymup1 and Jmu, Jmup1.
              _Tp J_mu, J_mup1;
              bessel_jn_mu_restricted(mu, x, J_mu, J_mup1, Y_mu, Y_mup1);
            }

          //  Forward recursion to get Ynu, Ynup1.
          _Tp Ynm1 = Y_mu;
          _Tp Yn   = Y_mup1;
          for (int n = 1; n <= Nmu; ++n)
            {
              _Tp Ynp1 = _Tp(2) * (mu + n) * Yn / x - Ynm1;
              Ynm1 = Yn;
              Yn   = Ynp1;
            }

          _Tp result  = Ynm1; // Y_nu

          return result;
        }
    }


    ///
    ///  Return the spherical Bessel function.
    ///
    template <typename _Tp>
    _Tp
    sph_bessel(const unsigned int n, const _Tp x)
    {

      if (std::isnan(x))
        return std::numeric_limits<_Tp>::quiet_NaN();

      if (x == _Tp(0))
        return std::numeric_limits<_Tp>::infinity();  //  TODO: This is wrong?

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in sph_bessel.");

      return std::sqrt(_Tp(_TR1_M_PI_2) / x)
           * cyl_bessel_j(_Tp(n) + _Tp(0.5L), x);
    }


    ///
    ///  Return the spherical Neumann function.
    ///
    template <typename _Tp>
    _Tp
    sph_neumann(const unsigned int n, const _Tp x)
    {

      if (std::isnan(x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (x == _Tp(0))
        return std::numeric_limits<_Tp>::infinity();  //  TODO: This is wrong?

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in sph_neumann.");

      return std::sqrt(_Tp(_TR1_M_PI_2) / x)
           * cyl_neumann_n(_Tp(n) + _Tp(0.5L), x);
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _TR1_BESSEL_FUNCTION_TCC
