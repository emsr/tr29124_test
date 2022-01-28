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

/** @file tr1/sf_mod_bessel.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

//
// ISO C++ 14882 TR1: 5.2  Special functions
//

// Written by Edward Smith-Rowland based on numerous mathematics books.

#ifndef _TR1_MODIFIED_BESSEL_FUNC_TCC
#define _TR1_MODIFIED_BESSEL_FUNC_TCC 1

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
    cyl_bessel_k_scaled_temme(const _Tp nu, const _Tp x,
                                _Tp & K_nu, _Tp & K_nup1, _Tp & Kp_nu)
    {
      const _Tp x2    = x / _Tp(2);
      const _Tp ln_x2 = std::log(x2);
      const _Tp x2_nu = std::exp(nu * ln_x2);
      const _Tp pi_nu   = _TR1_M_PI * nu;
      const _Tp sigma   = -nu * ln_x2;
      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp sinrat  = (std::abs(pi_nu) < eps ? _Tp(1) : pi_nu / std::sin(pi_nu));
      const _Tp sinhrat = (std::abs(sigma) < eps ? _Tp(1) : ::sinh(sigma) / sigma);
      const _Tp ex = std::exp(x);

      _Tp g_1pnu, g_1mnu, g1, g2;
      gamma_temme(nu, g_1pnu, g_1mnu, g1, g2);

      _Tp fk = sinrat * (::cosh(sigma) * g1
                           - sinhrat * ln_x2 * g2);
      _Tp pk = 0.5 / x2_nu * g_1pnu;
      _Tp qk = 0.5 * x2_nu * g_1mnu;
      _Tp hk = pk;
      _Tp ck = _Tp(1);
      _Tp sum0 = fk;
      _Tp sum1 = hk;

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
          const _Tp del0 = ck * fk;
          const _Tp del1 = ck * hk;
          sum0 += del0;
          sum1 += del1;
          if (std::abs(del0) < 0.5L * std::abs(sum0) * eps)
            break;
        }
  
      K_nu   = sum0 * ex;
      K_nup1 = sum1 * (_Tp(2) / x) * ex;
      Kp_nu  = - K_nup1 + (nu / x) * K_nu;

      return;
    }


    ///
    ///  @brief  Return the A series used in the asymptotic expansions
    ///          of the Bessel functions.
    ///
    template<typename _Tp>
    _Tp
    cyl_mod_bessel_a_series(const _Tp nu, const _Tp x,
                              const bool alternate)
    {

      _Tp sign = _Tp(1);
      if (alternate)
        sign = _Tp(-1);

      const _Tp fnu2 = _Tp(4) * nu * nu;
      const _Tp xx = _Tp(8) * x;

      _Tp term = _Tp(1);
      _Tp sum = _Tp(1);
      for (unsigned int i = 1; i < 20; ++i)
        {
          const _Tp i2 = _Tp(2 * i);
          const _Tp i4 = _Tp(4 * i);
          const _Tp b = i4 - _Tp(1);
          const _Tp numer = sign * (fnu2 - b * b);
          const _Tp denom = i2 * xx;
          term *= numer / denom;
          if (std::abs(term) < std::numeric_limits<_Tp>::epsilon())
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
    template<typename _Tp>
    _Tp
    cyl_bessel_i_asymp(const _Tp n, const _Tp x)
    {

      _Tp sum = cyl_mod_bessel_a_series(n, x, true);

      return std::exp(x) * sum / std::sqrt(2 * _TR1_M_PI * x);
    }


    ///
    /// x >> nu*nu+1
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_i_scaled_asymp(const _Tp nu, const _Tp x)
    {
      const _Tp mu = _Tp(4) * nu * nu;
      const _Tp mum1 = mu - _Tp(1);
      const _Tp mum9 = mu - _Tp(9);
      const _Tp pre = _Tp(1) / std::sqrt(_Tp(2) * _TR1_M_PI * x);
      const _Tp r = mu / x;
      const _Tp I = pre * (_Tp(1) - mum1 / (_Tp(8) * x)
                    + mum1*mum9 / (_Tp(128) * x * x));

      return I;
    }


    //
    //  Debye functions Abramowitz & Stegun, 9.3.9-10
    //
    template<typename _Tp>
    _Tp 
    debye_u1(const std::vector<_Tp> & tpow)
    {
      return (_Tp(3UL) * tpow[1]
            - _Tp(5UL) * tpow[3]) / _Tp(24UL);
    }

    template<typename _Tp>
    _Tp 
    debye_u2(const std::vector<_Tp> & tpow)
    {
      return (_Tp(81UL) * tpow[2]
           - _Tp(462UL) * tpow[4]
           + _Tp(385UL) * tpow[6]) / _Tp(1152UL);
    }

    template<typename _Tp>
    _Tp
    debye_u3(const std::vector<_Tp> & tpow)
    {
      return (_Tp(30375UL) * tpow[3]
           - _Tp(369603UL) * tpow[5]
           + _Tp(765765UL) * tpow[7]
           - _Tp(425425UL) * tpow[9]) / _Tp(414720UL);
    }

    template<typename _Tp>
    _Tp
    debye_u4(const std::vector<_Tp> & tpow)
    {
      return (_Tp(4465125UL) * tpow[4]
           - _Tp(94121676UL) * tpow[6]
          + _Tp(349922430UL) * tpow[8]
          - _Tp(446185740UL) * tpow[10]
          + _Tp(185910725UL) * tpow[12]) / _Tp(39813120UL);
    }

    template<typename _Tp>
    _Tp
    debye_u5(const std::vector<_Tp> & tpow)
    {
      return (_Tp(1519035525UL) * tpow[5]
           - _Tp(49286948607UL) * tpow[7]
          + _Tp(284499769554UL) * tpow[9]
          - _Tp(614135872350UL) * tpow[11]
          + _Tp(566098157625UL) * tpow[13]
          - _Tp(188699385875UL) * tpow[15]) / _Tp(6688604160UL);
    }


#if ! _GLIBCXX_USE_C99_MATH_TR1
    template<typename _Tp>
    _Tp
    xhypot(const _Tp x, const _Tp y)
    {
      const _Tp ax = std::abs(x);
      const _Tp ay = std::abs(y);
      if (ay > ax)
        {
          const _Tp rho = ax / ay;
          return ay * std::sqrt(_Tp(1) + rho * rho);
        }
      else
        {
          const _Tp rho = ay / ax;
          return ax * std::sqrt(_Tp(1) + rho * rho);
        }
    }
#endif


    ///
    ///  @brief  Return the uniform asymptotic approximation for the
    ///          modofied Bessel function @f$ K_\nu(x) @f$.
    ///
    ///  nu -> Inf; uniform in x > 0.  (See Abramowitz & Stegun, 9.7.7)
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_i_scaled_asymp_unif(const _Tp nu, const _Tp x)
    {
      const _Tp z = x / nu;
#if _GLIBCXX_USE_C99_MATH_TR1
      const _Tp hypot = std::tr1::hypot(_Tp(1), z);
#else
      const _Tp hypot = xhypot(_Tp(1), z);
#endif
      const _Tp pre = _Tp(1) / std::sqrt(_Tp(2) * _TR1_M_PI * nu * hypot);
      const _Tp eta = hypot + std::log(z / (_Tp(1) + hypot));
      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp root3_eps = std::exp(std::log(eps) / _Tp(3));
      const _Tp arg = ( z < _Tp(1) / root3_eps
                          ? nu * (-z + eta)
                          : -_Tp(0.5L) * nu / z * (_Tp(1) - _Tp(1) / (_Tp(12) * z * z)) );
      const _Tp ex = std::exp(arg);
      const _Tp t = _Tp(1) / hypot;
      std::vector<_Tp> tpow;
      tpow.push_back(_Tp(1));
      for (int i = 1; i < 16; ++i)
        tpow.push_back(t * tpow.back());
      _Tp sum = _Tp(1);
                + debye_u1(tpow) / (nu)
                + debye_u2(tpow) / (nu*nu)
                + debye_u3(tpow) / (nu*nu*nu)
                + debye_u4(tpow) / (nu*nu*nu*nu)
                + debye_u5(tpow) / (nu*nu*nu*nu*nu);
      const _Tp I = pre * ex * sum;
      return I;
      //On exponent failure: _Tp I = _Tp(0);
      //return I;
    }


    ///
    ///  @brief  Return the uniform asymptotic approximation for the
    ///          modofied Bessel function @f$ K_\nu(x) @f$.
    ///
    ///  nu -> Inf; uniform in x > 0.  (See Abramowitz & Stegun, 9.7.8)
    ///
    template<typename _Tp>
    _Tp
    gsl_sf_bessel_k_scaled_asymp_unif(const _Tp nu, const _Tp x)
    {
      const _Tp z = x / nu;
#if _GLIBCXX_USE_C99_MATH_TR1
      const _Tp hypot = std::tr1::hypot(_Tp(1), z);
#else
      const _Tp hypot = xhypot(_Tp(1), z);
#endif
      const _Tp pre = std::sqrt(_TR1_M_PI / (_Tp(2) * nu * hypot));
      const _Tp eta = hypot + std::log(z / (_Tp(1) + hypot));
      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp root3_eps = std::exp(std::log(eps) / _Tp(3));
      const _Tp arg = ( z < _Tp(1) / root3_eps
                          ? nu * (z - eta)
                          : _Tp(0.5L) * nu / z * (_Tp(1) + _Tp(1) / (_Tp(12) * z * z)) );
      const _Tp ex = std::exp(arg);
      const _Tp t = _Tp(1) / hypot;
      std::vector<_Tp> tpow;
      tpow.push_back(_Tp(1));
      for(int i = 1; i < 16; ++i)
        tpow.push_back(t * tpow.back());
      _Tp fact = _Tp(1);
      _Tp sum = _Tp(1)
                - debye_u1(tpow) / (nu)
                + debye_u2(tpow) / (nu*nu)
                - debye_u3(tpow) / (nu*nu*nu)
                + debye_u4(tpow) / (nu*nu*nu*nu)
                - debye_u5(tpow) / (nu*nu*nu*nu*nu);
      const _Tp K = pre * ex * sum;
      return K;
    }


    ///
    ///  @brief  Evaluate the continued fraction CF1 for @f$ I_{nu+1}/I_nu @f$
    ///          using Gautschi (Euler) equivalent series.
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_i_cont_frac_1_series(const _Tp nu, const _Tp x)
    {
      const int maxk = 20000;
      _Tp tk = _Tp(1);
      _Tp sum = _Tp(1);
      _Tp rhok = _Tp(0);

      for (int k = 1; k < maxk; ++k)
        {
          _Tp ak = (x / (nu + k))
                   * x / (_Tp(4) * (nu + _Tp(k + 1)));
          rhok = -ak * (_Tp(1) + rhok)
                 / (_Tp(1) + ak * (_Tp(1) + rhok));
          tk *= rhok;
          sum += tk;
          if (std::abs(tk / sum) < std::numeric_limits<_Tp>::epsilon())
            break;
        }

      _Tp ratio = x / (_Tp(2) * (nu + _Tp(1))) * sum;

      //if (k == maxk)
      //  throw_logic;

      return ratio;
    }


    /// Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
    /// to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
    ///
    /// This is unstable for small x; x > 2 is a good cutoff.
    /// Also requires |nu| < 1/2.
    template<typename _Tp>
    void
    cyl_bessel_k_scaled_steed_temme_cont_frac_2(const _Tp nu,
                                                  const _Tp x,
                                                  _Tp & K_nu,
                                                  _Tp & K_nup1,
                                                  _Tp & Kp_nu)
    {
      const int max_iter = 10000;

      _Tp bi = _Tp(2) * (_Tp(1) + x);
      _Tp di = _Tp(1) / bi;
      _Tp delhi = di;
      _Tp hi = di;

      _Tp qi   = _Tp(0);
      _Tp qip1 = _Tp(1);

      _Tp ai = -(_Tp(0.25L) - nu * nu);
      _Tp a1 = ai;
      _Tp ci = -ai;
      _Tp Qi = -ai;

      _Tp s = _Tp(1) + Qi * delhi;

      for (int i = 2; i <= max_iter; ++i)
        {
          ai -= _Tp(2 * (i - 1));
          ci  = -ai * ci / i;
          const _Tp tmp = (qi - bi * qip1) / ai;
          qi = qip1;
          qip1 = tmp;
          Qi += ci * qip1;
          bi += _Tp(2);
          di = _Tp(1) / (bi + ai * di);
          delhi = (bi * di - _Tp(1)) * delhi;
          hi += delhi;
          const _Tp dels = Qi * delhi;
          s += dels;
          if (std::abs(dels / s) < std::numeric_limits<_Tp>::epsilon())
            break;
        }

      hi *= -a1;

      K_nu = std::sqrt(_TR1_M_PI/(_Tp(2) * x)) / s;
      K_nup1 = K_nu * (nu + x + _Tp(0.5L) - hi) / x;
      Kp_nu = - K_nup1 + nu / x * K_nu;
      //if(i == max_iter)
      //  GSL_ERROR ("error", GSL_EMAXITER);

      return;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ I_{\nu} \f$.
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_i_scaled(const _Tp nu, const _Tp x)
    {

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_i_scaled.");

      if (std::isnan(x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (x == _Tp(0))
        {
          if (nu == _Tp(0))
            return _Tp(1);
          else
            return _Tp(0);
        }

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp root3_eps = std::exp(std::log(eps) / _Tp(3));

      if (x * x < _Tp(10) * (nu + _Tp(1)))
        {
          const _Tp ex = exp(-x);
          const _Tp b = cyl_bessel_ij_series(nu, x, _Tp(1), 100);
          return b * ex;
        }
      else if (_Tp(0.5L) / (nu * nu + x * x) < root3_eps)
        {
      //  return cyl_bessel_i_asymp(nu, x);
          return cyl_bessel_i_scaled_asymp_unif(nu, x);
        }
      else
        {
          int Nmu = static_cast<int>(nu + _Tp(0.5L));
          _Tp mu = nu - Nmu;  // -1/2 <= mu <= 1/2
          _Tp K_mu, K_mup1, Kp_mu;
          _Tp K_nu, K_nup1, K_num1;
          _Tp I_nu_ratio;

          //  Obtain K_mu, K_mup1
          if (x < _Tp(2))
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
              K_nup1 = _Tp(2) * (mu + _Tp(n + 1)) * K_nu / x + K_num1;
            }

          //  Calculate I_{nu+1}/I_nu
          I_nu_ratio = cyl_bessel_i_cont_frac_1_series(nu, x);

          // solve for I_nu
          _Tp I_nu = _Tp(1) / (x * (K_nup1 + I_nu_ratio * K_nu));
          I_nu = eps * (_Tp(0.5L) * Nmu + _Tp(2)) * std::abs(I_nu);
          return I_nu;
       }
    }


    ///
    ///
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_i(_Tp nu, _Tp x)
    {
      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_i.");

      _Tp I_scaled = cyl_bessel_i_scaled(nu, x);
      return I_scaled * std::exp(x);
    }


    ///
    ///  @brief  Return the asymptotic solution of the modified
    ///          Bessel function \f$ K_{\nu} \f$.
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_k_asymp(const _Tp n, const _Tp x)
    {

      _Tp sum = cyl_mod_bessel_a_series(n, x, false);

      return std::sqrt(_TR1_M_PI_2 / x) * std::exp(-x) * sum;
    }


    ///
    ///  x >> nu*nu + 1
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_k_scaled_asymp(const _Tp nu, const _Tp x)
    {
      const _Tp mu = _Tp(4) * nu * nu;
      const _Tp mum1 = mu - _Tp(1);
      const _Tp mum9 = mu - _Tp(9);
      const _Tp pre = std::sqrt(_TR1_M_PI / (_Tp(2) * x));
      const _Tp r = nu / x;
      const _Tp K = pre * (_Tp(1) + mum1 / (_Tp(8) * x) + mum1 * mum9 / (_Tp(128) * x * x));

      return K;
    }


    ///
    ///  @brief  Return the scaled modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu} \f$.
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_k_scaled(const _Tp nu, const _Tp x)
    {

      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_k_scaled.");

      int Nmu = static_cast<int>(nu + _Tp(0.5L));
      _Tp mu = nu - Nmu;  // -1/2 <= mu <= 1/2
      _Tp K_mu, K_mup1, Kp_mu;

      if (x < _Tp(2))
        cyl_bessel_k_scaled_temme(mu, x, K_mu, K_mup1, Kp_mu);
      else
        cyl_bessel_k_scaled_steed_temme_cont_frac_2(mu, x, K_mu, K_mup1, Kp_mu);

      //  Recurse forward to obtain K_num1, K_nu.

      _Tp K_nu   = K_mu;
      _Tp K_nup1 = K_mup1;

      for (int n = 0; n < Nmu; ++n)
        {
          const _Tp K_num1 = K_nu;
          K_nu   = K_nup1;
          K_nup1 = _Tp(2) * (mu + n + 1) / x * K_nu + K_num1;
        }

      return K_nu;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu} \f$.
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_k(const _Tp nu, const _Tp x)
    {
      if (std::isnan(nu) || std::isnan(x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_k.");

      _Tp K_scaled = cyl_bessel_k_scaled(nu, x);
      return K_scaled * std::exp(-x);
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _TR1_MODIFIED_BESSEL_FUNC_TCC
