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
    ///  @brief  Compute the modified Bessel functions @f$ I_\nu(x) @f$ and
    ///          @f$ K_\nu(x) @f$ and their first derivatives
    ///          @f$ I'_\nu(x) @f$ and @f$ K'_\nu(x) @f$ respectively.
    ///          These four functions are computed together for numerical
    ///          stability.
    ///
    template <typename _Tp>
    void
    bessel_ik(const _Tp nu, const _Tp x,
                _Tp & I_nu, _Tp & K_nu, _Tp & Ip_nu, _Tp & Kp_nu)
    {

      if (x < _Tp(0) || nu < _Tp(0))
        throw_runtime_error(N("Bad arguments in "
                                  "bessel_ik."));

      if (std::isnan(nu) || std::isnan(x))
        {
          I_nu = std::numeric_limits<_Tp>::quiet_NaN();
          K_nu = std::numeric_limits<_Tp>::quiet_NaN();
          Ip_nu = std::numeric_limits<_Tp>::quiet_NaN();
          Kp_nu = std::numeric_limits<_Tp>::quiet_NaN();
          return;
        }

      if (x == _Tp(0))
        {
          if (nu == _Tp(0))
            {
              I_nu = _Tp(1);
              K_nu = std::numeric_limits<_Tp>::infinity();
              Ip_nu = _Tp(0);
              Kp_nu = -std::numeric_limits<_Tp>::infinity();
            }
          else
            {
              I_nu = _Tp(0);
              K_nu = std::numeric_limits<_Tp>::infinity();
              //Ip_nu = ???
              //Kp_nu = ???
            }
          return;
        }

      const _Tp eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp fp_min = _Tp(10) * std::numeric_limits<_Tp>::epsilon();
      const int max_iter = 15000;
      const _Tp x_min = _Tp(2);

      const int nl = static_cast<int>(nu + _Tp(0.5L));

#define _TR1_MU_WART 1
#if _TR1_MU_WART
      //  Bump |mu| away from 0.5.  This is a wart.
      _Tp mu = nu - nl;
      if (mu == -_Tp(0.5L))
        mu += 10 * eps;
      else if (mu == _Tp(0.5L))
        mu -= 10 * eps;
#else
      const _Tp mu = nu - nl;
#endif
      const _Tp mu2 = mu * mu;
      const _Tp xi = _Tp(1) / x;
      const _Tp xi2 = _Tp(2) * xi;
      _Tp h = nu * xi;
      if ( h < fp_min )
        h = fp_min;
      _Tp b = xi2 * nu;
      _Tp d = _Tp(0);
      _Tp c = h;
      int i;
      for ( i = 1; i <= max_iter; ++i )
        {
          b += xi2;
          d = _Tp(1) / (b + d);
          c = b + _Tp(1) / c;
          const _Tp del = c * d;
          h = del * h;
          if (std::abs(del - _Tp(1)) < eps)
            break;
        }
      if (i > max_iter)
        throw_runtime_error(N("Argument x too large in bessel_ik; "
                                  "try asymptotic expansion."));
      _Tp I_nul = fp_min;
      _Tp Ip_nul = h * I_nul;
      _Tp I_nul1 = I_nul;
      _Tp Ip_nu1 = Ip_nul;
      _Tp fact = nu * xi;
      for (int l = nl; l >= 1; --l)
        {
          const _Tp I_nutemp = fact * I_nul + Ip_nul;
          fact -= xi;
          Ip_nul = fact * I_nutemp + I_nul;
          I_nul = I_nutemp;
        }
      _Tp f = Ip_nul / I_nul;
      _Tp K_mu, K_nu1;
      if (x < x_min)
        {
          const _Tp x2 = x / _Tp(2);
          const _Tp pimu = PI * mu;
          const _Tp fact = (std::abs(pimu) < eps
                            ? _Tp(1) : pimu / std::sin(pimu));
          _Tp d = -std::log(x2);
          _Tp e = mu * d;
          const _Tp fact2 = (std::abs(e) < eps
                            ? _Tp(1) : std::sinh(e) / e);
          _Tp gam1, gam2, gampl, gammi;
          gamma_temme(mu, gampl, gammi, gam1, gam2);
          _Tp ff = fact * (gam1 * std::cosh(e)
                             + gam2 * fact2 * d);
          _Tp sum = ff;
          e = std::exp(e);
          _Tp p = e / (_Tp(2) * gampl);
          _Tp q = _Tp(1) / (_Tp(2) * e * gammi);
          _Tp c = _Tp(1);
          d = x2 * x2;
          _Tp sum1 = p;
          int i;
          for (i = 1; i <= max_iter; ++i)
            {
              ff = (i * ff + p + q) / (i * i - mu2);
              c *= d / i;
              p /= i - mu;
              q /= i + mu;
              const _Tp del = c * ff;
              sum += del; 
              const _Tp del1 = c * (p - i * ff);
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
          _Tp b = _Tp(2) * (_Tp(1) + x);
          _Tp d = _Tp(1) / b;
          _Tp delh = d;
          _Tp h = delh;
          _Tp q1 = _Tp(0);
          _Tp q2 = _Tp(1);
          _Tp a1 = 0.25L - mu2;
          _Tp q = c = a1;
          _Tp a = -a1;
          _Tp s = _Tp(1) + q * delh;
          int i;
          for (i = 2; i <= max_iter; ++i)
            {
              a -= 2 * (i - 1);
              c = -a * c / i;
              const _Tp qnew = (q1 - b * q2) / a;
              q1 = q2;
              q2 = qnew;
              q += c * qnew;
              b += _Tp(2);
              d = _Tp(1) / (b + a * d);
              delh = (b * d - _Tp(1)) * delh;
              h += delh;
              const _Tp dels = q * delh;
              s += dels;
              if ( std::abs(dels / s) < eps )
                break;
            }
          if (i > max_iter)
            throw_runtime_error(N("Steed's method failed "
                                      "in bessel_ik."));
          h = a1 * h;
          K_mu = std::sqrt(PI / (_Tp(2) * x)) * std::exp(-x) / s;
          K_nu1 = K_mu * (mu + x + _Tp(0.5L) - h) * xi;
        }

      _Tp K_mup = mu * xi * K_mu - K_nu1;
      _Tp I_mu = xi / (f * K_mu - K_mup);
      I_nu = I_mu * I_nul1 / I_nul;
      Ip_nu = I_mu * Ip_nu1 / I_nul;
      for ( i = 1; i <= nl; ++i )
        {
          const _Tp K_nutemp = (mu + i) * xi2 * K_nu1 + K_mu;
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
    template<typename _Tp>
    _Tp
    cyl_bessel_i(const _Tp nu, const _Tp x)
    {
      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_i.");

      if (std::isnan(nu) || std::isnan(x))
        return std::numeric_limits<_Tp>::quiet_NaN();

      _Tp I_nu, K_nu, Ip_nu, Kp_nu;
      bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);

      return I_nu;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu}(x) \f$.
    ///
    template<typename _Tp>
    _Tp
    cyl_bessel_k(const _Tp nu, const _Tp x)
    {
      if (x < _Tp(0))
        throw std::domain_error("Bad argument in cyl_bessel_k.");

      if (std::isnan(nu) || std::isnan(x))
        return std::numeric_limits<_Tp>::quiet_NaN();

      _Tp I_nu, K_nu, Ip_nu, Kp_nu;
      bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);

      return K_nu;
    }


    ///
    ///  @brief  Compute the spherical modified Bessel functions
    ///          @f$ i_n(x) @f$ and @f$ k_n(x) @f$ and their first
    ///          derivatives @f$ i'_n(x) @f$ and @f$ k'_n(x) @f$
    ///          respectively.
    ///
    template <typename _Tp>
    void
    sph_bessel_ik(const unsigned int n, const _Tp x,
                    _Tp & i_n, _Tp & k_n, _Tp & ip_n, _Tp & kp_n)
    {

      if (x < _Tp(0))
        throw_runtime_error(N("Bad arguments "
                                  "in sph_bessel."));

      const _Tp nu = _Tp(n) + _Tp(0.5L);

      _Tp I_nu, Ip_nu, K_nu, Kp_nu;
      bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);

      const _Tp factor = _Tp(SQRTPIO2) / std::sqrt(x);

      i_n = factor * I_nu;
      k_n = factor * K_nu;
      ip_n = factor * Ip_nu - i_n / (_Tp(2) * x);
      kp_n = factor * Kp_nu - k_n / (_Tp(2) * x);

      return;
    }


    ///
    ///  @brief  Compute the Airy functions
    ///          @f$ Ai(x) @f$ and @f$ Bi(x) @f$ and their first
    ///          derivatives @f$ Ai'(x) @f$ and @f$ Bi(x) @f$
    ///          respectively.
    ///
    template <typename _Tp>
    void
    airy(const _Tp x,
           _Tp & Ai, _Tp & Bi, _Tp & Aip, _Tp & Bip)
    {
      const _Tp absx = std::abs(x);
      const _Tp rootx = std::sqrt(absx);
      const _Tp z = _Tp(2) * absx * rootx / _Tp(3);
      if (x > _Tp(0))
        {
          _Tp I_nu, Ip_nu, K_nu, Kp_nu;

          bessel_jn(_Tp(1)/_Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu);
          Ai = rootx * K_nu / (SQRT3 * PI);
          Bi = rootx * (K_nu / PI + _Tp(2) * I_nu / SQRT3);

          bessel_jn(_Tp(2)/_Tp(3), z, I_nu, K_nu, Ip_nu, Kp_nu);
          Aip = -x * K_nu / (SQRT3 * PI);
          Bip = x * (K_nu / PI + _Tp(2) * I_nu / SQRT3);
        }
      else if (x < _Tp(0))
        {
          _Tp J_nu, Jp_nu, N_nu, Np_nu;

          bessel_jn(_Tp(1)/_Tp(3), z, J_nu, N_nu, Jp_nu, Np_nu);
          Ai = rootx * (J_nu - N_nu / SQRT3) / _Tp(2);
          Bi = -rootx * (N_nu + J_nu / SQRT3) / _Tp(2);

          bessel_jn(_Tp(2)/_Tp(3), z, J_nu, N_nu, Jp_nu, Np_nu);
          Aip = absx * (N_nu / SQRT3 + J_nu) / _Tp(2);
          Bip = absx * (J_nu / SQRT3 - N_nu) / _Tp(2);
        }
      else
        {
          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
          // The number is Ai(0) = 3^{-2/3}/\Gamma(2/3).
          Ai = _Tp(0.35502805388781723926L);
          Bi = Ai * SQRT3;

          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
          // The number is Ai'(0) = -3^{-1/3}/\Gamma(1/3).
          Aip = -_Tp(0.25881940379280679840L);
          Bip = -Aip * SQRT3;
        }

      return;
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _TR1_MODIFIED_BESSEL_FUNC_TCC
