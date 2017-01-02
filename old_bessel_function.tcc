// Special functions -*- C++ -*-

// Copyright (C) 2006-2017
// Free Software Foundation, Inc.
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

#include "specfun_util.h"

namespace std
{
_GLIBCXX_BEGIN_NAMESPACE(tr1)

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
    __cyl_neumann_temme(const _Tp __nu, const _Tp __x,
                        _Tp & __Ynu, _Tp & __Ynup1)
    {
      const _Tp __x2 = __x / _Tp(2);
      const _Tp __ln_x2 = std::log(__x2);
      const _Tp __x2_nu = std::exp(__nu * __ln_x2);
      const _Tp __pi_nu   = _TR1_M_PI * __nu;
      const _Tp __alpha   = __pi_nu / _Tp(2);
      const _Tp __sigma   = -__nu * __ln_x2;
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __sinrat  = (std::abs(__pi_nu) < __eps
                            ? _Tp(1) : __pi_nu / std::sin(__pi_nu));
      const _Tp __sinhrat = (std::abs(__sigma) < __eps
                            ? _Tp(1) : std::sinh(__sigma) / __sigma);
      const _Tp __sinhalf = (std::abs(__alpha) < __eps
                            ? _Tp(1) : std::sin(__alpha) / __alpha);
      const _Tp __sin_sqr = __nu * _TR1_M_PI * _TR1_M_PI
                          * __sinhalf * __sinhalf / _Tp(2);

      _Tp __g_1pnu, __g_1mnu, __g1, __g2;
      __gamma_temme(__nu, __g_1pnu, __g_1mnu, __g1, __g2);

      _Tp __fk = (_Tp(2) / _TR1_M_PI) * __sinrat
               * (std::cosh(__sigma) * __g1 - __sinhrat * __ln_x2 * __g2);
      _Tp __pk = _Tp(1) / _TR1_M_PI / __x2_nu * __g_1pnu;
      _Tp __qk = _Tp(1) / _TR1_M_PI * __x2_nu * __g_1mnu;
      _Tp __hk = __pk;
      _Tp __ck = _Tp(1);

      _Tp __sum0 = __fk + __sin_sqr * __qk;
      _Tp __sum1 = __pk;

      const unsigned int __max_iter = 15000;
      unsigned int __k = 0;
      while (__k < __max_iter)
        {
          ++__k;
          __fk  = (__k * __fk + __pk + __qk) / (__k * __k - __nu * __nu);
          __ck *= -__x2 * __x2 / __k;
          __pk /= (__k - __nu);
          __qk /= (__k + __nu);
          const _Tp __gk  = __fk + __sin_sqr * __qk;
          __hk  = -__k * __gk + __pk; 
          const _Tp __del0 = __ck * __gk;
          const _Tp __del1 = __ck * __hk;
          __sum0 += __del0;
          __sum1 += __del1;
          if (std::abs(__del0) < 0.5L * (_Tp(1) + std::abs(__sum0)) * __eps)
            break;
        }

      __Ynu   = -__sum0;
      __Ynup1 = -__sum1 * _Tp(2) / __x;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Bessel function
    ///         of order nu: \f$J_{\nu}\f$.
    ///
    template <typename _Tp>
    _Tp
    __cyl_bessel_j_asymp(const _Tp __nu, const _Tp __x)
    {
      const _Tp __coef = std::sqrt(_Tp(2)/(_Tp(_TR1_M_PI) * __x));
      const _Tp __mu   = _Tp(4) * __nu * __nu;
      const _Tp __mum1 = __mu - _Tp(1);
      const _Tp __mum9 = __mu - _Tp(9);
      const _Tp __mum25 = __mu - _Tp(25);
      const _Tp __P = _Tp(1) - __mum1 * __mum9 / (_Tp(128) * __x * __x);
      const _Tp __Q = __mum1/(_Tp(8) * __x)
                    * (_Tp(1) - __mum9 * __mum25 / (_Tp(384) * __x * __x));

      const _Tp __arg = __x - (__nu + _Tp(0.5L)) * _Tp(_TR1_M_PI_2);
      const _Tp __c = std::cos(__arg);
      const _Tp __s = std::sin(__arg);

      const _Tp __Jnu = __coef * (__c * __P - __s * __Q);

      return __Jnu;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Bessel function
    ///         amplitude and phase functions of order nu: \f$ M_{\nu} \f$
    ///         and \f$ \theta_{\nu} \f$
    ///
    template <typename _Tp>
    void
    __cyl_bessel_asymp_amp_phase(const _Tp __nu, const _Tp __x,
                                 _Tp & __amp, _Tp & __phase)
    {
      const _Tp __r  = _Tp(2) * __nu / __x;
      const _Tp __rr = __r * __r;
      const _Tp __xx = __x * __x;

      const _Tp __amp_term1 = (__rr - _Tp(1) / __xx) / _Tp(8);
      const _Tp __amp_term2 = __amp_term1
                            * (__rr - _Tp(9) / __xx) * _Tp(3) / _Tp(16);
      __amp = (_Tp(2) / _Tp(_TR1_M_PI))
            * (_Tp(1) + __amp_term1 + __amp_term2);
      __amp = std::sqrt(__amp) / std::sqrt(__x);

      const _Tp __ph_term1 = __x * (__rr - _Tp(1) / __xx) / _Tp(8);
      const _Tp __ph_term2 = __ph_term1
                           * (__rr - _Tp(25) / __xx) / _Tp(48);
      __phase = -_Tp(_TR1_M_PI) / _Tp(4) + __ph_term1 + __ph_term2;

      return;
    }


    ///
    ///  @brief This routine computes the asyptotic cylindrical Neumann
    ///         function of order \f$ \nu \f$: \f$ N_{\nu} \f$.
    ///
    template <typename _Tp>
    _Tp
    __cyl_neumann_asymp(const _Tp __nu, const _Tp __x)
    {
      _Tp __amp, __phase;
      __cyl_bessel_asymp_amp_phase(__nu, __x, __amp, __phase);
      const _Tp __sin_term = std::sin(__x) * std::cos(__phase)
                           + std::cos(__x) * std::sin(__phase);
      const _Tp __N_nu = __amp * __sin_term;

      return __N_nu;
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
    __cyl_bessel_ij_series(const _Tp __nu, const _Tp __x, const _Tp __sgn,
                           const unsigned int __max_iter)
    {

      const _Tp __x2 = __x / _Tp(2);
      _Tp __fact = __nu * std::log(__x2);
#if _GLIBCXX_USE_C99_MATH_TR1
      __fact -= std::tr1::lgamma(__nu + _Tp(1));
#else
      __fact -= __log_gamma(__nu + _Tp(1));
#endif
      __fact = std::exp(__fact);
      const _Tp __xx4 = __sgn * __x2 * __x2;
      _Tp __Jn = _Tp(1);
      _Tp __term = _Tp(1);

      for (unsigned int __i = 1; __i < __max_iter; ++__i)
        {
          __term *= __xx4 / (_Tp(__i) * (__nu + _Tp(__i)));
          __Jn += __term;
          if (std::abs(__term / __Jn) < std::numeric_limits<_Tp>::epsilon())
            break;
        }

      return __fact * __Jn;
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
    __bessel_j_cont_frac_1(const _Tp __nu, const _Tp __x,
                           _Tp & __ratio, _Tp & __sgn)
    {
      const _Tp __big = std::numeric_limits<_Tp>::max();
      const int __maxiter = 10000;
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      int __n = 1;
      _Tp __Anm2 = _Tp(1);
      _Tp __Bnm2 = _Tp(0);
      _Tp __Anm1 = _Tp(0);
      _Tp __Bnm1 = _Tp(1);
      _Tp __a1 = __x / (_Tp(2) * (__nu + _Tp(1)));
      _Tp __An = __Anm1 + __a1 * __Anm2;
      _Tp __Bn = __Bnm1 + __a1 * __Bnm2;
      _Tp __an;
      _Tp __dn = __a1;
      __ratio = __An / __Bn;
      __sgn  = _Tp(1);

      while (__n < __maxiter)
        {
          ++__n;
          __Anm2 = __Anm1;
          __Bnm2 = __Bnm1;
          __Anm1 = __An;
          __Bnm1 = __Bn;
          __an = -__x * __x / (_Tp(4)
               * (__nu + _Tp(__n - 1)) * (__nu + _Tp(__n)));
          __An = __Anm1 + __an * __Anm2;
          __Bn = __Bnm1 + __an * __Bnm2;

          if (std::abs(__An) > __big || std::abs(__Bn) > __big)
            {
              __An /= __big;
              __Bn /= __big;
              __Anm1 /= __big;
              __Bnm1 /= __big;
              __Anm2 /= __big;
              __Bnm2 /= __big;
            }

          const _Tp __old_ratio = __ratio;
          __ratio = __An / __Bn;
          const _Tp __del = __old_ratio / __ratio;

          __dn = _Tp(1) / (_Tp(2) * (__nu + __n) / __x - __dn);
          if (__dn < _Tp(0))
            __sgn = -__sgn;

          if (std::abs(__del - _Tp(1)) < _Tp(2) * __eps)
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
    __bessel_jn_steed_cont_frac_2(const _Tp __nu, const _Tp __x,
                                  _Tp & __p, _Tp & __q)
    {
      const int __max_iter = 10000;
      const _Tp _small = std::sqrt(std::numeric_limits<_Tp>::min());
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();

      _Tp __x_inv = _Tp(1) / __x;
      _Tp __a = _Tp(0.25L) - __nu * __nu;
      __p = -__x_inv / _Tp(2);
      __q = _Tp(1);
      _Tp __br = _Tp(2) * __x;
      _Tp __bi = _Tp(2);
      _Tp __fact = __a * __x_inv / (__p * __p + __q * __q);
      _Tp __cr = __br + __q * __fact;
      _Tp __ci = __bi + __p * __fact;
      _Tp __den = __br * __br + __bi * __bi;
      _Tp __dr = __br / __den;
      _Tp __di = -__bi / __den;
      _Tp __dlr = __cr * __dr - __ci * __di;
      _Tp __dli = __cr * __di + __ci * __dr;
      _Tp __temp = __p * __dlr - __q * __dli;
      __q = __p * __dli + __q * __dlr;
      __p = __temp;
      for (int __i = 2; __i <= __max_iter; ++__i)
        {
          __a  += _Tp(2 * (__i - 1));
          __bi += _Tp(2);
          __dr = __a * __dr + __br;
          __di = __a * __di + __bi;
          if (std::abs(__dr) + std::abs(__di) < _small)
            __dr = _small;
          __fact = __a / (__cr * __cr + __ci * __ci);
          __cr = __br + __cr * __fact;
          __ci = __bi - __ci * __fact;
          if (std::abs(__cr) + std::abs(__ci) < _small)
            __cr = _small;
          __den = __dr * __dr + __di * __di;
          __dr /= __den;
          __di /= -__den;
          __dlr = __cr * __dr - __ci * __di;
          __dli = __cr * __di + __ci * __dr;
          __temp = __p * __dlr - __q * __dli;
          __q = __p * __dli + __q * __dlr;
          __p = __temp;
          if (std::abs(__dlr - _Tp(1)) + std::abs(__dli) < __eps)
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
    __cyl_bessel_a_series(const _Tp __nu, const _Tp __x)
    {

      const _Tp __fnu2 = _Tp(4) * __nu * __nu;
      const _Tp __xx = _Tp(64) * __x * __x;

      _Tp __eterm = _Tp(1);
      _Tp __esum = _Tp(1);
      _Tp __oterm = (__fnu2 - _Tp(1)) / (_Tp(8) * __x);
      _Tp __osum = __oterm;
      for (unsigned int __i = 1; __i < 10; ++__i)
        {
          const _Tp __i2 = _Tp(2 * __i);
          const _Tp __i4 = _Tp(4 * __i);
          const _Tp __a = __i4 - _Tp(3);
          const _Tp __b = __i4 - _Tp(1);
          const _Tp __c = __i4 + _Tp(1);
          const _Tp __enumer = -(__fnu2 - __a * __a) * (__fnu2 - __b * __b);
          const _Tp __edenom = (__i2 - _Tp(1)) * (__i2) * __xx;
          __eterm *= __enumer / __edenom;
          const _Tp __onumer = -(__fnu2 - __b * __b) * (__fnu2 - __c * __c);
          const _Tp __odenom = (__i2) * (__i2 + _Tp(1)) * __xx;
          __oterm *= __onumer / __odenom;
          if (std::abs(__eterm) < std::numeric_limits<_Tp>::epsilon()
           && std::abs(__oterm) < std::numeric_limits<_Tp>::epsilon())
            {
              break;
            }
          __esum += __eterm;
          __osum += __oterm;
        }

      return std::make_pair(__esum, __osum);
    }


    ///
    ///  @brief This routine returns the cylindrical Bessel function
    ///         of order \f$ \nu \f$: \f$ J_{\nu} \f$
    ///         by asymptotic expansion.  Olver p ~ 130.
    ///
    template <typename _Tp>
    _Tp
    __cyl_bessel_j_olver_asymp(const _Tp __nu, const _Tp __x)
    {

      std::pair<_Tp,_Tp> __sum = __cyl_bessel_a_series(__nu, __x);

      const _Tp __arg = __x - (__nu + _Tp(0.5L)) * _Tp(_TR1_M_PI_2);
      _Tp __Jn = std::sqrt(_Tp(2) / (_Tp(_TR1_M_PI) * __x))
               * (std::cos(__arg) * __sum.first
                - std::sin(__arg) * __sum.second);

      return __Jn;
    }


    ///
    ///  Return the cylindrical Bessel function of order \f$ \nu \f$:
    ///  \f$ J_{\nu} \f$.
    ///
    template <typename _Tp>
    _Tp
    __cyl_bessel_j(const _Tp __nu, const _Tp __x)
    {

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument in __cyl_bessel_j."));

      if (std::isnan(__nu) || std::isnan(__x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        {
          if (__nu == _Tp(0))
            return _Tp(1);
          else
            return _Tp(0);
        }

      if (__x * __x < _Tp(10) * (__nu + _Tp(1)))
        return __cyl_bessel_ij_series(__nu, __x, _Tp(-1), 100);
      else if (__nu > _Tp(50))
        return __cyl_bessel_j_olver_asymp(__nu, __x);
      else if (__x > _Tp(1000))
        return __cyl_bessel_j_asymp(__nu, __x);
      else
        {
          // -1/2 <= mu <= 1/2
          int __Nmu = static_cast<int>(__nu + 0.5L);
          _Tp __mu = __nu - _Tp(__Nmu);

          //  Get J ratio.
          _Tp __Jnup1_Jnu, __sgn_Jnu;
          __bessel_j_cont_frac_1(__nu, __x, __Jnup1_Jnu, __sgn_Jnu);

          if (__x < _Tp(2))
            {
              // Determine Y_mu, Y_mup1 directly and recurse forward to nu.
              // Then use the CF1 information to solve for J_nu and J_nup1.
              _Tp __Y_mu, __Y_mup1;
              __cyl_neumann_temme(__mu, __x, __Y_mu, __Y_mup1);

              _Tp __Ynm1 = __Y_mu;
              _Tp __Yn   = __Y_mup1;
              _Tp __Ynp1 = _Tp(0);
              for (int __n = 1; __n < __Nmu; ++__n)
                {
                  __Ynp1 = _Tp(2) * (__mu + __n) * __Yn / __x - __Ynm1;
                  __Ynm1 = __Yn;
                  __Yn = __Ynp1;
                }

              const _Tp __result = _Tp(2) / (_TR1_M_PI * __x)
                                 / (__Jnup1_Jnu * __Yn - __Ynp1);
              return __result;
            }
          else
            {
              //  Recurse backward from nu to mu, determining the J ratio
              //  at mu. Use this together with a Steed method CF2 to
              //  determine the actual J_mu, and thus obtain the normalization.
              _Tp __Jmu;
              _Tp __Jmup1_Jmu;
              _Tp __sgn_Jmu;
              _Tp __Jmuprime_Jmu;
              _Tp __P, __Q;
              __bessel_jn_steed_cont_frac_2(__mu, __x, __P, __Q);
              _Tp __gamma;
              const _Tp __sqrt_min = std::sqrt(std::numeric_limits<_Tp>::min());
 
              _Tp __Jnp1 = __sgn_Jnu * __sqrt_min * __Jnup1_Jnu;
              _Tp __Jn   = __sgn_Jnu * __sqrt_min;
              _Tp __Jnm1;
              for (int __n = __Nmu; __n > 0; --__n)
                {
                  __Jnm1 = _Tp(2) * (__mu + __n) * __Jn / __x - __Jnp1;
                  __Jnp1 = __Jn;
                  __Jn = __Jnm1;
                }
              __Jmup1_Jmu = __Jnp1 / __Jn;
              __sgn_Jmu = (__Jn < _Tp(0) ? _Tp(-1) : _Tp(1));
              __Jmuprime_Jmu = __mu / __x - __Jmup1_Jmu;

              __gamma = (__P - __Jmuprime_Jmu) / __Q;
              __Jmu = __sgn_Jmu
                  * std::sqrt(2 / (_TR1_M_PI * __x) / (__Q + __gamma * (__P - __Jmuprime_Jmu)));

              const _Tp __result = __Jmu * (__sgn_Jnu * __sqrt_min) / __Jn;

              return __result;
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
    __cyl_neumann_n_olver_asymp(const _Tp __nu, const _Tp __x)
    {

      std::pair<_Tp,_Tp> __sum = __cyl_bessel_a_series(__nu, __x);

      const _Tp __arg = __x - (__nu + _Tp(0.5L)) * _Tp(_TR1_M_PI_2);
      const _Tp __N_nu = std::sqrt(_Tp(2) / (_Tp(_TR1_M_PI) * __x))
                       * (std::sin(__arg) * __sum.first
                        + std::cos(__arg) * __sum.second);

      return __N_nu;
    }


    ///
    ///  Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
    ///
    template <typename _Tp>
    void
    __bessel_jn_mu_restricted(const _Tp __mu, const _Tp __x,
                              _Tp & __Jmu, _Tp & __Jmup1,
                              _Tp & __Ymu, _Tp & __Ymup1)
    {

      if (__x < _Tp(0) || std::abs(__mu) > _Tp(0.5L))
        {
          __Jmu = _Tp(0);
          __Jmup1 = _Tp(0);
          __Ymu = _Tp(0);
          __Ymup1 = _Tp(0);
        }
      else if (__x == _Tp(0))
        {
          if (__mu == _Tp(0))
            __Jmu = _Tp(1);
          else
            __Jmu = _Tp(0);
          __Jmup1 = _Tp(0);
          __Ymu = _Tp(0);
          __Ymup1 = _Tp(0);
        }
      else {

        if (__x < _Tp(2))
          {
            //  Use Taylor series for J and the Temme series for Y.
            //  The Taylor series for J requires nu > 0, so we shift
            //  up one and use the recursion relation to get Jmu, in
            //  case mu < 0.
            _Tp __Jmup2;
            __Jmup1 = __cyl_bessel_ij_series(__mu + _Tp(1), __x, _Tp(-1), 100);
            __Jmup2 = __cyl_bessel_ij_series(__mu + _Tp(2), __x, _Tp(-1), 100);
            __Jmu  = _Tp(2) * (__mu + _Tp(1)) * __Jmup1 / __x - __Jmup2;
            __cyl_neumann_temme(__mu, __x, __Ymu, __Ymup1);
          }
        else if (__x < _Tp(1000))
          {
            _Tp __J_ratio, __J_sgn;
            __bessel_j_cont_frac_1(__mu, __x, __J_ratio, __J_sgn);
            _Tp __P, __Q;
            __bessel_jn_steed_cont_frac_2(__mu, __x, __P, __Q);
            _Tp __Jprime_J_ratio = __mu / __x - __J_ratio;
            _Tp __gamma = (__P - __Jprime_J_ratio) / __Q;
            __Jmu = __J_sgn * std::sqrt(_Tp(2) / (_TR1_M_PI * __x)
                  / (__Q + __gamma * (__P - __Jprime_J_ratio)));
            __Jmup1 = __J_ratio * __Jmu;
            __Ymu = __gamma * __Jmu;
            __Ymup1 = __Ymu * (__mu / __x - __P - __Q / __gamma);
            return;
          }
        else
          {
            //  Use asymptotics for large argument.
            __Jmu = __cyl_bessel_j_asymp(__mu, __x );
            __Jmup1 = __cyl_bessel_j_asymp(__mu + _Tp(1), __x);
            __Ymu = __cyl_neumann_asymp(__mu, __x);
            __Ymup1 = __cyl_neumann_asymp(__mu + _Tp(1), __x);
            return;
          }
      }
    }


    ///
    ///  Return the cylindrical Neumann function of order nu: \f$ N_{\nu}(x) \f$.
    ///
    template <typename _Tp>
    _Tp
    __cyl_neumann_n(const _Tp __nu, const _Tp __x)
    {

      if (std::isnan(__nu) || std::isnan(__x))
        return - std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        return std::numeric_limits<_Tp>::infinity();

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument in __cyl_neumann_n."));

      if (__nu > _Tp(50))
        {
          return __cyl_neumann_n_olver_asymp(__nu, __x);
        }
      else
        {
          int __Nmu = static_cast<int>(__nu + _Tp(0.5L));
          // -1/2 <= mu <= 1/2
          _Tp __mu = __nu - __Nmu;

          _Tp __Y_mu, __Y_mup1;

          if (__x < _Tp(2))
            {
              //  Determine Ymu, Ymup1 directly.
              __cyl_neumann_temme(__mu, __x, __Y_mu, __Y_mup1);
            }
          else
            {
              //  Determine Ymu, Ymup1 and Jmu, Jmup1.
              _Tp __J_mu, __J_mup1;
              __bessel_jn_mu_restricted(__mu, __x, __J_mu, __J_mup1, __Y_mu, __Y_mup1);
            }

          //  Forward recursion to get Ynu, Ynup1.
          _Tp __Ynm1 = __Y_mu;
          _Tp __Yn   = __Y_mup1;
          for (int __n = 1; __n <= __Nmu; ++__n)
            {
              _Tp __Ynp1 = _Tp(2) * (__mu + __n) * __Yn / __x - __Ynm1;
              __Ynm1 = __Yn;
              __Yn   = __Ynp1;
            }

          _Tp __result  = __Ynm1; // Y_nu

          return __result;
        }
    }


    ///
    ///  Return the spherical Bessel function.
    ///
    template <typename _Tp>
    _Tp
    __sph_bessel(const unsigned int __n, const _Tp __x)
    {

      if (std::isnan(__x))
        return std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        return std::numeric_limits<_Tp>::infinity();  //  TODO: This is wrong?

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument in __sph_bessel."));

      return std::sqrt(_Tp(_TR1_M_PI_2) / __x)
           * __cyl_bessel_j(_Tp(__n) + _Tp(0.5L), __x);
    }


    ///
    ///  Return the spherical Neumann function.
    ///
    template <typename _Tp>
    _Tp
    __sph_neumann(const unsigned int __n, const _Tp __x)
    {

      if (std::isnan(__x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        return std::numeric_limits<_Tp>::infinity();  //  TODO: This is wrong?

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument in __sph_neumann."));

      return std::sqrt(_Tp(_TR1_M_PI_2) / __x)
           * __cyl_neumann_n(_Tp(__n) + _Tp(0.5L), __x);
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE
}

#endif // _TR1_BESSEL_FUNCTION_TCC
