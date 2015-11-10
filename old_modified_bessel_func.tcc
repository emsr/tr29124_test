// Special functions -*- C++ -*-

// Copyright (C) 2006-2007
// Free Software Foundation, Inc.
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
    ///
    ///
    template <typename _Tp>
    void
    __cyl_bessel_k_scaled_temme(const _Tp __nu, const _Tp __x,
                                _Tp & __K_nu, _Tp & __K_nup1, _Tp & __Kp_nu)
    {
      const _Tp __x2    = __x / _Tp(2);
      const _Tp __ln_x2 = std::log(__x2);
      const _Tp __x2_nu = std::exp(__nu * __ln_x2);
      const _Tp __pi_nu   = _TR1_M_PI * __nu;
      const _Tp __sigma   = -__nu * __ln_x2;
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __sinrat  = (std::abs(__pi_nu) < __eps ? _Tp(1) : __pi_nu / std::sin(__pi_nu));
      const _Tp __sinhrat = (std::abs(__sigma) < __eps ? _Tp(1) : ::sinh(__sigma) / __sigma);
      const _Tp __ex = std::exp(__x);

      _Tp __g_1pnu, __g_1mnu, __g1, __g2;
      __gamma_temme(__nu, __g_1pnu, __g_1mnu, __g1, __g2);

      _Tp __fk = __sinrat * (::cosh(__sigma) * __g1
                           - __sinhrat * __ln_x2 * __g2);
      _Tp __pk = 0.5 / __x2_nu * __g_1pnu;
      _Tp __qk = 0.5 * __x2_nu * __g_1mnu;
      _Tp __hk = __pk;
      _Tp __ck = _Tp(1);
      _Tp __sum0 = __fk;
      _Tp __sum1 = __hk;

      const unsigned int __max_iter = 15000;
      unsigned int __k = 0;
      while(__k < __max_iter)
        {
          ++__k;
          __fk  = (__k * __fk + __pk + __qk) / (__k * __k - __nu * __nu);
          __ck *= __x2 * __x2 / __k;
          __pk /= (__k - __nu);
          __qk /= (__k + __nu);
          __hk  = -__k * __fk + __pk;
          const _Tp __del0 = __ck * __fk;
          const _Tp __del1 = __ck * __hk;
          __sum0 += __del0;
          __sum1 += __del1;
          if (std::abs(__del0) < 0.5L * std::abs(__sum0) * __eps)
            break;
        }
  
      __K_nu   = __sum0 * __ex;
      __K_nup1 = __sum1 * (_Tp(2) / __x) * __ex;
      __Kp_nu  = - __K_nup1 + (__nu / __x) * __K_nu;

      return;
    }


    ///
    ///  @brief  Return the A series used in the asymptotic expansions
    ///          of the Bessel functions.
    ///
    template<typename _Tp>
    _Tp
    __cyl_mod_bessel_a_series(const _Tp __nu, const _Tp __x,
                              const bool __alternate)
    {

      _Tp __sign = _Tp(1);
      if (__alternate)
        __sign = _Tp(-1);

      const _Tp __fnu2 = _Tp(4) * __nu * __nu;
      const _Tp __xx = _Tp(8) * __x;

      _Tp __term = _Tp(1);
      _Tp __sum = _Tp(1);
      for (unsigned int __i = 1; __i < 20; ++__i)
        {
          const _Tp __i2 = _Tp(2 * __i);
          const _Tp __i4 = _Tp(4 * __i);
          const _Tp __b = __i4 - _Tp(1);
          const _Tp __numer = __sign * (__fnu2 - __b * __b);
          const _Tp __denom = __i2 * __xx;
          __term *= __numer / __denom;
          if (std::abs(__term) < std::numeric_limits<_Tp>::epsilon())
            {
              break;
            }
          __sum += __term;
        }

      return __sum;
    }


    ///
    ///  @brief  Return the asymptotic solution of the Bessel function.
    ///
    template<typename _Tp>
    inline _Tp
    __cyl_bessel_i_asymp(const _Tp __n, const _Tp __x)
    {

      _Tp __sum = __cyl_mod_bessel_a_series(__n, __x, true);

      return std::exp(__x) * __sum / std::sqrt(2 * _TR1_M_PI * __x);
    }


    ///
    /// x >> nu*nu+1
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_i_scaled_asymp(const _Tp __nu, const _Tp __x)
    {
      const _Tp __mu = _Tp(4) * __nu * __nu;
      const _Tp __mum1 = __mu - _Tp(1);
      const _Tp __mum9 = __mu - _Tp(9);
      const _Tp __pre = _Tp(1) / std::sqrt(_Tp(2) * _TR1_M_PI * __x);
      const _Tp __r = __mu / __x;
      const _Tp __I = __pre * (_Tp(1) - __mum1 / (_Tp(8) * __x)
                    + __mum1*__mum9 / (_Tp(128) * __x * __x));

      return __I;
    }


    //
    //  Debye functions Abramowitz & Stegun, 9.3.9-10
    //
    template<typename _Tp>
    inline _Tp 
    __debye_u1(const std::vector<_Tp> & __tpow)
    {
      return (_Tp(3UL) * __tpow[1]
            - _Tp(5UL) * __tpow[3]) / _Tp(24UL);
    }

    template<typename _Tp>
    inline _Tp 
    __debye_u2(const std::vector<_Tp> & __tpow)
    {
      return (_Tp(81UL) * __tpow[2]
           - _Tp(462UL) * __tpow[4]
           + _Tp(385UL) * __tpow[6]) / _Tp(1152UL);
    }

    template<typename _Tp>
    inline _Tp
    __debye_u3(const std::vector<_Tp> & __tpow)
    {
      return (_Tp(30375UL) * __tpow[3]
           - _Tp(369603UL) * __tpow[5]
           + _Tp(765765UL) * __tpow[7]
           - _Tp(425425UL) * __tpow[9]) / _Tp(414720UL);
    }

    template<typename _Tp>
    inline _Tp
    __debye_u4(const std::vector<_Tp> & __tpow)
    {
      return (_Tp(4465125UL) * __tpow[4]
           - _Tp(94121676UL) * __tpow[6]
          + _Tp(349922430UL) * __tpow[8]
          - _Tp(446185740UL) * __tpow[10]
          + _Tp(185910725UL) * __tpow[12]) / _Tp(39813120UL);
    }

    template<typename _Tp>
    inline _Tp
    __debye_u5(const std::vector<_Tp> & __tpow)
    {
      return (_Tp(1519035525UL) * __tpow[5]
           - _Tp(49286948607UL) * __tpow[7]
          + _Tp(284499769554UL) * __tpow[9]
          - _Tp(614135872350UL) * __tpow[11]
          + _Tp(566098157625UL) * __tpow[13]
          - _Tp(188699385875UL) * __tpow[15]) / _Tp(6688604160UL);
    }


#if ! _GLIBCXX_USE_C99_MATH_TR1
    template<typename _Tp>
    _Tp
    __xhypot(const _Tp __x, const _Tp __y)
    {
      const _Tp __ax = std::abs(__x);
      const _Tp __ay = std::abs(__y);
      if (__ay > __ax)
        {
          const _Tp __rho = __ax / __ay;
          return __ay * std::sqrt(_Tp(1) + __rho * __rho);
        }
      else
        {
          const _Tp __rho = __ay / __ax;
          return __ax * std::sqrt(_Tp(1) + __rho * __rho);
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
    __cyl_bessel_i_scaled_asymp_unif(const _Tp __nu, const _Tp __x)
    {
      const _Tp __z = __x / __nu;
#if _GLIBCXX_USE_C99_MATH_TR1
      const _Tp __hypot = std::tr1::hypot(_Tp(1), __z);
#else
      const _Tp __hypot = __xhypot(_Tp(1), __z);
#endif
      const _Tp __pre = _Tp(1) / std::sqrt(_Tp(2) * _TR1_M_PI * __nu * __hypot);
      const _Tp __eta = __hypot + std::log(__z / (_Tp(1) + __hypot));
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __root3_eps = std::exp(std::log(__eps) / _Tp(3));
      const _Tp __arg = ( __z < _Tp(1) / __root3_eps
                          ? __nu * (-__z + __eta)
                          : -_Tp(0.5L) * __nu / __z * (_Tp(1) - _Tp(1) / (_Tp(12) * __z * __z)) );
      const _Tp __ex = std::exp(__arg);
      const _Tp __t = _Tp(1) / __hypot;
      std::vector<_Tp> __tpow;
      __tpow.push_back(_Tp(1));
      for (int __i = 1; __i < 16; ++__i)
        __tpow.push_back(__t * __tpow.back());
      _Tp __sum = _Tp(1);
                + __debye_u1(__tpow) / (__nu)
                + __debye_u2(__tpow) / (__nu*__nu)
                + __debye_u3(__tpow) / (__nu*__nu*__nu)
                + __debye_u4(__tpow) / (__nu*__nu*__nu*__nu)
                + __debye_u5(__tpow) / (__nu*__nu*__nu*__nu*__nu);
      const _Tp __I = __pre * __ex * __sum;
      return __I;
      //On exponent failure: _Tp __I = _Tp(0);
      //return __I;
    }


    ///
    ///  @brief  Return the uniform asymptotic approximation for the
    ///          modofied Bessel function @f$ K_\nu(x) @f$.
    ///
    ///  nu -> Inf; uniform in x > 0.  (See Abramowitz & Stegun, 9.7.8)
    ///
    template<typename _Tp>
    _Tp
    gsl_sf_bessel_k_scaled_asymp_unif(const _Tp __nu, const _Tp __x)
    {
      const _Tp __z = __x / __nu;
#if _GLIBCXX_USE_C99_MATH_TR1
      const _Tp __hypot = std::tr1::hypot(_Tp(1), __z);
#else
      const _Tp __hypot = __xhypot(_Tp(1), __z);
#endif
      const _Tp __pre = std::sqrt(_TR1_M_PI / (_Tp(2) * __nu * __hypot));
      const _Tp __eta = __hypot + std::log(__z / (_Tp(1) + __hypot));
      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __root3_eps = std::exp(std::log(__eps) / _Tp(3));
      const _Tp __arg = ( __z < _Tp(1) / __root3_eps
                          ? __nu * (__z - __eta)
                          : _Tp(0.5L) * __nu / __z * (_Tp(1) + _Tp(1) / (_Tp(12) * __z * __z)) );
      const _Tp __ex = std::exp(__arg);
      const _Tp __t = _Tp(1) / __hypot;
      std::vector<_Tp> __tpow;
      __tpow.push_back(_Tp(1));
      for(int __i = 1; __i < 16; ++__i)
        __tpow.push_back(__t * __tpow.back());
      _Tp __fact = _Tp(1);
      _Tp __sum = _Tp(1)
                - __debye_u1(__tpow) / (__nu)
                + __debye_u2(__tpow) / (__nu*__nu)
                - __debye_u3(__tpow) / (__nu*__nu*__nu)
                + __debye_u4(__tpow) / (__nu*__nu*__nu*__nu)
                - __debye_u5(__tpow) / (__nu*__nu*__nu*__nu*__nu);
      const _Tp __K = __pre * __ex * __sum;
      return __K;
    }


    ///
    ///  @brief  Evaluate the continued fraction CF1 for @f$ I_{nu+1}/I_nu @f$
    ///          using Gautschi (Euler) equivalent series.
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_i_cont_frac_1_series(const _Tp __nu, const _Tp __x)
    {
      const int __maxk = 20000;
      _Tp __tk = _Tp(1);
      _Tp __sum = _Tp(1);
      _Tp __rhok = _Tp(0);

      for (int __k = 1; __k < __maxk; ++__k)
        {
          _Tp __ak = (__x / (__nu + __k))
                   * __x / (_Tp(4) * (__nu + _Tp(__k + 1)));
          __rhok = -__ak * (_Tp(1) + __rhok)
                 / (_Tp(1) + __ak * (_Tp(1) + __rhok));
          __tk *= __rhok;
          __sum += __tk;
          if (std::abs(__tk / __sum) < std::numeric_limits<_Tp>::epsilon())
            break;
        }

      _Tp __ratio = __x / (_Tp(2) * (__nu + _Tp(1))) * __sum;

      //if (__k == __maxk)
      //  __throw_logic;

      return __ratio;
    }


    /// Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
    /// to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
    ///
    /// This is unstable for small x; x > 2 is a good cutoff.
    /// Also requires |nu| < 1/2.
    template<typename _Tp>
    void
    __cyl_bessel_k_scaled_steed_temme_cont_frac_2(const _Tp __nu,
                                                  const _Tp __x,
                                                  _Tp & __K_nu,
                                                  _Tp & __K_nup1,
                                                  _Tp & __Kp_nu)
    {
      const int __max_iter = 10000;

      _Tp __bi = _Tp(2) * (_Tp(1) + __x);
      _Tp __di = _Tp(1) / __bi;
      _Tp __delhi = __di;
      _Tp __hi = __di;

      _Tp __qi   = _Tp(0);
      _Tp __qip1 = _Tp(1);

      _Tp __ai = -(_Tp(0.25L) - __nu * __nu);
      _Tp __a1 = __ai;
      _Tp __ci = -__ai;
      _Tp __Qi = -__ai;

      _Tp __s = _Tp(1) + __Qi * __delhi;

      for (int __i = 2; __i <= __max_iter; ++__i)
        {
          __ai -= _Tp(2 * (__i - 1));
          __ci  = -__ai * __ci / __i;
          const _Tp __tmp = (__qi - __bi * __qip1) / __ai;
          __qi = __qip1;
          __qip1 = __tmp;
          __Qi += __ci * __qip1;
          __bi += _Tp(2);
          __di = _Tp(1) / (__bi + __ai * __di);
          __delhi = (__bi * __di - _Tp(1)) * __delhi;
          __hi += __delhi;
          const _Tp __dels = __Qi * __delhi;
          __s += __dels;
          if (std::abs(__dels / __s) < std::numeric_limits<_Tp>::epsilon())
            break;
        }

      __hi *= -__a1;

      __K_nu = std::sqrt(_TR1_M_PI/(_Tp(2) * __x)) / __s;
      __K_nup1 = __K_nu * (__nu + __x + _Tp(0.5L) - __hi) / __x;
      __Kp_nu = - __K_nup1 + __nu / __x * __K_nu;
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
    __cyl_bessel_i_scaled(const _Tp __nu, const _Tp __x)
    {

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument in "
                                      "__cyl_bessel_i_scaled."));

      if (std::isnan(__x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (__x == _Tp(0))
        {
          if (__nu == _Tp(0))
            return _Tp(1);
          else
            return _Tp(0);
        }

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __root3_eps = std::exp(std::log(__eps) / _Tp(3));

      if (__x * __x < _Tp(10) * (__nu + _Tp(1)))
        {
          const _Tp __ex = exp(-__x);
          const _Tp __b = __cyl_bessel_ij_series(__nu, __x, _Tp(1), 100);
          return __b * __ex;
        }
      else if (_Tp(0.5L) / (__nu * __nu + __x * __x) < __root3_eps)
        {
      //  return __cyl_bessel_i_asymp(__nu, __x);
          return __cyl_bessel_i_scaled_asymp_unif(__nu, __x);
        }
      else
        {
          int __Nmu = static_cast<int>(__nu + _Tp(0.5L));
          _Tp __mu = __nu - __Nmu;  // -1/2 <= mu <= 1/2
          _Tp __K_mu, __K_mup1, __Kp_mu;
          _Tp __K_nu, __K_nup1, __K_num1;
          _Tp __I_nu_ratio;

          //  Obtain K_mu, K_mup1
          if (__x < _Tp(2))
            {
              __cyl_bessel_k_scaled_temme(__mu, __x, __K_mu, __K_mup1, __Kp_mu);
            }
          else
            {
              __cyl_bessel_k_scaled_steed_temme_cont_frac_2(__mu, __x, __K_mu, __K_mup1, __Kp_mu);
            }

          //  Recurse forward to obtain K_num1, K_nu
          __K_nu   = __K_mu;
          __K_nup1 = __K_mup1;

          for(int __n = 0; __n < __Nmu; ++__n)
            {
              __K_num1 = __K_nu;
              __K_nu   = __K_nup1;
              __K_nup1 = _Tp(2) * (__mu + _Tp(__n + 1)) * __K_nu / __x + __K_num1;
            }

          //  Calculate I_{nu+1}/I_nu
          __I_nu_ratio = __cyl_bessel_i_cont_frac_1_series(__nu, __x);

          // solve for I_nu
          _Tp __I_nu = _Tp(1) / (__x * (__K_nup1 + __I_nu_ratio * __K_nu));
          __I_nu = __eps * (_Tp(0.5L) * __Nmu + _Tp(2)) * std::abs(__I_nu);
          return __I_nu;
       }
    }


    ///
    ///
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_i(_Tp __nu, _Tp __x)
    {
      if (std::isnan(__nu) || std::isnan(__x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument "
                                      "in __cyl_bessel_i."));

      _Tp __I_scaled = __cyl_bessel_i_scaled(__nu, __x);
      return __I_scaled * std::exp(__x);
    }


    ///
    ///  @brief  Return the asymptotic solution of the modified
    ///          Bessel function \f$ K_{\nu} \f$.
    ///
    template<typename _Tp>
    inline _Tp
    __cyl_bessel_k_asymp(const _Tp __n, const _Tp __x)
    {

      _Tp __sum = __cyl_mod_bessel_a_series(__n, __x, false);

      return std::sqrt(_TR1_M_PI_2 / __x) * std::exp(-__x) * __sum;
    }


    ///
    ///  x >> nu*nu + 1
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_k_scaled_asymp(const _Tp __nu, const _Tp __x)
    {
      const _Tp __mu = _Tp(4) * __nu * __nu;
      const _Tp __mum1 = __mu - _Tp(1);
      const _Tp __mum9 = __mu - _Tp(9);
      const _Tp __pre = std::sqrt(_TR1_M_PI / (_Tp(2) * __x));
      const _Tp __r = __nu / __x;
      const _Tp __K = __pre * (_Tp(1) + __mum1 / (_Tp(8) * __x) + __mum1 * __mum9 / (_Tp(128) * __x * __x));

      return __K;
    }


    ///
    ///  @brief  Return the scaled modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu} \f$.
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_k_scaled(const _Tp __nu, const _Tp __x)
    {

      if (std::isnan(__nu) || std::isnan(__x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument "
                                      "in __cyl_bessel_k_scaled."));

      int __Nmu = static_cast<int>(__nu + _Tp(0.5L));
      _Tp __mu = __nu - __Nmu;  // -1/2 <= mu <= 1/2
      _Tp __K_mu, __K_mup1, __Kp_mu;

      if (__x < _Tp(2))
        __cyl_bessel_k_scaled_temme(__mu, __x, __K_mu, __K_mup1, __Kp_mu);
      else
        __cyl_bessel_k_scaled_steed_temme_cont_frac_2(__mu, __x, __K_mu, __K_mup1, __Kp_mu);

      //  Recurse forward to obtain K_num1, K_nu.

      _Tp __K_nu   = __K_mu;
      _Tp __K_nup1 = __K_mup1;

      for (int __n = 0; __n < __Nmu; ++__n)
        {
          const _Tp __K_num1 = __K_nu;
          __K_nu   = __K_nup1;
          __K_nup1 = _Tp(2) * (__mu + __n + 1) / __x * __K_nu + __K_num1;
        }

      return __K_nu;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu} \f$.
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_k(const _Tp __nu, const _Tp __x)
    {
      if (std::isnan(__nu) || std::isnan(__x))
        return -std::numeric_limits<_Tp>::quiet_NaN();

      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument "
                                      "in __cyl_bessel_k."));

      _Tp __K_scaled = __cyl_bessel_k_scaled(__nu, __x);
      return __K_scaled * std::exp(-__x);
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE
}

#endif // _TR1_MODIFIED_BESSEL_FUNC_TCC
