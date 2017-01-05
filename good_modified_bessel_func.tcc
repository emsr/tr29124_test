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
    ///  @brief  Compute the modified Bessel functions @f$ I_\nu(x) @f$ and
    ///          @f$ K_\nu(x) @f$ and their first derivatives
    ///          @f$ I'_\nu(x) @f$ and @f$ K'_\nu(x) @f$ respectively.
    ///          These four functions are computed together for numerical
    ///          stability.
    ///
    template <typename _Tp>
    void
    __bessel_ik(const _Tp __nu, const _Tp __x,
                _Tp & __I_nu, _Tp & __K_nu, _Tp & __Ip_nu, _Tp & __Kp_nu)
    {

      if (__x < _Tp(0) || __nu < _Tp(0))
        __throw_runtime_error(__N("Bad arguments in "
                                  "__bessel_ik."));

      if (std::isnan(__nu) || std::isnan(__x))
        {
          __I_nu = std::numeric_limits<_Tp>::quiet_NaN();
          __K_nu = std::numeric_limits<_Tp>::quiet_NaN();
          __Ip_nu = std::numeric_limits<_Tp>::quiet_NaN();
          __Kp_nu = std::numeric_limits<_Tp>::quiet_NaN();
          return;
        }

      if (__x == _Tp(0))
        {
          if (__nu == _Tp(0))
            {
              __I_nu = _Tp(1);
              __K_nu = std::numeric_limits<_Tp>::infinity();
              __Ip_nu = _Tp(0);
              __Kp_nu = -std::numeric_limits<_Tp>::infinity();
            }
          else
            {
              __I_nu = _Tp(0);
              __K_nu = std::numeric_limits<_Tp>::infinity();
              //__Ip_nu = ???
              //__Kp_nu = ???
            }
          return;
        }

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __fp_min = _Tp(10) * std::numeric_limits<_Tp>::epsilon();
      const int __max_iter = 15000;
      const _Tp __x_min = _Tp(2);

      const int __nl = static_cast<int>(__nu + _Tp(0.5L));

      const _Tp __mu = __nu - __nl;
      const _Tp __mu2 = __mu * __mu;
      const _Tp __xi = _Tp(1) / __x;
      const _Tp __xi2 = _Tp(2) * __xi;
      _Tp __h = __nu * __xi;
      if ( __h < __fp_min )
        __h = __fp_min;
      _Tp __b = __xi2 * __nu;
      _Tp __d = _Tp(0);
      _Tp __c = __h;
      int __i;
      for ( __i = 1; __i <= __max_iter; ++__i )
        {
          __b += __xi2;
          __d = _Tp(1) / (__b + __d);
          __c = __b + _Tp(1) / __c;
          const _Tp __del = __c * __d;
          __h = __del * __h;
          if (std::abs(__del - _Tp(1)) < __eps)
            break;
        }
      if (__i > __max_iter)
        __throw_runtime_error(__N("Argument x too large in __bessel_ik; "
                                  "try asymptotic expansion."));
      _Tp __I_nul = __fp_min;
      _Tp __Ip_nul = __h * __I_nul;
      _Tp __I_nul1 = __I_nul;
      _Tp __Ip_nu1 = __Ip_nul;
      _Tp __fact = __nu * __xi;
      for (int __l = __nl; __l >= 1; --__l)
        {
          const _Tp __I_nutemp = __fact * __I_nul + __Ip_nul;
          __fact -= __xi;
          __Ip_nul = __fact * __I_nutemp + __I_nul;
          __I_nul = __I_nutemp;
        }
      _Tp __f = __Ip_nul / __I_nul;
      _Tp __K_mu, __K_nu1;
      if (__x < __x_min)
        {
          const _Tp __x2 = __x / _Tp(2);
          const _Tp __pimu = __PI * __mu;
          const _Tp __fact = (std::abs(__pimu) < __eps
                            ? _Tp(1) : __pimu / std::sin(__pimu));
          _Tp __d = -std::log(__x2);
          _Tp __e = __mu * __d;
          const _Tp __fact2 = (std::abs(__e) < __eps
                            ? _Tp(1) : sinh(__e) / __e);
          _Tp __gam1, __gam2, __gampl, __gammi;
          __gamma_temme(__mu, __gampl, __gammi, __gam1, __gam2);
          _Tp __ff = __fact * (__gam1 * cosh(__e) + __gam2 * __fact2 * __d);
          _Tp __sum = __ff;
          __e = std::exp(__e);
          _Tp __p = __e / (_Tp(2) * __gampl);
          _Tp __q = _Tp(1) / (_Tp(2) * __e * __gammi);
          _Tp __c = _Tp(1);
          __d = __x2 * __x2;
          _Tp __sum1 = __p;
          int __i;
          for (__i = 1; __i <= __max_iter; ++__i)
            {
              __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
              __c *= __d / __i;
              __p /= __i - __mu;
              __q /= __i + __mu;
              const _Tp __del = __c * __ff;
              __sum += __del; 
              const _Tp __del1 = __c * (__p - __i * __ff);
              __sum1 += __del1;
              if (std::abs(__del) < __eps * std::abs(__sum))
                break;
            }
          if (__i > __max_iter)
            __throw_runtime_error(__N("Bessel k series failed to converge "
                                      "in __bessel_ik."));
          __K_mu = __sum;
          __K_nu1 = __sum1 * __xi2;
        }
      else
        {
          _Tp __b = _Tp(2) * (_Tp(1) + __x);
          _Tp __d = _Tp(1) / __b;
          _Tp __delh = __d;
          _Tp __h = __delh;
          _Tp __q1 = _Tp(0);
          _Tp __q2 = _Tp(1);
          _Tp __a1 = 0.25L - __mu2;
          _Tp __q = __c = __a1;
          _Tp __a = -__a1;
          _Tp __s = _Tp(1) + __q * __delh;
          int __i;
          for (__i = 2; __i <= __max_iter; ++__i)
            {
              __a -= 2 * (__i - 1);
              __c = -__a * __c / __i;
              const _Tp __qnew = (__q1 - __b * __q2) / __a;
              __q1 = __q2;
              __q2 = __qnew;
              __q += __c * __qnew;
              __b += _Tp(2);
              __d = _Tp(1) / (__b + __a * __d);
              __delh = (__b * __d - _Tp(1)) * __delh;
              __h += __delh;
              const _Tp __dels = __q * __delh;
              __s += __dels;
              if ( std::abs(__dels / __s) < __eps )
                break;
            }
          if (__i > __max_iter)
            __throw_runtime_error(__N("Steed's method failed "
                                      "in __bessel_ik."));
          __h = __a1 * __h;
          __K_mu = std::sqrt(__PI / (_Tp(2) * __x)) * std::exp(-__x) / __s;
          __K_nu1 = __K_mu * (__mu + __x + _Tp(0.5L) - __h) * __xi;
        }

      _Tp __K_mup = __mu * __xi * __K_mu - __K_nu1;
      _Tp __I_mu = __xi / (__f * __K_mu - __K_mup);
      __I_nu = __I_mu * __I_nul1 / __I_nul;
      __Ip_nu = __I_mu * __Ip_nu1 / __I_nul;
      for ( __i = 1; __i <= __nl; ++__i )
        {
          const _Tp __K_nutemp = (__mu + __i) * __xi2 * __K_nu1 + __K_mu;
          __K_mu = __K_nu1;
          __K_nu1 = __K_nutemp;
        }
      __K_nu = __K_mu;
      __Kp_nu = __nu * __xi * __K_mu - __K_nu1;

      return;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ I_{\nu}(x) \f$.
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_i(const _Tp __nu, const _Tp __x)
    {
      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument "
                                      "in __cyl_bessel_i."));

      if (std::isnan(__nu) || std::isnan(__x))
        return std::numeric_limits<_Tp>::quiet_NaN();

      _Tp __I_nu, __K_nu, __Ip_nu, __Kp_nu;
      __bessel_ik(__nu, __x, __I_nu, __K_nu, __Ip_nu, __Kp_nu);

      return __I_nu;
    }


    ///
    ///  @brief  Return the modified Bessel function of order \f$ \nu \f$:
    ///          \f$ K_{\nu}(x) \f$.
    ///
    template<typename _Tp>
    _Tp
    __cyl_bessel_k(const _Tp __nu, const _Tp __x)
    {
      if (__x < _Tp(0))
        std::__throw_domain_error(__N("Bad argument "
                                      "in __cyl_bessel_k."));

      if (std::isnan(__nu) || std::isnan(__x))
        return std::numeric_limits<_Tp>::quiet_NaN();

      _Tp __I_nu, __K_nu, __Ip_nu, __Kp_nu;
      __bessel_ik(__nu, __x, __I_nu, __K_nu, __Ip_nu, __Kp_nu);

      return __K_nu;
    }


    ///
    ///  @brief  Compute the spherical modified Bessel functions
    ///          @f$ i_n(x) @f$ and @f$ k_n(x) @f$ and their first
    ///          derivatives @f$ i'_n(x) @f$ and @f$ k'_n(x) @f$
    ///          respectively.
    ///
    template <typename _Tp>
    void
    __sph_bessel_ik(const unsigned int __n, const _Tp __x,
                    _Tp & __i_n, _Tp & __k_n, _Tp & __ip_n, _Tp & __kp_n)
    {

      if (__x < _Tp(0))
        __throw_runtime_error(__N("Bad arguments "
                                  "in sph_bessel."));

      const _Tp __nu = _Tp(__n) + _Tp(0.5L);

      _Tp __I_nu, __Ip_nu, __K_nu, __Kp_nu;
      __bessel_ik(__nu, __x, __I_nu, __K_nu, __Ip_nu, __Kp_nu);

      const _Tp __factor = _Tp(__SQRTPIO2) / std::sqrt(__x);

      __i_n = __factor * __I_nu;
      __k_n = __factor * __K_nu;
      __ip_n = __factor * __Ip_nu - __i_n / (_Tp(2) * __x);
      __kp_n = __factor * __Kp_nu - __k_n / (_Tp(2) * __x);

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
    __airy(const _Tp __x,
           _Tp & __Ai, _Tp & __Bi, _Tp & __Aip, _Tp & __Bip)
    {
      const _Tp __absx = std::abs(__x);
      const _Tp __rootx = std::sqrt(__absx);
      const _Tp __z = _Tp(2) * __absx * __rootx / _Tp(3);
      if (__x > _Tp(0))
        {
          _Tp __I_nu, __Ip_nu, __K_nu, __Kp_nu;

          __bessel_jn(_Tp(1)/_Tp(3), __z, __I_nu, __K_nu, __Ip_nu, __Kp_nu);
          __Ai = __rootx * __K_nu / (__SQRT3 * __PI);
          __Bi = __rootx * (__K_nu / __PI + _Tp(2) * __I_nu / __SQRT3);

          __bessel_jn(_Tp(2)/_Tp(3), __z, __I_nu, __K_nu, __Ip_nu, __Kp_nu);
          __Aip = -__x * __K_nu / (__SQRT3 * __PI);
          __Bip = __x * (__K_nu / __PI + _Tp(2) * __I_nu / __SQRT3);
        }
      else if (__x < _Tp(0))
        {
          _Tp __J_nu, __Jp_nu, __N_nu, __Np_nu;

          __bessel_jn(_Tp(1)/_Tp(3), __z, __J_nu, __N_nu, __Jp_nu, __Np_nu);
          __Ai = __rootx * (__J_nu - __N_nu / __SQRT3) / _Tp(2);
          __Bi = -__rootx * (__N_nu + __J_nu / __SQRT3) / _Tp(2);

          __bessel_jn(_Tp(2)/_Tp(3), __z, __J_nu, __N_nu, __Jp_nu, __Np_nu);
          __Aip = __absx * (__N_nu / __SQRT3 + __J_nu) / _Tp(2);
          __Bip = __absx * (__J_nu / __SQRT3 - __N_nu) / _Tp(2);
        }
      else
        {
          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
          // The number is Ai(0) = 3^{-2/3}/\Gamma(2/3).
          __Ai = _Tp(0.35502805388781723926L);
          __Bi = __Ai * __SQRT3;

          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
          // The number is Ai'(0) = -3^{-1/3}/\Gamma(1/3).
          __Aip = -_Tp(0.25881940379280679840L);
          __Bip = -__Aip * __SQRT3;
        }

      return;
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE
}

#endif // _TR1_MODIFIED_BESSEL_FUNC_TCC
