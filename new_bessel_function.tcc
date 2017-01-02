// Special functions -*- C++ -*-

// Copyright (C) 2006-2017
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
    ///
    ///
    template <typename _Tp>
    void
    __bessel_jn(const _Tp __nu, const _Tp __x,
                _Tp & __J_nu, _Tp & __N_nu, _Tp & __Jp_nu, _Tp & __Np_nu)
    {

      if (std::isnan(__nu) || std::isnan(__x))
        {
          __J_nu = std::numeric_limits<_Tp>::quiet_NaN();
          __N_nu = std::numeric_limits<_Tp>::quiet_NaN();
          __Jp_nu = std::numeric_limits<_Tp>::quiet_NaN();
          __Np_nu = std::numeric_limits<_Tp>::quiet_NaN();
          return;
        }

      if (__x == _Tp(0))
        {
          if (__nu == _Tp(0))
            {
              __J_nu = _Tp(1);
              __N_nu = -std::numeric_limits<_Tp>::infinity();
              __Jp_nu = _Tp(0);
              __Np_nu = std::numeric_limits<_Tp>::infinity();
            }
          else
            {
              __J_nu = _Tp(0);
              __N_nu = -std::numeric_limits<_Tp>::infinity();
              //__Jp_nu = ???
              //__Np_nu = ???
            }
          return;
        }

      const _Tp __eps = std::numeric_limits<_Tp>::epsilon();
      const _Tp __fp_min = _Tp(10) * std::numeric_limits<_Tp>::min();
      const int __max_iter = 15000;
      const _Tp __x_min = _Tp(2);
      const _Tp __PI = _Tp(3.1415926535897932384626433832795029L);

      if (__x < _Tp(0) || __nu < _Tp(0))
        __throw_runtime_error(__N("Bad arguments in __bessel_jn."));

      const int __nl = (__x < __x_min
                    ? static_cast<int>(__nu + _Tp(0.5L))
                    : std::max(0, static_cast<int>(__nu - __x + _Tp(1.5L))));

      const _Tp __mu = __nu - __nl;
      const _Tp __mu2 = __mu * __mu;
      const _Tp __xi = _Tp(1) / __x;
      const _Tp __xi2 = _Tp(2) * __xi;
      _Tp __w = __xi2 / __PI;
      int __isign = 1;
      _Tp __h = __nu * __xi;
      if (__h < __fp_min)
        __h = __fp_min;
      _Tp __b = __xi2 * __nu;
      _Tp __d = _Tp(0);
      _Tp __c = __h;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
        {
          __b += __xi2;
          __d = __b - __d;
          if (std::abs(__d) < __fp_min)
            __d = __fp_min;
          __c = __b - _Tp(1) / __c;
          if (std::abs(__c) < __fp_min)
            __c = __fp_min;
          __d = _Tp(1) / __d;
          const _Tp __del = __c * __d;
          __h = __del * __h;
          if (__d < _Tp(0))
            __isign = -__isign;
          if (std::abs(__del - _Tp(1)) < __eps)
            break;
        }
      if (__i > __max_iter)
        __throw_runtime_error(__N("Argument x too large in __bessel_jn; "
                                  "try asymptotic expansion."));
      _Tp __J_nul = __isign * __fp_min;
      _Tp __Jp_nul = __h * __J_nul;
      _Tp __J_nul1 = __J_nul;
      _Tp __Jp_nu1 = __Jp_nul;
      _Tp __fact = __nu * __xi;
      for ( int __l = __nl; __l >= 1; --__l )
        {
          const _Tp __J_nutemp = __fact * __J_nul + __Jp_nul;
          __fact -= __xi;
          __Jp_nul = __fact * __J_nutemp - __J_nul;
          __J_nul = __J_nutemp;
        }
      if (__J_nul == _Tp(0))
        __J_nul = __eps;
      _Tp __f= __Jp_nul / __J_nul;
      _Tp __N_mu, __N_nu1, __N_mup, __J_mu;
      if (__x < __x_min)
        {
          const _Tp __x2 = __x / _Tp(2);
          const _Tp __pimu = __PI * __mu;
          _Tp __fact = (std::abs(__pimu) < __eps ? _Tp(1)
                      : __pimu / std::sin(__pimu));
          _Tp __d = -std::log(__x2);
          _Tp __e = __mu * __d;
          _Tp __fact2 = (std::abs(__e) < __eps
                       ? _Tp(1) : std::sinh(__e) / __e);
          _Tp __gam1, __gam2, __gampl, __gammi;
          __gamma_temme(__mu, __gampl, __gammi, __gam1, __gam2);
          _Tp __ff = (_Tp(2) / __PI) * __fact
                   * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
          __e = std::exp(__e);
          _Tp __p = __e / (__PI * __gampl);
          _Tp __q = _Tp(1) / (__e * __PI * __gammi);
          const _Tp __pimu2 = __pimu / _Tp(2);
          _Tp __fact3 = (std::abs(__pimu2) < __eps
                       ? _Tp(1) : std::sin(__pimu2) / __pimu2 );
          _Tp __r = __PI * __pimu2 * __fact3 * __fact3;
          _Tp __c = _Tp(1);
          __d = -__x2 * __x2;
          _Tp __sum = __ff + __r * __q;
          _Tp __sum1 = __p;
          for (__i = 1; __i <= __max_iter; ++__i)
            {
              __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
              __c *= __d / _Tp(__i);
              __p /= _Tp(__i) - __mu;
              __q /= _Tp(__i) + __mu;
              const _Tp __del = __c * (__ff + __r * __q);
              __sum += __del; 
              const _Tp __del1 = __c * __p - __i * __del;
              __sum1 += __del1;
              if ( std::abs(__del) < __eps * (_Tp(1) + std::abs(__sum)) )
                break;
            }
          if ( __i > __max_iter )
            __throw_runtime_error(__N("Bessel y series failed to converge "
                                      "in __bessel_jn."));
          __N_mu = -__sum;
          __N_nu1 = -__sum1 * __xi2;
          __N_mup = __mu * __xi * __N_mu - __N_nu1;
          __J_mu = __w / (__N_mup - __f * __N_mu);
        }
      else
        {
          _Tp __a = _Tp(0.25L) - __mu2;
          _Tp __q = _Tp(1);
          _Tp __p = -__xi / _Tp(2);
          _Tp __br = _Tp(2) * __x;
          _Tp __bi = _Tp(2);
          _Tp __fact = __a * __xi / (__p * __p + __q * __q);
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
          int __i;
          for (__i = 2; __i <= __max_iter; ++__i)
            {
              __a += _Tp(2 * (__i - 1));
              __bi += _Tp(2);
              __dr = __a * __dr + __br;
              __di = __a * __di + __bi;
              if (std::abs(__dr) + std::abs(__di) < __fp_min)
                __dr = __fp_min;
              __fact = __a / (__cr * __cr + __ci * __ci);
              __cr = __br + __cr * __fact;
              __ci = __bi - __ci * __fact;
              if (std::abs(__cr) + std::abs(__ci) < __fp_min)
                __cr = __fp_min;
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
          if (__i > __max_iter)
            __throw_runtime_error(__N("Lentz's method failed "
                                      "in __bessel_jn."));
          const _Tp __gam = (__p - __f) / __q;
          __J_mu = std::sqrt(__w / ((__p - __f) * __gam + __q));
          __J_mu = ::copysign( __J_mu, __J_nul );
          __N_mu = __J_mu * __gam;
          __N_mup = __N_mu * (__p + __q / __gam);
          __N_nu1 = __mu * __xi * __N_mu - __N_mup;
      }
      __fact = __J_mu / __J_nul;
      __J_nu = __J_nul1 * __fact;
      __Jp_nu = __Jp_nu1 * __fact;
      for (__i = 1; __i <= __nl; ++__i)
        {
          const _Tp __N_nutemp = (__mu + __i) * __xi2 * __N_nu1 - __N_mu;
          __N_mu = __N_nu1;
          __N_nu1 = __N_nutemp;
        }
      __N_nu = __N_mu;
      __Np_nu = __nu * __xi * __N_mu - __N_nu1;

      return;
    }


    ///
    ///
    ///
    template <typename _Tp>
    void
    __sph_bessel_jn(const int __n, const _Tp __x,
                    _Tp & __j_n, _Tp & __n_n, _Tp & __jp_n, _Tp & __np_n)
    {

      if ( __n < 0 || __x < _Tp(0) )
        __throw_runtime_error(__N("Bad arguments in sph_bessel."));

      const _Tp __nu = _Tp(__n) + _Tp(0.5L);

      _Tp __J_nu, __Jp_nu, __N_nu, __Np_nu;
      __bessel_jn(__nu, __x, __J_nu, __N_nu, __Jp_nu, __Np_nu);

      const _Tp __SQRTPIO2 = _Tp(1.2533141373155002512078826424055226L);
      const _Tp __factor = __SQRTPIO2 / std::sqrt(__x);

      __j_n = __factor * __J_nu;
      __n_n = __factor * __N_nu;
      __jp_n = __factor * __Jp_nu - __j_n / (_Tp(2) * __x);
      __np_n = __factor * __Np_nu - __n_n / (_Tp(2) * __x);

      return;
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE
}

#endif // _TR1_BESSEL_FUNCTION_TCC

