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
    template <typename Tp>
    void
    __bessel_ik(const Tp __nu, const Tp __x,
                Tp & __I_nu, Tp & __K_nu, Tp & __Ip_nu, Tp & __Kp_nu)
    {

      if (std::isnan(__nu) || std::isnan(__x))
        {
          __I_nu = std::numeric_limits<Tp>::quiet_NaN();
          __K_nu = std::numeric_limits<Tp>::quiet_NaN();
          __Ip_nu = std::numeric_limits<Tp>::quiet_NaN();
          __Kp_nu = std::numeric_limits<Tp>::quiet_NaN();
          return;
        }

      if (__x == Tp(0))
        {
          if (__nu == Tp(0))
            {
              __I_nu = Tp(1);
              __K_nu = std::numeric_limits<Tp>::infinity();
              __Ip_nu = Tp(0);
              __Kp_nu = -std::numeric_limits<Tp>::infinity();
            }
          else
            {
              __I_nu = Tp(0);
              __K_nu = std::numeric_limits<Tp>::infinity();
              //__Ip_nu = ???
              //__Kp_nu = ???
            }
          return;
        }

      if (__x < Tp(0) || __nu < Tp(0))
        __throw_runtime_error(__N("Bad arguments in "
                                  "__bessel_ik."));

      const Tp __eps = std::numeric_limits<Tp>::epsilon();
      const Tp __fp_min = Tp(10) * std::numeric_limits<Tp>::epsilon();
      const int __max_iter = 15000;
      const Tp __x_min = Tp(2);
      const Tp __PI = Tp(3.1415926535897932384626433832795029L);

      const int __nl = static_cast<int>(__nu + Tp(0.5L));

      const Tp __mu = __nu - __nl;
      const Tp __mu2 = __mu * __mu;
      const Tp __xi = Tp(1) / __x;
      const Tp __xi2 = Tp(2) * __xi;
      Tp __h = __nu * __xi;
      if ( __h < __fp_min )
        __h = __fp_min;
      Tp __b = __xi2 * __nu;
      Tp __d = Tp(0);
      Tp __c = __h;
      int __i;
      for ( __i = 1; __i <= __max_iter; ++__i )
        {
          __b += __xi2;
          __d = Tp(1) / (__b + __d);
          __c = __b + Tp(1) / __c;
          const Tp __del = __c * __d;
          __h = __del * __h;
          if (std::abs(__del - Tp(1)) < __eps)
            break;
        }
      if (__i > __max_iter)
        __throw_runtime_error(__N("Argument x too large in __bessel_jn; "
                                  "try asymptotic expansion."));
      Tp __I_nul = __fp_min;
      Tp __Ip_nul = __h * __I_nul;
      Tp __I_nul1 = __I_nul;
      Tp __Ip_nu1 = __Ip_nul;
      Tp __fact = __nu * __xi;
      for (int __l = __nl; __l >= 1; --__l)
        {
          const Tp __I_nutemp = __fact * __I_nul + __Ip_nul;
          __fact -= __xi;
          __Ip_nul = __fact * __I_nutemp + __I_nul;
          __I_nul = __I_nutemp;
        }
      Tp __f = __Ip_nul / __I_nul;
      Tp __K_mu, __K_nu1;
      if (__x < __x_min)
        {
          const Tp __x2 = __x / Tp(2);
          const Tp __pimu = __PI * __mu;
          const Tp __fact = (std::abs(__pimu) < __eps
                            ? Tp(1) : __pimu / std::sin(__pimu));
          Tp __d = -std::log(__x2);
          Tp __e = __mu * __d;
          const Tp __fact2 = (std::abs(__e) < __eps
                            ? Tp(1) : sinh(__e) / __e);
          Tp __gam1, __gam2, __gampl, __gammi;
          __gamma_temme(__mu, __gampl, __gammi, __gam1, __gam2);
          Tp __ff = __fact * (__gam1 * cosh(__e) + __gam2 * __fact2 * __d);
          Tp __sum = __ff;
          __e = std::exp(__e);
          Tp __p = __e / (Tp(2) * __gampl);
          Tp __q = Tp(1) / (Tp(2) * __e * __gammi);
          Tp __c = Tp(1);
          __d = __x2 * __x2;
          Tp __sum1 = __p;
          int __i;
          for (__i = 1; __i <= __max_iter; ++__i)
            {
              __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
              __c *= __d / __i;
              __p /= __i - __mu;
              __q /= __i + __mu;
              const Tp __del = __c * __ff;
              __sum += __del; 
              const Tp __del1 = __c * (__p - __i * __ff);
              __sum1 += __del1;
              if (std::abs(__del) < __eps * std::abs(__sum))
                break;
            }
          if (__i > __max_iter)
            __throw_runtime_error(__N("Bessel k series failed to converge "
                                      "in __bessel_jn."));
          __K_mu = __sum;
          __K_nu1 = __sum1 * __xi2;
        }
      else
        {
          Tp __b = Tp(2) * (Tp(1) + __x);
          Tp __d = Tp(1) / __b;
          Tp __delh = __d;
          Tp __h = __delh;
          Tp __q1 = Tp(0);
          Tp __q2 = Tp(1);
          Tp __a1 = 0.25L - __mu2;
          Tp __q = __c = __a1;
          Tp __a = -__a1;
          Tp __s = Tp(1) + __q * __delh;
          int __i;
          for (__i = 2; __i <= __max_iter; ++__i)
            {
              __a -= 2 * (__i - 1);
              __c = -__a * __c / __i;
              const Tp __qnew = (__q1 - __b * __q2) / __a;
              __q1 = __q2;
              __q2 = __qnew;
              __q += __c * __qnew;
              __b += Tp(2);
              __d = Tp(1) / (__b + __a * __d);
              __delh = (__b * __d - Tp(1)) * __delh;
              __h += __delh;
              const Tp __dels = __q * __delh;
              __s += __dels;
              if ( std::abs(__dels / __s) < __eps )
                break;
            }
          if (__i > __max_iter)
            __throw_runtime_error(__N("Steed's method failed "
                                      "in __bessel_jn."));
          __h = __a1 * __h;
          __K_mu = std::sqrt(__PI / (Tp(2) * __x)) * std::exp(-__x) / __s;
          __K_nu1 = __K_mu * (__mu + __x + Tp(0.5L) - __h) * __xi;
        }

      Tp __K_mup = __mu * __xi * __K_mu - __K_nu1;
      Tp __I_mu = __xi / (__f * __K_mu - __K_mup);
      __I_nu = __I_mu * __I_nul1 / __I_nul;
      __Ip_nu = __I_mu * __Ip_nu1 / __I_nul;
      for ( __i = 1; __i <= __nl; ++__i )
        {
          const Tp __K_nutemp = (__mu + __i) * __xi2 * __K_nu1 + __K_mu;
          __K_mu = __K_nu1;
          __K_nu1 = __K_nutemp;
        }
      __K_nu = __K_mu;
      __Kp_nu = __nu * __xi * __K_mu - __K_nu1;

      return;
    }


    ///
    ///
    ///
    template <typename Tp>
    void
    __airy(const Tp __x,
           Tp & __Ai, Tp & __Bi, Tp & __Aip, Tp & __Bip)
    {
      const Tp __PI = Tp(3.1415926535897932384626433832795029L);
      const Tp __SQRT3 = std::sqrt(Tp(3));
      const Tp __absx = std::abs(__x);
      const Tp __rootx = std::sqrt(__absx);
      const Tp __z = Tp(2) * __absx * __rootx / Tp(3);
      if (__x > Tp(0))
        {
          Tp __I_nu, __Ip_nu, __K_nu, __Kp_nu;

          __bessel_jn(Tp(1)/Tp(3), __z, __I_nu, __K_nu, __Ip_nu, __Kp_nu);
          __Ai = __rootx * __K_nu / (__SQRT3 * __PI);
          __Bi = __rootx * (__K_nu / __PI + Tp(2) * __I_nu / __SQRT3);

          __bessel_jn(Tp(2)/Tp(3), __z, __I_nu, __K_nu, __Ip_nu, __Kp_nu);
          __Aip = -__x * __K_nu / (__SQRT3 * __PI);
          __Bip = __x * (__K_nu / __PI + Tp(2) * __I_nu / __SQRT3);
        }
      else if (__x < Tp(0))
        {
          Tp __J_nu, __Jp_nu, __N_nu, __Np_nu;

          __bessel_jn(Tp(1)/Tp(3), __z, __J_nu, __N_nu, __Jp_nu, __Np_nu);
          __Ai = __rootx * (__J_nu - __N_nu / __SQRT3) / Tp(2);
          __Bi = -__rootx * (__N_nu + __J_nu / __SQRT3) / Tp(2);

          __bessel_jn(Tp(2)/Tp(3), __z, __J_nu, __N_nu, __Jp_nu, __Np_nu);
          __Aip = __absx * (__N_nu / __SQRT3 + __J_nu) / Tp(2);
          __Bip = __absx * (__J_nu / __SQRT3 - __N_nu) / Tp(2);
        }
      else
        {
          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
          // The number is Ai(0) = 3^{-2/3}/\Gamma(2/3).
          __Ai = Tp(0.35502805388781723926L);
          __Bi = __Ai * __SQRT3;

          // Reference:
          //   Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
          // The number is Ai'(0) = -3^{-1/3}/\Gamma(1/3).
          __Aip = -Tp(0.25881940379280679840L);
          __Bip = -__Aip * __SQRT3;
        }

      return;
    }


    ///
    ///
    ///
    template <typename Tp>
    void
    __sph_bessel_ik(const int __n, const Tp __x,
                    Tp & __si, Tp & __sk, Tp & __sip, Tp & __skp)
    {

      if ( __n < 0 || __x < Tp(0) )
        __throw_runtime_error(__N("Bad arguments "
                                  "in sph_bessel."));

      const Tp __nu = Tp(__n) + Tp(0.5L);

      Tp __I_nu, __Ip_nu, __K_nu, __Kp_nu;
      __bessel_ik(__nu, __x, __I_nu, __K_nu, __Ip_nu, __Kp_nu);

      const Tp __SQRTPIO2 = Tp(1.2533141373155002512078826424055226L);
      const Tp __factor = __SQRTPIO2 / std::sqrt(__x);

      __si = __factor * __I_nu;
      __sk = __factor * __K_nu;
      __sip = __factor * __Ip_nu - __si / (Tp(2) * __x);
      __skp = __factor * __Kp_nu - __sk / (Tp(2) * __x);

      return;
    }

  } // namespace std::tr1::__detail

  /* @} */ // group tr1_math_spec_func

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _TR1_BESSEL_FUNCTION_TCC
