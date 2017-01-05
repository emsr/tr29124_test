// Special functions -*- C++ -*-

// Copyright (C) 2006-2017 Free Software Foundation, Inc.
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

/** @file bits/sf_psi.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
//   (1) Handbook of Mathematical Functions,
//       ed. Milton Abramowitz and Irene A. Stegun,
//       Dover Publications,
//       Section 9, pp. 355-434, Section 10 pp. 435-478
//   (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
//   (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//       W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//       2nd ed, pp. 240-245

#ifndef _GLIBCXX_BITS_SF_PSI_TCC
#define _GLIBCXX_BITS_SF_PSI_TCC 1

#include <limits>
//#include <bits/specfun_util.h>
#include "specfun_util.h"

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief  Return the Hurwitz zeta function @f$ \zeta(x,s) @f$
   *  	for all s != 1 and x > -1.
   * 
   *  The Hurwitz zeta function is defined by:
   *  @f[
   *    \zeta(x,s) = \sum_{n=0}^{\infty} \frac{1}{(n + x)^s}
   *  @f]
   *  The Riemann zeta function is a special case:
   *  @f[
   *    \zeta(s) = \zeta(1,s)
   *  @f]
   * 
   *  This functions uses the double sum that converges for s != 1
   *  and x > -1:
   *  @f[
   *    \zeta(x,s) = \frac{1}{s-1}
   *  	     \sum_{n=0}^{\infty} \frac{1}{n + 1}
   *  	     \sum_{k=0}^{n} (-1)^k \frac{n!}{(n-k)!k!} (x+k)^{-s}
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __hurwitz_zeta_glob(const _Tp __a, const _Tp __s)
    {
      _Tp __zeta = _Tp(0);

      //  Max e exponent before overflow.
      const _Tp __max_bincoeff = std::numeric_limits<_Tp>::max_exponent10
                               * std::log(_Tp(10)) - _Tp(1);

      const unsigned int __maxit = 10000;
      for (unsigned int __i = 0; __i < __maxit; ++__i)
	{
          bool __punt = false;
          _Tp __sgn = _Tp(1);
          _Tp __term = _Tp(0);
          for (unsigned int __j = 0; __j <= __i; ++__j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              _Tp __bincoeff =  std::lgamma(_Tp(1 + __i))
                              - std::lgamma(_Tp(1 + __j))
                              - std::lgamma(_Tp(1 + __i - __j));
#else
              _Tp __bincoeff =  __log_gamma(_Tp(1 + __i))
                              - __log_gamma(_Tp(1 + __j))
                              - __log_gamma(_Tp(1 + __i - __j));
#endif
              if (__bincoeff > __max_bincoeff)
                {
                  //  This only gets hit for x << 0.
                  __punt = true;
                  break;
                }
              __bincoeff = std::exp(__bincoeff);
              __term += __sgn * __bincoeff * std::pow(_Tp(__a + __j), -__s);
              __sgn *= _Tp(-1);
            }
          if (__punt)
            break;
          __term /= _Tp(__i + 1);
          if (std::abs(__term / __zeta) < std::numeric_limits<_Tp>::epsilon())
            break;
          __zeta += __term;
        }

      __zeta /= __s - _Tp(1);

      return __zeta;
    }


  /**
   *   @brief  Return the Hurwitz zeta function @f$ \zeta(x,s) @f$
   *           for all s != 1 and x > -1.
   * 
   *   The Hurwitz zeta function is defined by:
   *   @f[
   *     \zeta(x,s) = \sum_{n=0}^{\infty} \frac{1}{(n + x)^s}
   *   @f]
   *   The Riemann zeta function is a special case:
   *   @f[
   *     \zeta(s) = \zeta(1,s)
   *   @f]
   */
  template<typename _Tp>
    _Tp
    __hurwitz_zeta(const _Tp __a, const _Tp __s)
    {
      return __hurwitz_zeta_glob(__a, __s);
    }


  /**
   *   @brief  Return the digamma function.
   *   The digamma or @f$ psi(x) @f$ function is defined by
   *   @f[
   *     psi(x) = \frac{Gamma'(x)}{\Gamma(x)}
   *   @f]
   */
  template<typename _Tp>
    _Tp
    __psi(const _Tp __x)
    {
      ///  @todo Finish me!!!
    }


  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __psi_1(const _Tp __x)
    {
      int __n = std::nearbyint(__x);

      if (__x == _Tp(__n))
        {
          _Tp __psi = -__numeric_constants<_Tp>::euler();
          for (int __i = 1; __i <= __n; ++__i )
            __psi += _Tp(1) / __i;
          return __psi;
        }
      else
        {
        ///  @todo Finish me!!!
        }
    }


  /**
   *   @brief  Return the polygamma function @f$ psi^{(n)}(x) @f$.
   * 
   *   The polygamma function is related to the Hurwitz zeta function:
   *   @f[
   *     psi^{(n)}(x) = (-1)^{n+1} m! \zeta(m+1,x)
   *   @f]
   */
  template<typename _Tp>
    _Tp
    __psi(const unsigned int __n, const _Tp __x)
    {
      if (__x <= _Tp(0))
        __throw_domain_error("__psi: srgument out of range");
      else if (__n == 0)
        return __psi(__x);
      else if (__n == 1)
        return __psi_1(__x);
      else
        {
          const _Tp __hzeta = __hurwitz_zeta(_Tp(__n + 1), __x);
#if _GLIBCXX_USE_C99_MATH_TR1
          const _Tp __ln_nfact = std::lgamma(_Tp(__n + 1));
#else
          const _Tp __ln_nfact = __log_gamma(_Tp(__n + 1));
#endif
          _Tp __result = std::exp(__ln_nfact) * __hzeta;
          if (__n % 2 == 1)
            __result = -__result;
          return __result;
        }
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_BITS_SF_PSI_TCC
