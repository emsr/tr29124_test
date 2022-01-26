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

/** @file emsr/sf_polygamma.tcc
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

#ifndef _GLIBCXX_BITS_SF_DIGAMMA_TCC
#define _GLIBCXX_BITS_SF_DIGAMMA_TCC 1

#include <limits>
#include <emsr/fp_type_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

// Implementation-space details.
namespace __detail
{
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
    hurwitz_zeta_glob(const _Tp a, const _Tp s)
    {
      _Tp zeta = _Tp(0);

      //  Max e exponent before overflow.
      const _Tp max_binom = std::numeric_limits<_Tp>::max_exponent10
                            * std::log(_Tp(10)) - _Tp(1);

      const unsigned int maxit = 10000;
      for (unsigned int i = 0; i < maxit; ++i)
	{
          bool punt = false;
          _Tp sgn = _Tp(1);
          _Tp term = _Tp(0);
          for (unsigned int j = 0; j <= i; ++j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              _Tp binom = std::lgamma(_Tp(1 + i))
                          - std::lgamma(_Tp(1 + j))
                          - std::lgamma(_Tp(1 + i - j));
#else
              _Tp binom = log_gamma(_Tp(1 + i))
                          - log_gamma(_Tp(1 + j))
                          - log_gamma(_Tp(1 + i - j));
#endif
              if (binom > max_binom)
                {
                  //  This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              binom = std::exp(binom);
              term += sgn * binom * std::pow(_Tp(a + j), -s);
              sgn *= _Tp(-1);
            }
          if (punt)
            break;
          term /= _Tp(i + 1);
          if (std::abs(term / zeta) < std::numeric_limits<_Tp>::epsilon())
            break;
          zeta += term;
        }

      zeta /= s - _Tp(1);

      return zeta;
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
    hurwitz_zeta(const _Tp a, const _Tp s)
    {
      return hurwitz_zeta_glob(a, s);
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
    digamma(const _Tp x)
    {
      ///  @todo Finish me!!!
    }


  /**
   * 
   */
  template<typename _Tp>
    _Tp
    digamma_1(const _Tp x)
    {
      int n = std::nearbyint(x);

      if (x == _Tp(n))
        {
          _Tp digamma = -numeric_constants<_Tp>::euler();
          for (int i = 1; i <= n; ++i )
            digamma += _Tp(1) / i;
          return digamma;
        }
      else
        {
        ///  @todo Finish me!!!
        }
    }


  /**
   *   @brief  Return the polygamma function @f$ polygamma^{(n)}(x) @f$.
   * 
   *   The polygamma function is related to the Hurwitz zeta function:
   *   @f[
   *     psi^{(n)}(x) = (-1)^{n+1} m! \zeta(m+1,x)
   *   @f]
   */
  template<typename _Tp>
    _Tp
    polygamma(const unsigned int n, const _Tp x)
    {
      if (x <= _Tp(0))
        throw_domain_error("polygamma: srgument out of range");
      else if (n == 0)
        return polygamma(x);
      else if (n == 1)
        return digamma_1(x);
      else
        {
          const _Tp hzeta = hurwitz_zeta(_Tp(n + 1), x);
#if _GLIBCXX_USE_C99_MATH_TR1
          const _Tp ln_nfact = std::lgamma(_Tp(n + 1));
#else
          const _Tp ln_nfact = log_gamma(_Tp(n + 1));
#endif
          _Tp result = std::exp(ln_nfact) * hzeta;
          if (n % 2 == 1)
            result = -result;
          return result;
        }
    }
} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _GLIBCXX_BITS_SF_DIGAMMA_TCC
