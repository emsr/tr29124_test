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

/** @file tr1/riemann_zeta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{tr1/cmath}
 */

//
// ISO C++ 14882 TR1: 5.2  Special functions
//

// Written by Edward Smith-Rowland based on:
//   (1) Handbook of Mathematical Functions,
//       Ed. by Milton Abramowitz and Irene A. Stegun,
//       Dover Publications, New-York, Section 5, pp. 807-808.
//   (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
//   (3) Gamma, Exploring Euler's Constant, Julian Havil,
//       Princeton, 2003.

#ifndef _GLIBCXX_TR1_RIEMANN_ZETA_TCC
#define _GLIBCXX_TR1_RIEMANN_ZETA_TCC 1

#include <special_function_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

#if __STDCPP_WANT_MATH_SPEC_FUNCS__
# define _GLIBCXX_MATH_NS ::std
#elif defined(_GLIBCXX_TR1_CMATH)
namespace tr1
{
# define _GLIBCXX_MATH_NS ::std::tr1
#else
# error do not include this header directly, use <cmath> or <tr1/cmath>
#endif
  // [5.2] Special functions

  // Implementation-space details.
  namespace __detail
  {
    /**
     *   @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
     *           by summation for s > 1.
     * 
     *   The Riemann zeta function is defined by:
     *    \f[
     *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
     *    \f]
     *   For s < 1 use the reflection formula:
     *    \f[
     *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
     *    \f]
     */
    template<typename Tp>
    Tp
    riemann_zeta_sum(Tp s)
    {
      //  A user shouldn't get to this.
      if (s < Tp(1))
        throw std::domain_error("Bad argument in zeta sum.");

      const unsigned int max_iter = 10000;
      Tp zeta = Tp(0);
      for (unsigned int k = 1; k < max_iter; ++k)
        {
          Tp term = std::pow(static_cast<Tp>(k), -s);
          if (term < std::numeric_limits<Tp>::epsilon())
            {
              break;
            }
          zeta += term;
        }

      return zeta;
    }


    /**
     *   @brief  Evaluate the Riemann zeta function @f$ \zeta(s) @f$
     *           by an alternate series for s > 0.
     * 
     *   The Riemann zeta function is defined by:
     *    \f[
     *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
     *    \f]
     *   For s < 1 use the reflection formula:
     *    \f[
     *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
     *    \f]
     */
    template<typename Tp>
    Tp
    riemann_zeta_alt(Tp s)
    {
      Tp sgn = Tp(1);
      Tp zeta = Tp(0);
      for (unsigned int i = 1; i < 10000000; ++i)
        {
          Tp term = sgn / std::pow(i, s);
          if (std::abs(term) < std::numeric_limits<Tp>::epsilon())
            break;
          zeta += term;
          sgn *= Tp(-1);
        }
      zeta /= Tp(1) - std::pow(Tp(2), Tp(1) - s);

      return zeta;
    }


    /**
     *   @brief  Evaluate the Riemann zeta function by series for all s != 1.
     *           Convergence is great until largish negative numbers.
     *           Then the convergence of the > 0 sum gets better.
     *
     *   The series is:
     *    \f[
     *      \zeta(s) = \frac{1}{1-2^{1-s}}
     *                 \sum_{n=0}^{\infty} \frac{1}{2^{n+1}}
     *                 \sum_{k=0}^{n} (-1)^k \frac{n!}{(n-k)!k!} (k+1)^{-s}
     *    \f]
     *   Havil 2003, p. 206.
     *
     *   The Riemann zeta function is defined by:
     *    \f[
     *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
     *    \f]
     *   For s < 1 use the reflection formula:
     *    \f[
     *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
     *    \f]
     */
    template<typename Tp>
    Tp
    riemann_zeta_glob(Tp s)
    {
      Tp zeta = Tp(0);

      const Tp eps = std::numeric_limits<Tp>::epsilon();
      //  Max e exponent before overflow.
      const Tp max_binom = std::exp(std::numeric_limits<Tp>::max_exponent10
                            * std::log(Tp(10)) - Tp(1));

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (s < Tp(0))
        {
#if _GLIBCXX_USE_C99_MATH_TR1
          if (_GLIBCXX_MATH_NS::fmod(s, Tp(2)) == Tp(0))
            return Tp(0);
          else
#endif
            {
              Tp zeta = riemann_zeta_glob(Tp(1) - s);
              zeta *= std::pow(Tp(2)
                     * numeric_constants<Tp>::pi(), s)
                     * std::sin(numeric_constants<Tp>::pi_2() * s)
#if _GLIBCXX_USE_C99_MATH_TR1
                     * std::exp(_GLIBCXX_MATH_NS::lgamma(Tp(1) - s))
#else
                     * std::exp(log_gamma(Tp(1) - s))
#endif
                     / numeric_constants<Tp>::pi();
              return zeta;
            }
        }

      Tp num = Tp(0.25L);
      const unsigned int maxit = 10000;
      zeta = Tp(0.5L);
      // This for loop starts at 1 because we already calculated the
      // value of the zeroeth order in zeta above
      for (unsigned int i = 1; i < maxit; ++i)
        {
          bool punt = false;
          Tp term = Tp(1.0L);
          Tp binom = Tp(1.0L);
          // This for loop starts at 1 because we already calculated the value
          // of the zeroeth order in term above.
          for (unsigned int j = 1; j <= i; ++j)
            {
              Tp incr = Tp(i - j + 1) / Tp(j);
              binom *= -incr;
              if(std::abs(binom) > max_binom )
                {
                  // This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              term += binom * std::pow(Tp(1 + j), -s);
            }
          if (punt)
            break;
          term *= num;
          zeta += term;
          if (std::abs(term/zeta) < eps)
            break;
          num *= Tp(0.5L);
        }

      zeta /= Tp(1) - std::pow(Tp(2), Tp(1) - s);

      return zeta;
    }


    /**
     *   @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
     *           using the product over prime factors.
     *    \f[
     *      \zeta(s) = \Pi_{i=1}^\infty \frac{1}{1 - p_i^{-s}}
     *    \f]
     *    where @f$ {p_i} @f$ are the prime numbers.
     * 
     *   The Riemann zeta function is defined by:
     *    \f[
     *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
     *    \f]
     *   For s < 1 use the reflection formula:
     *    \f[
     *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
     *    \f]
     */
    template<typename Tp>
    Tp
    riemann_zeta_product(Tp s)
    {
      static const Tp prime[] = {
        Tp(2), Tp(3), Tp(5), Tp(7), Tp(11), Tp(13), Tp(17), Tp(19),
        Tp(23), Tp(29), Tp(31), Tp(37), Tp(41), Tp(43), Tp(47),
        Tp(53), Tp(59), Tp(61), Tp(67), Tp(71), Tp(73), Tp(79),
        Tp(83), Tp(89), Tp(97), Tp(101), Tp(103), Tp(107), Tp(109)
      };
      static const unsigned int num_primes = sizeof(prime) / sizeof(Tp);

      Tp zeta = Tp(1);
      for (unsigned int i = 0; i < num_primes; ++i)
        {
          const Tp fact = Tp(1) - std::pow(prime[i], -s);
          zeta *= fact;
          if (Tp(1) - fact < std::numeric_limits<Tp>::epsilon())
            break;
        }

      zeta = Tp(1) / zeta;

      return zeta;
    }


    /**
     *   @brief  Return the Riemann zeta function @f$ \zeta(s) @f$.
     * 
     *   The Riemann zeta function is defined by:
     *    \f[
     *      \zeta(s) = \sum_{k=1}^{\infty} k^{-s} for s > 1
     *                 \frac{(2\pi)^s}{pi} sin(\frac{\pi s}{2})
     *                 \Gamma (1 - s) \zeta (1 - s) for s < 1
     *    \f]
     *   For s < 1 use the reflection formula:
     *    \f[
     *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
     *    \f]
     */
    template<typename Tp>
    Tp
    riemann_zeta(Tp s)
    {
      if (std::isnan(s))
        return std::numeric_limits<Tp>::quiet_NaN();
      else if (s == Tp(1))
        return std::numeric_limits<Tp>::infinity();
      else if (s < -Tp(19))
        {
          Tp zeta = riemann_zeta_product(Tp(1) - s);
          zeta *= std::pow(Tp(2) * numeric_constants<Tp>::pi(), s)
                 * std::sin(numeric_constants<Tp>::pi_2() * s)
#if _GLIBCXX_USE_C99_MATH_TR1
                 * std::exp(_GLIBCXX_MATH_NS::lgamma(Tp(1) - s))
#else
                 * std::exp(log_gamma(Tp(1) - s))
#endif
                 / numeric_constants<Tp>::pi();
          return zeta;
        }
      else if (s < Tp(20))
        {
          //  Global double sum or McLaurin?
          bool glob = true;
          if (glob)
            return riemann_zeta_glob(s);
          else
            {
              if (s > Tp(1))
                return riemann_zeta_sum(s);
              else
                {
                  Tp zeta = std::pow(Tp(2)
                                * numeric_constants<Tp>::pi(), s)
                         * std::sin(numeric_constants<Tp>::pi_2() * s)
#if _GLIBCXX_USE_C99_MATH_TR1
                             * _GLIBCXX_MATH_NS::tgamma(Tp(1) - s)
#else
                             * std::exp(log_gamma(Tp(1) - s))
#endif
                             * riemann_zeta_sum(Tp(1) - s);
                  return zeta;
                }
            }
        }
      else
        return riemann_zeta_product(s);
    }


    /**
     *   @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
     *           for all s != 1 and a > -1.
     * 
     *   The Hurwitz zeta function is defined by:
     *   @f[
     *     \zeta(s,a) = \sum_{n=0}^{\infty} \frac{1}{(n + a)^s}
     *   @f]
     *   The Riemann zeta function is a special case:
     *   @f[
     *     \zeta(s) = \zeta(s,1)
     *   @f]
     * 
     *   This functions uses the double sum that converges for s != 1
     *   and a > -1:
     *   @f[
     *     \zeta(s,a) = \frac{1}{s-1}
     *                \sum_{n=0}^{\infty} \frac{1}{n + 1}
     *                \sum_{k=0}^{n} (-1)^k \frac{n!}{(n-k)!k!} (k+a)^{-s}
     *   @f]
     */
    template<typename Tp>
    Tp
    hurwitz_zeta_glob(Tp s, Tp a)
    {
      Tp zeta = Tp(0);

      const Tp eps = std::numeric_limits<Tp>::epsilon();
      //  Max e exponent before overflow.
      const Tp max_binom = std::numeric_limits<Tp>::max_exponent10
                            * std::log(Tp(10)) - Tp(1);

      const unsigned int maxit = 10000;
      for (unsigned int i = 0; i < maxit; ++i)
        {
          bool punt = false;
          Tp sgn = Tp(1);
          Tp term = Tp(0);
          for (unsigned int j = 0; j <= i; ++j)
            {
#if _GLIBCXX_USE_C99_MATH_TR1
              Tp binom = _GLIBCXX_MATH_NS::lgamma(Tp(1 + i))
                          - _GLIBCXX_MATH_NS::lgamma(Tp(1 + j))
                          - _GLIBCXX_MATH_NS::lgamma(Tp(1 + i - j));
#else
              Tp binom = log_gamma(Tp(1 + i))
                          - log_gamma(Tp(1 + j))
                          - log_gamma(Tp(1 + i - j));
#endif
              if (binom > max_binom)
                {
                  //  This only gets hit for x << 0.
                  punt = true;
                  break;
                }
              binom = std::exp(binom);
              term += sgn * binom * std::pow(Tp(j + a), -s);
              sgn *= Tp(-1);
            }
          if (punt)
            break;
          term /= Tp(i + 1);
          if (std::abs(term / zeta) < eps)
            break;
          zeta += term;
        }

      zeta /= s - Tp(1);

      return zeta;
    }


    /**
     *   @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
     *           for all s != 1 and a > -1.
     * 
     *   The Hurwitz zeta function is defined by:
     *   @f[
     *     \zeta(s,a) = \sum_{n=0}^{\infty} \frac{1}{(n + a)^s}
     *   @f]
     *   The Riemann zeta function is a special case:
     *   @f[
     *     \zeta(s) = \zeta(s,1)
     *   @f]
     */
    template<typename Tp>
    Tp
    hurwitz_zeta(Tp s, Tp a)
    { return hurwitz_zeta_glob(s, a); }
  } // namespace __detail
#undef _GLIBCXX_MATH_NS
#if ! __STDCPP_WANT_MATH_SPEC_FUNCS__ && defined(_GLIBCXX_TR1_CMATH)
} // namespace tr1
#endif

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _GLIBCXX_TR1_RIEMANN_ZETA_TCC
