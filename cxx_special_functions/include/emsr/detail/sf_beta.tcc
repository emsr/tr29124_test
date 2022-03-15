
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

/** @file bits/sf_beta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 6, pp. 253-266
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 213-216
// (4) Gamma, Exploring Euler's Constant, Julian Havil,
//     Princeton, 2003.

#ifndef SF_BETA_TCC
#define SF_BETA_TCC 1

#include <stdexcept>

#include <emsr/sf_gamma.h>

namespace emsr
{
namespace detail
{
  /**
   * @brief  Return the beta function: @f$ B(a,b) @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param a The first argument of the beta function.
   * @param b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename Tp>
    Tp
    beta_gamma(Tp a, Tp b)
    {
      auto na = int(std::nearbyint(a));
      auto nb = int(std::nearbyint(b));
      auto nab = int(std::nearbyint(a + b));
      if (nab <= 0 && nab == a + b)
	{
	  if (na != a || na > 0)
	    return Tp{0};
	  else if (nb != b || nb > 0)
	    return Tp{0};
	  else
	    return emsr::quiet_NaN<Tp>();
	}
      else
	{
	  Tp bet;
	  if (std::abs(b) > std::abs(a))
	    {
	      bet = gamma(b) / gamma(a + b);
	      bet *= gamma(a);
	    }
	  else
	    {
	      bet = gamma(a) / gamma(a + b);
	      bet *= gamma(b);
	    }

	  return bet;
	}
    }

  /**
   * @brief  Return the beta function @f$B(a,b)@f$ using
   *	     the log gamma functions.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param a The first argument of the beta function.
   * @param b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename Tp>
    Tp
    beta_lgamma(Tp a, Tp b)
    {
      auto na = int(std::nearbyint(a));
      auto nb = int(std::nearbyint(b));
      auto nab = int(std::nearbyint(a + b));
      if (nab <= 0 && nab == a + b)
	{
	  if (na != a || na > 0)
	    return Tp{0};
	  else if (nb != b || nb > 0)
	    return Tp{0};
	  else
	    return emsr::quiet_NaN<Tp>(a);
	}
      else
	{
	  auto bet = log_gamma(a)
		     + log_gamma(b)
		     - log_gamma(a + b);
	  auto sign = log_gamma_sign(a)
		      * log_gamma_sign(b)
		      * log_gamma_sign(a + b);

	  if (bet > emsr::log_max<Tp>())
            return sign * emsr::infinity<Tp>(a);
	  else
	    return sign * std::exp(bet);
	}
    }


  /**
   * @brief  Return the beta function @f$B(x,y)@f$ using
   *	     the product form.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   * Here, we employ the product form:
   * @f[
   *   B(a,b) = \frac{a + b}{a b} \prod_{k=1}^{\infty}
   *            \frac{1 + (a + b) / k}{(1 + a / k) (1 + b / k)}
   *          = \frac{a + b}{ab} \prod_{k=1}^{\infty}
   *            \left[1 - \frac{ab}{(a + k)(b + k)}\right]
   * @f]
   *
   * @param a The first argument of the beta function.
   * @param b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename Tp>
    Tp
    beta_product(Tp a, Tp b)
    {
      const auto s_eps = emsr::epsilon<Tp>();
      const auto ab = a * b;
      auto bet = (a + b) / ab;

      const unsigned int s_max_iter = 100000; // Could need 1 / sqrt(s_eps)

      const auto apk = a, apb = b;
      for (unsigned int k = 1; k < s_max_iter; ++k)
	{
	  auto term = Tp{1} - ab / (++apk) / (++apb);
	  bet *= term;
	  if (std::abs(Tp{1} - term) < s_eps)
	    break;
	}

      return bet;
    }


  /**
   * @brief  Return the beta function @f$ B(a,b) @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param a The first argument of the beta function.
   * @param b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename Tp>
    Tp
    beta(Tp a, Tp b)
    {
      if (std::isnan(a) || std::isnan(b))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(a) < s_num_factorials<Tp>
	    && std::abs(b) < s_num_factorials<Tp>
	    && std::abs(a + b) < s_num_factorials<Tp>)
	return beta_gamma(a, b);
      else
	return beta_lgamma(a, b);
    }


  /**
   * Return the regularized incomplete beta function, @f$ I_x(a,b) @f$,
   * of arguments @c a, @c b, and @c x.
   *
   *
   * @param a The first parameter
   * @param b The second parameter
   * @param x The argument
   */
  template<typename Tp>
    Tp
    ibeta_cont_frac(Tp a, Tp b, Tp x)
    {
      constexpr unsigned int s_itmax = 100;
      const auto s_fpmin = 1000 * emsr::lim_min<Tp>();
      const auto s_eps = emsr::epsilon<Tp>();

      auto apb = a + b;
      auto ap1 = a + Tp{1};
      auto am1 = a - Tp{1};
      auto c = Tp{1};
      auto d = Tp{1} - apb * x / ap1;
      if (std::abs(d) < s_fpmin)
	d = s_fpmin;
      d = Tp{1} / d;
      auto h = d;
      for (unsigned int m = 1; m <= s_itmax; ++m)
	{
	  auto m2 = 2 * m;

	  //  Even step of the recurrence.
	  auto aa = Tp(m) * (b - Tp(m)) * x
		    / ((am1 + Tp(m2)) * (a + Tp(m2)));
	  d = Tp{1} + aa * d;
	  if (std::abs(d) < s_fpmin)
	    d = s_fpmin;
	  c = Tp{1} + aa / c;
	  if (std::abs(c) < s_fpmin)
	    c = s_fpmin;
	  d = Tp{1} / d;
	  h *= d * c;

	  //  Odd step of the recurrence.
	  aa = -(a + Tp(m)) * (apb + Tp(m)) * x
	       / ((a + Tp(m2)) * (ap1 + Tp(m2)));
	  d = Tp{1} + aa * d;
	  if (std::abs(d) < s_fpmin)
	    d = s_fpmin;
	  c = Tp{1} + aa / c;
	  if (std::abs(c) < s_fpmin)
	    c = s_fpmin;
	  d = Tp{1} / d;
	  auto del = d * c;
	  h *= del;

	  if (std::abs(del - Tp{1}) < s_eps)
	    return h;
	}
      throw std::runtime_error("ibeta_cont_frac: continued fractions failed to converge");
    }

  /**
   * Return the regularized incomplete beta function, @f$ I_x(a,b) @f$,
   * of arguments @c a, @c b, and @c x.
   *
   * The regularized incomplete beta function is defined by:
   * @f[
   *   I_x(a,b) = \frac{B_x(a,b)}{B(a,b)}
   * @f]
   * where
   * @f[
   *   B_x(a,b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt
   * @f]
   * is the non-regularized beta function and @f$ B(a,b) @f$
   * is the usual beta function.
   *
   * @param a The first parameter
   * @param b The second parameter
   * @param x The argument
   */
  template<typename Tp>
    Tp
    beta_inc(Tp a, Tp b, Tp x)
    {
      const auto s_NaN = emsr::quiet_NaN(x);


      if (x < Tp{0} || x > Tp{1})
	throw std::domain_error("beta_inc: argument out of range");
      else if (std::isnan(x) || std::isnan(a) || std::isnan(b))
	return s_NaN;
      else if (a == Tp{0} && b == Tp{0})
	return s_NaN;
      else if (a == Tp{0})
	{
	  if (x > Tp{0})
	    return Tp{1};
	  else
	    return Tp{0};
	}
      else if (b == Tp{0})
	{
	  if (x < Tp{1})
	    return Tp{0};
	  else
	    return Tp{1};
	}
      else
	{
	  auto sign = log_gamma_sign(a + b)
		      * log_gamma_sign(a) * log_gamma_sign(b);
	  auto fact = sign * std::exp(log_gamma(a + b)
		      - log_gamma(a) - log_gamma(b)
		      + a * std::log(x) + b * std::log(Tp{1} - x));

	  if (x < (a + Tp{1}) / (a + b + Tp{2}))
	    return fact * ibeta_cont_frac(a, b, x) / a;
	  else
	    return Tp{1}
		 - fact * ibeta_cont_frac(b, a, Tp{1} - x) / b;
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_BETA_TCC
