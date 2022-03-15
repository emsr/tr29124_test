
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

/** @file bits/sf_expint.tcc
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
//     Ed. by Milton Abramowitz and Irene A. Stegun,
//     Dover Publications, New-York, Section 5, pp. 228-251.
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 222-225.
//

#ifndef SF_EXPINT_TCC
#define SF_EXPINT_TCC 1

#include <stdexcept>

#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{

  template<typename Tp> Tp expint_E1(Tp);

  /**
   * @brief Return the exponential integral @f$ E_1(x) @f$ by series summation.
   * 	    This should be good for @f$ x < 1 @f$.
   *
   * The exponential integral is given by
   *  @f[
   *    E_1(x) = \int_{1}^{\infty} \frac{e^{-xt}}{t} dt
   *  @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_E1_series(Tp x)
    {
      const auto s_eps = emsr::epsilon(x);
      auto term = Tp{1};
      auto esum = Tp{0};
      auto osum = Tp{0};
      const unsigned int max_iter = 1000;
      for (unsigned int i = 1; i < max_iter; ++i)
	{
	  term *= - x / i;
	  if (std::abs(term)
		 < s_eps * std::min(std::abs(esum), std::abs(osum)))
	    break;
	  if (term >= Tp{0})
	    esum += term / i;
	  else
	    osum += term / i;
	}

      return - esum - osum
	     - emsr::egamma_v<Tp> - std::log(x);
    }


  /**
   * @brief Return the exponential integral @f$ E_1(x) @f$
   * 	    by asymptotic expansion.
   *
   * The exponential integral is given by
   * @f[
   *   E_1(x) = \int_{1}^\infty \frac{e^{-xt}}{t} dt
   * @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_E1_asymp(Tp x)
    {
      auto term = Tp{1};
      auto esum = Tp{1};
      auto osum = Tp{0};
      const unsigned int max_iter = 1000;
      for (unsigned int i = 1; i < max_iter; ++i)
	{
	  auto prev = term;
	  term *= - i / x;
	  if (std::abs(term) > std::abs(prev))
	    break;
	  if (term >= Tp{0})
	    esum += term;
	  else
	    osum += term;
	}

      return std::exp(-x) * (esum + osum) / x;
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$ by series summation.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  n  The order of the exponential integral function.
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_En_series(unsigned int n, Tp x)
    {
      const unsigned int s_max_iter = 1000;
      const auto s_eps = emsr::epsilon(x);
      const int nm1 = n - 1;
      const auto s_gamma_E = emsr::egamma_v<Tp>;
      const auto logx = std::log(x);
      Tp sum = (nm1 != 0
		? Tp{1} / nm1
		: -logx - s_gamma_E);
      Tp fact = Tp{1};
      for (unsigned int i = 1; i <= s_max_iter; ++i)
	{
	  fact *= -x / Tp(i);
	  Tp term;
	  if (int(i) != nm1)
	    term = -fact / Tp(i - nm1);
	  else
	    {
	      Tp psi = -s_gamma_E;
	      for (int ii = 1; ii <= nm1; ++ii)
		psi += Tp{1} / Tp(ii);
	      term = fact * (psi - logx);
	    }
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    return sum;
	}
      throw std::runtime_error("expint_En_series: series summation failed");
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$
   * 	    by continued fractions.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  n  The order of the exponential integral function.
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_En_cont_frac(unsigned int n, Tp x)
    {
      const unsigned int s_max_iter = 1000;
      const auto s_eps = emsr::epsilon(x);
      const auto s_fp_min = Tp{4} * emsr::lim_min(x);
      const int nm1 = n - 1;
      auto b = x + Tp(n);
      auto c = Tp{1} / s_fp_min;
      auto d = Tp{1} / b;
      auto h = d;
      for ( unsigned int i = 1; i <= s_max_iter; ++i )
	{
	  auto a = -Tp(i * (nm1 + i));
	  b += Tp{2};
	  d = Tp{1} / (a * d + b);
	  if (std::abs(d) < s_fp_min)
	    d = std::copysign(s_fp_min, d);
	  c = b + a / c;
	  if (std::abs(c) < s_fp_min)
	    c = std::copysign(s_fp_min, c);
	  const auto del = c * d;
	  h *= del;
	  if (std::abs(del - Tp{1}) < s_eps)
	    return h * std::exp(-x);
	}
      throw std::runtime_error("expint_En_cont_frac: continued fraction failed");
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$ by recursion.
   * 	    Use upward recursion for @f$ x < n @f$
   * 	    and downward recursion (Miller's algorithm) otherwise.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  n  The order of the exponential integral function.
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_En_recursion(unsigned int n, Tp x)
    {
      Tp En;
      Tp E1 = expint_E1(x);
      if (x < Tp(n))
	{
	  // Forward recursion is stable only for n < x.
	  En = E1;
	  for (unsigned int j = 2; j < n; ++j)
	    En = (std::exp(-x) - x * En) / Tp(j - 1);
	}
      else
	{
	  // Backward recursion is stable only for n >= x.
	  En = Tp{1};
	  /// @todo Find a principled starting number
	  /// for the @f$ E_n(x) @f$ downward recursion.
	  const int N = n + 20;
	  Tp save = Tp{0};
	  for (int j = N; j > 0; --j)
	    {
	      En = (std::exp(-x) - j * En) / x;
	      if (j == n)
		save = En;
	    }
	    Tp norm = En / E1;
	    En /= norm;
	}

      return En;
    }

  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$ by series summation.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_Ei_series(Tp x)
    {
      Tp term = Tp{1};
      Tp sum = Tp{0};
      const auto s_eps = emsr::epsilon(x);
      const unsigned int max_iter = 1000;
      for (unsigned int i = 1; i < max_iter; ++i)
	{
	  term *= x / i;
	  sum += term / i;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}

      return emsr::egamma_v<Tp>
	   + sum + std::log(x);
    }


  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$
   * 	    by asymptotic expansion.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_Ei_asymp(Tp x)
    {
      Tp term = Tp{1};
      Tp sum = Tp{1};
      const auto s_eps = emsr::epsilon(x);
      const unsigned int max_iter = 1000;
      for (unsigned int i = 1; i < max_iter; ++i)
	{
	  Tp prev = term;
	  term *= i / x;
	  if (std::abs(term) >= std::abs(prev))
	    break;
	  sum += term;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	}

      return std::exp(x) * sum / x;
    }


  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_Ei(Tp x)
    {
      const auto s_eps = emsr::epsilon(x);
      if (x < Tp{0})
	return -expint_E1(-x);
      else if (x < -std::log(s_eps))
	return expint_Ei_series(x);
      else
	return expint_Ei_asymp(x);
    }


  /**
   * @brief Return the exponential integral @f$ E_1(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   E_1(x) = \int_{1}^\infty \frac{e^{-xt}}{t} dt
   * @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_E1(Tp x)
    {
      if (x < Tp{0})
	return -expint_Ei(-x);
      else if (x < Tp{1})
	return expint_E1_series(x);
      else if (x < Tp{100})
	/// @todo Find a good asymptotic switch point in @f$ E_1(x) @f$.
	return expint_En_cont_frac(1, x);
      else
	return expint_E1_asymp(x);
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$
   * 	    for large argument.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  n  The order of the exponential integral function.
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_En_asymp(unsigned int n, Tp x)
    {
      auto term = Tp{1};
      auto sum = Tp{1};
      for (unsigned int i = 1; i <= n; ++i)
	{
	  auto prev = term;
	  term *= -Tp(n - i + 1) / x;
	  if (std::abs(term) > std::abs(prev))
	    break;
	  sum += term;
	}

      return std::exp(-x) * sum / x;
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$
   * 	    for large order.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  n  The order of the exponential integral function.
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint_En_large_n(unsigned int n, Tp x)
    {
      const auto xpn = x + n;
      const auto xpn2 = xpn * xpn;
      const auto s_eps = emsr::epsilon(x);
      auto term = Tp{1};
      auto sum = Tp{1};
      for (unsigned int i = 1; i <= n; ++i)
	{
	  term *= (n - 2 * (i - 1) * x) / xpn2;
	  if (std::abs(term) < s_eps * std::abs(sum))
	    break;
	  sum += term;
	}

      return std::exp(-x) * sum / xpn;
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  n  The order of the exponential integral function.
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (n <= 1 && x == Tp{0})
	return emsr::infinity(x);
      else
	{
	  if (n == 0)
	    return std::exp(-x) / x;
	  else if (n == 1)
	    return expint_E1(x);
	  else if (x == Tp{0})
	    return Tp{1} / static_cast<Tp>(n - 1);
	  else if (x < Tp{1})
	    return expint_En_series(n, x);
	  else if (n > 50000)
	    /// @todo Study arbitrary switch to large-n @f$ E_n(x) @f$.
	    return expint_En_large_n(n, x);
	  else if (x > Tp{100})
	    /// @todo Find a good asymptotic switch point in @f$ E_n(x) @f$.
	    return expint_En_asymp(n, x);
	  else
	    return expint_En_cont_frac(n, x);
	}
    }


  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename Tp>
    Tp
    expint(Tp x)
    {
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return expint_Ei(x);
    }

  /**
   * @brief Return the logarithmic integral @f$ li(x) @f$.
   *
   * The logarithmic integral is given by
   * @f[
   *   li(x) = Ei(\log(x))
   * @f]
   *
   * @param  x  The argument of the logarithmic integral function.
   * @return  The logarithmic integral.
   */
  template<typename Tp>
    Tp
    logint(const Tp x)
    {
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (std::abs(x) == Tp{1})
	return emsr::infinity(x);
      else
	return expint(std::log(x));
    }

  /**
   * @brief Return the hyperbolic cosine integral @f$ Chi(x) @f$.
   *
   * The hyperbolic cosine integral is given by
   * @f[
   *   Chi(x) = (Ei(x) - E_1(x))/ 2 = (Ei(x) + Ei(-x))/2
   * @f]
   *
   * @param  x  The argument of the hyperbolic cosine integral function.
   * @return  The hyperbolic cosine integral.
   */
  template<typename Tp>
    Tp
    coshint(const Tp x)
    {
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == Tp{0})
	return Tp{0};
      else
	return (expint_Ei(x) - expint_E1(x)) / Tp{2};
    }

  /**
   * @brief Return the hyperbolic sine integral @f$ Shi(x) @f$.
   *
   * The hyperbolic sine integral is given by
   * @f[
   *   Shi(x) = (Ei(x) + E_1(x))/2 = (Ei(x) - Ei(-x))/2
   * @f]
   *
   * @param  x  The argument of the hyperbolic sine integral function.
   * @return  The hyperbolic sine integral.
   */
  template<typename Tp>
    Tp
    sinhint(const Tp x)
    {
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return (expint_Ei(x) + expint_E1(x)) / Tp{2};
    }

} // namespace detail
} // namespace emsr

#endif // SF_EXPINT_TCC
