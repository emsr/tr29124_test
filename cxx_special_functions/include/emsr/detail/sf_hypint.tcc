
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_hypint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_HYPINT_TCC
#define SF_HYPINT_TCC 1

#include <stdexcept>
#include <complex>

namespace emsr
{
namespace detail
{

  /**
   *  @brief This function computes the hyperbolic cosine @f$ Chi(x) @f$
   *    and hyperbolic sine @f$ Shi(x) @f$ integrals
   *    by continued fraction for positive argument.
   */
  ////FIXME!!!!
  template<typename Tp>
    void
    chshint_cont_frac(Tp t, Tp& Chi, Tp& Shi)
    {
      const unsigned int s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};

      // Evaluate Chi and Shi by Lentz's modified method of continued fracions.
      std::complex<Tp> b(Tp{1}, t);
      std::complex<Tp> c(Tp{1} / s_fp_min);
      std::complex<Tp> d(Tp{1} / b);
      std::complex<Tp> h(d);
      unsigned int i = 2;
      while (true)
	{
	  auto a = -Tp(i - 1) * Tp(i - 1);
	  b += Tp{2};
	  d = Tp{1} / (a * d + b);
	  c = b + a / c;
	  auto del = c * d;
	  h *= del;
	  if (std::abs(del.real() - Tp{1}) + std::abs(del.imag()) < s_eps)
	    break;
	  if (i > s_max_iter)
	    throw std::runtime_error("chshint_cont_frac: Continued fraction evaluation failed");
	  ++i;
	}
      h *= std::polar(Tp{1}, -t);
      Chi = -h.real();
      Shi = s_pi_2 + h.imag();

      return;
    }


  /**
   *  @brief This function computes the hyperbolic cosine @f$ Chi(x) @f$
   *    and hyperbolic sine @f$ Shi(x) @f$ integrals
   *    by series summation for positive argument.
   */
  template<typename Tp>
    void
    chshint_series(Tp t, Tp& Chi, Tp& Shi)
    {
      const auto s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_gamma_e = emsr::egamma_v<Tp>;

      // Evaluate Chi and Shi by series simultaneously.
      Tp Csum(0), Ssum(0);
      if (t * t < s_fp_min)
	{
	  // Avoid underflow.
	  Csum = Tp{0};
	  Ssum = t;
	}
      else
	{
	  // Evaluate Shi and Chi by series expansion.
	  Tp sum(0);
	  Tp fact(1);
	  auto odd = true;
	  auto k = 1;
	  while (true)
	    {
	      fact *= t / k;
	      Tp term = fact / k;
	      sum += term;
	      Tp err = term / std::abs(sum);
	      if (odd)
		{
		  Ssum = sum;
		  sum = Csum;
		}
	      else
		{
		  Csum = sum;
		  sum = Ssum;
		}
	      if (err < s_eps)
		break;
	      odd = !odd;
	      ++k;
	      if (k > s_max_iter)
		throw std::runtime_error("chshint_series: Series evaluation failed");
	    }
	}
      Chi = s_gamma_e + std::log(t) + Csum;
      Shi = Ssum;

      return;
    }


  /**
   *  @brief This function returns the hyperbolic cosine @f$ Ci(x) @f$
   *    and hyperbolic sine @f$ Si(x) @f$ integrals as a pair.
   *
   *  The hyperbolic cosine integral is defined by:
   *  @f[
   *      Chi(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cosh(t) - 1}{t}
   *  @f]
   *
   *  The hyperbolic sine integral is defined by:
   *  @f[
   *      Shi(x) = \int_0^x dt \frac{\sinh(t)}{t}
   *  @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    chshint(Tp x, Tp& Chi, Tp& Shi)
    {
      const auto s_NaN = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return std::make_pair(s_NaN, s_NaN);

      auto t = std::abs(x);
      if (t == Tp{0})
	{
	  Chi = -emsr::infinity(x);
	  Shi = Tp{0};
	}
      else if (t > Tp{2})
	chshint_cont_frac(t, Chi, Shi);
      else
	chshint_series(t, Chi, Shi);

      if (x < Tp{0})
	Shi = -Shi;

      return std::make_pair(Chi, Shi);
    }

} // namespace detail
} // namespace emsr

#endif // SF_HYPINT_TCC
