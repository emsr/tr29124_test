
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

/** @file bits/sf_trigint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_TRIGINT_TCC
#define SF_TRIGINT_TCC 1

#include <stdexcept>
#include <complex>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{

  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by continued fraction
   *         for positive argument.
   */
  template<typename Tp>
    void
    sincosint_cont_frac(Tp t, Tp& Si, Tp& Ci)
    {
      const auto s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};

      // Evaluate Ci and Si by Lentz's modified method of continued fractions.
      std::complex<Tp> b(Tp{1}, t);
      std::complex<Tp> c(Tp{1} / s_fp_min);
      std::complex<Tp> d(Tp{1} / b);
      std::complex<Tp> h(d);
      int i = 2;
      while (true)
	{
	  auto a = -Tp(i - 1) * Tp(i - 1);
	  b += Tp{2};
	  d = Tp{1} / (a * d + b);
	  c = b + a / c;
	  std::complex<Tp> del = c * d;
	  h *= del;
	  if (std::abs(del - Tp{1}) < s_eps)
	    break;
	  if (i > s_max_iter)
	    throw std::runtime_error("sincosint_cont_frac: continued fraction evaluation failed");
	  ++i;
	}
      h *= std::polar(Tp{1}, -t);
      Ci = -h.real();
      Si = s_pi_2 + h.imag();

      return;
    }


  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by series summation
   *         for positive argument.
   */
  template<typename Tp>
    void
    sincosint_series(Tp t, Tp& Si, Tp& Ci)
    {
      const auto s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_gamma_e = emsr::egamma_v<Tp>;

      // Evaluate Ci and Si by series simultaneously.
      Tp sumc(0), sums(0);
      if (t * t < s_fp_min)
	{
	  // Avoid underflow.
	  sumc = Tp{0};
	  sums = t;
	}
      else
	{
	  // Evaluate Si and Ci by series expansion.
	  Tp sum(0);
	  Tp sign(1), fact(1);
	  bool odd = true;
	  unsigned int k = 1;
	  while (true)
	    {
	      fact *= t / k;
	      Tp term = fact / k;
	      sum += sign * term;
	      Tp err = term / std::abs(sum);
	      if (odd)
		{
		  sign = -sign;
		  sums = sum;
		  sum = sumc;
		}
	      else
		{
		  sumc = sum;
		  sum = sums;
		}
	      if (err < s_eps)
		break;
	      odd = !odd;
	      ++k;
	      if (k > s_max_iter)
		throw std::runtime_error("sincosint_series: series evaluation failed");
	    }
	}
      Si = sums;
      Ci = s_gamma_e + std::log(t) + sumc;

      return;
    }


  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by asymptotic series summation
   *         for positive argument.
   *
   *   The asymptotic series is very good for x > 50.
   */
  template<typename Tp>
    void
    sincosint_asymp(Tp t, Tp& Si, Tp& Ci)
    {
      const auto s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(t);
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};

      auto invt = Tp{1} / t;
      auto term = Tp{1}; // 0!
      auto sume = Tp{term};
      term *= invt; // 1! / t
      auto sumo = Tp{term};
      auto sign = Tp{1};
      auto even = true;
      auto k = 2;
      while (true)
	{
	  term *= k * invt;

	  if (even)
	    {
	      sign = -sign;
	      sume += sign * term;

	    }
	  else
	    {
	      sumo += sign * term;
	      if (term / std::abs(sumo) < s_eps)
		break;
	    }

	  even = !even;

	  if (k > s_max_iter)
	    throw std::runtime_error("sincosint_asymp: Series evaluation failed");
	  ++k;
	}

      auto sint = std::sin(t);
      auto cost = std::cos(t);
      Si = s_pi_2
	  - cost * invt * sume
	  - sint * invt * sumo;
      Ci = sint * invt * sume
	  - cost * invt * sumo;

      return;
    }


  /**
   *  @brief This function returns the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals as a @c pair.
   *
   *  The sine integral is defined by:
   *  @f[
   *      Si(x) = \int_0^x dt \frac{\sin(t)}{t}
   *  @f]
   *
   *  The cosine integral is defined by:
   *  @f[
   *      Ci(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cos(t) - 1}{t}
   *  @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    sincosint(Tp x)
    {
      const auto s_NaN = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return std::make_pair(s_NaN, s_NaN);

      auto t = std::abs(x);
      Tp Ci, Si;
      if (t == Tp{0})
	{
	  Si = Tp{0};
	  Ci = -emsr::infinity(x);
	}
      else if (t > Tp{1000}) // Check this!
	sincosint_asymp(t, Si, Ci);
      else if (t > Tp{2})
	sincosint_cont_frac(t, Si, Ci);
      else
	sincosint_series(t, Si, Ci);

      if (x < Tp{0})
	Si = -Si;

      return std::make_pair(Si, Ci);
    }

} // namespace detail
} // namespace emsr

#endif // SF_TRIGINT_TCC
