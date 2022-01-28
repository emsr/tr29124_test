
// Copyright (C) 2017-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_bernoulli.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */


#ifndef SF_BERNOULLI_TCC
#define SF_BERNOULLI_TCC 1

#include <emsr/numeric_limits.h>

namespace emsr
{
namespace detail
{
  /**
   * @brief This returns Bernoulli numbers from a table or by summation
   * for larger values.
   * @f[
   *    B_{2n} = (-1)^{n+1} 2\frac{(2n)!}{(2\pi)^{2n}} \zeta(2n)
   * @f]
   *
   * Note that
   * @f[
   *    \zeta(2n) - 1 = (-1)^{n+1} \frac{(2\pi)^{2n}}{(2n)!} B_{2n} - 2
   * @f]
   * are small and rapidly decreasing finctions of n.
   *
   * @param n the order n of the Bernoulli number.
   * @return  The Bernoulli number of order n.
   */
  template<typename Tp>
    constexpr Tp
    bernoulli_series(unsigned int n)
    {
      constexpr unsigned long s_num_bern_tab = 12;
      constexpr Tp
      s_bernoulli_2n[s_num_bern_tab]
      {
	 Tp{1ULL},
	 Tp{1ULL}             / Tp{6ULL},
	-Tp{1ULL}             / Tp{30ULL},
	 Tp{1ULL}             / Tp{42ULL},
	-Tp{1ULL}             / Tp{30ULL},
	 Tp{5ULL}             / Tp{66ULL},
	-Tp{691ULL}           / Tp{2730ULL},
	 Tp{7ULL}             / Tp{6ULL},
	-Tp{3617ULL}          / Tp{510ULL},
	 Tp{43867ULL}         / Tp{798ULL},
	-Tp{174611ULL}        / Tp{330ULL},
	 Tp{854513ULL}        / Tp{138ULL}
      };

      if (n == 0)
	return Tp{1};
      else if (n == 1)
	return -Tp{1} / Tp{2};
      // Take care of the rest of the odd ones.
      else if (n % 2 == 1)
	return Tp{0};
      // Take care of some small evens that are painful for the series.
      else if (n / 2 < s_num_bern_tab)
	return s_bernoulli_2n[n / 2];
      else
	{
	  constexpr auto s_2pi = emsr::tau_v<Tp>;
	  auto fact = Tp{1};
	  if ((n / 2) % 2 == 0)
	    fact *= -Tp{1};
	  for (unsigned int k = 1; k <= n; ++k)
	    fact *= k / s_2pi;
	  fact *= Tp{2};

	 // Riemann zeta function minus-1 for even integer argument.
	  auto sum = Tp{0};
	  for (unsigned int i = 2; i < 1000; ++i)
	    {
	      auto term = std::pow(Tp(i), -Tp(n));
	      sum += term;
	      if (term < emsr::epsilon<Tp>() * sum)
		break;
	    }

	  return fact + fact * sum;
	}
    }

  /**
   *   @brief This returns Bernoulli number @f$ B_n @f$.
   *
   *   @param n the order n of the Bernoulli number.
   *   @return  The Bernoulli number of order n.
   */
  template<typename Tp>
    constexpr Tp
    bernoulli(unsigned int n)
    { return bernoulli_series<Tp>(n); }

  /**
   * @brief This returns Bernoulli number @f$ B_2n @f$ at even integer
   * arguments @f$ 2n @f$.
   *
   * @param n the half-order n of the Bernoulli number.
   * @return  The Bernoulli number of order 2n.
   */
  template<typename Tp>
    constexpr Tp
    bernoulli_2n(unsigned int n)
    { return bernoulli_series<Tp>(2 * n); }

  /**
   * Return the Bernoulli polynomial @f$ B_n(x) @f$ of order n at argument x.
   *
   * The values at 0 and 1 are equal to the corresponding Bernoulli number:
   * @f[
   *   B_n(0) = B_n(1) = B_n
   * @f]
   *
   * The derivative is proportional to the previous polynomial:
   * @f[
   *   B_n'(x) = n * B_{n-1}(x)
   * @f]
   *
   * The series expansion is:
   * @f[
   *   B_n(x) = \sum_{k=0}^{n} B_k binom{n}{k} x^{n-k}
   * @f]
   *
   * A useful argument promotion is:
   * @f[
   *   B_n(x+1) - B_n(x) = n * x^{n-1}
   * @f]
   */
  template<typename Tp>
    Tp
    bernoulli(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto _B_n = bernoulli<Tp>(0);
	  auto binomial = Tp{1};
	  for (auto k = 1u; k <= n; ++k)
	    {
	      binomial *= Tp(n + 1 - k) / Tp(k);
	      _B_n = x * _B_n + binomial
		   * bernoulli<Tp>(k);
	    }
	  return _B_n;
	}
    }
} // namespace detail
} // namespace emsr

#endif // SF_BERNOULLI_TCC
