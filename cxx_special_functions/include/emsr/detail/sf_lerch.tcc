
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

/** @file bits/sf_hyperg.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) 

#ifndef SF_LERCH_TCC
#define SF_LERCH_TCC 1

#include <emsr/summation.h>

namespace emsr
{
namespace detail
{

  /**
   * A functor for a vanWijnGaarden compressor.
   * vanWijnGaarden requires:
   *   Tp operator()(int) that returns a term in the original defining series.
   */
  template<typename Tp>
    class lerch_term
    {
    public:

      using value_type = Tp;

      lerch_term(value_type z, value_type s, value_type a)
      : m_z{z}, m_s{s}, m_a{a}
      { }

      value_type
      operator()(std::size_t i) const
      {
	return std::pow(m_z, value_type(i))
	     / std::pow(m_a + value_type(i), m_s);
      }

    private:

      value_type m_z;
      value_type m_s;
      value_type m_a;
    };

  /**
   * This function blows up on nonpositive integeral parameter a.
   *
   * @param z The argument.
   * @param s The order @f$ s != 1 @f$.
   * @param a The scale parameter @f$ a > -1 @f$.
   */
  template<typename Tp>
    Tp
    lerch_sum(Tp z, Tp s, Tp a)
    {
      const auto s_nan = emsr::quiet_NaN(s);
      const auto s_eps = emsr::epsilon(s);

      const auto aint = emsr::fp_is_integer(a);
      if (aint && aint() <= 0)
	return s_nan;
      else if (std::abs(std::abs(z) - Tp{1}) < s_eps
		&& std::real(s) <= Tp{1} + s_eps)
	return s_nan;
      else if (std::abs(z) > Tp{1} + s_eps)
	return s_nan;
      else
	{
	  constexpr auto s_maxit = 100000u;
	  auto zpow = Tp{1};
	  auto sum = std::pow(a, -s);
	  for (auto k = 1u; k < s_maxit; ++k)
	    {
	      zpow *= z;
	      auto term = zpow * std::pow(a + k, -s);
	      sum += term;
	      if (std::abs(term / sum) < s_eps)
		break;
	    }
	  return sum;
	}
    }

  /**
   * Try the WenigerDelta<MonotoneVanWijngaarden> composition.
   *
   * @param z The argument.
   * @param s The order @f$ s != 1 @f$.
   * @param a The scale parameter @f$ a > -1 @f$.
   */
  template<typename Tp>
    Tp
    lerch_delta_vanwijngaarden_sum(Tp z, Tp s, Tp a)
    {
      const auto s_eps = emsr::epsilon(s);
      constexpr auto s_maxit = 1000u;

      emsr::WenigerDeltaSum<emsr::VanWijngaardenSum<Tp>> WDvW;
      if (z >= Tp{0})
	{
	  using lerch_t = lerch_term<Tp>;
	  using lerch_cmp_t = emsr::VanWijngaardenCompressor<lerch_t>;
	  auto VwT = lerch_cmp_t(lerch_t(z, s, a));
	  for (auto k = 0u; k < s_maxit; ++k)
	    {
	      auto term = VwT[k];
	      WDvW += term;
	      if (std::abs(term) < s_eps * std::abs(WDvW()))
		break;
	    }
	  return WDvW();
	}
      else
	{
	  auto LT = lerch_term<Tp>(z, s, a);
	  for (auto k = 0u; k < s_maxit; ++k)
	    {
	      auto term = LT(k);
	      WDvW += term;
	      if (std::abs(term) < s_eps * std::abs(WDvW()))
		break;
	    }
	  return WDvW();
	}
    }

  /**
   * Return the Lerch transcendent @f$ \Phi(z,s,a) @f$.
   *
   * The series is:
   * @f[   *
   *   \Phi(z,s,a) = \sum_{k=0}^{\infty}\frac{z^k}{(a+k^s}
   * @f]
   *
   * This function blows up on nonpositive integeral parameter a.
   *
   * @param z The argument.
   * @param s The order @f$ s != 1 @f$.
   * @param a The scale parameter @f$ a > -1 @f$.
   */
  template<typename Tp>
    Tp
    lerch_phi(Tp z, Tp s, Tp a)
    {
      const auto s_nan = emsr::quiet_NaN(s);
      const auto s_eps = emsr::epsilon(s);

      if (std::isnan(z) || std::isnan(s) || std::isnan(a))
	return s_nan;
      else if (std::abs(std::abs(z) - Tp{1}) < s_eps
		&& std::real(s) <= Tp{1} + s_eps)
	return s_nan;
      else if (std::abs(z) > Tp{1} + s_eps)
	return s_nan;
      else
	{
	  const auto aint = emsr::fp_is_integer(a);

	  const auto sint = emsr::fp_is_integer(s);
	  const bool tinyz = std::abs(z) < s_eps; // s_min?
	  const bool smallz = !tinyz && (std::abs(z) < Tp{0.5});

	  if (aint && aint() <= 0)
	    return s_nan;
	  else if (a < Tp{0})
	    {
	      if (sint)
		{
		  int sign = sint() % 2 == 0 ? +1 : -1;
		  if (tinyz)
		    return sign * Tp{1} / std::pow(std::abs(a), s);
		  else
		    {
		      const auto m = -int(std::floor(a));
		      const auto a1 = a + Tp(m);
		      auto sum1 = Tp{0};
		      for (int i = 0; i < m; ++i)
			{
			  sum1 += sign * std::pow(std::abs(z), i)
				 / std::pow(std::abs(a + i), Tp(sint()));
			  if (z < Tp{0})
			    sign = -sign;
			}
		      auto sum = Tp{0};
		      if (smallz)
			sum = lerch_sum(z, s, a1);
		      else
			sum
			  = lerch_delta_vanwijngaarden_sum(z, s, a1);
		      sign = 1;
		      if (z < Tp{0} && m % 2 != 0)
			sign = -1;
		      return sum1
			   + sum * sign * std::pow(std::abs(z), m);
		    }
		}
	      else // s is not an integer - Phi is complex.
		return s_nan;
	    }
	  else if (tinyz)
	    return Tp{1} / std::pow(a, s);
	  else // a > 0
	    {
	      if (smallz)
		return lerch_sum(z, s, a);
	      else
		return lerch_delta_vanwijngaarden_sum(z, s, a);
	    }
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_LERCH_TCC
