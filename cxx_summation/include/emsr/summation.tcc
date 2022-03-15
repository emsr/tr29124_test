
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

/** @file bits/summation.tcc
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SUMMATION_TCC
#define SUMMATION_TCC 1

#include <vector>
#include <emsr/numeric_limits.h>
#include <utility> // For exchange.

namespace emsr
{

  /**
   * Add a new term to a vanWijnGaarden sum.
   */
  template<typename Tp>
    VanWijngaardenSum<Tp>&
    VanWijngaardenSum<Tp>::operator+=(value_type term)
    {
      if (std::isnan(term))
	throw std::runtime_error("VanWijngaardenSum: bad term");
      if (std::isinf(term))
	throw std::runtime_error("VanWijngaardenSum: infinite term");
      if (this->m_num_terms > 1 && this->m_term * term > value_type{0})
	throw std::runtime_error("VanWijngaardenSum: terms not alternating in sign");

      ++this->m_num_terms;
      this->m_term = term;

      if (this->m_num_terms <= this->m_start_term)
	this->m_sum += term;
      else if (this->m_delta.size() == 0)
	{
	  this->m_delta.push_back(term);
	  this->m_sum += value_type{0.5L} * this->m_delta.back();
	}
      else
	{
	  auto temp = this->m_delta[0];
	  this->m_delta[0] = term;
	  auto n = this->m_delta.size();
	  for (auto j = 0ull; j < n - 1; ++j)
	    temp = std::exchange(this->m_delta[j + 1],
			     value_type{0.5L} * (this->m_delta[j] + temp));
	  auto next = value_type{0.5L} * (this->m_delta.back() + temp);
	  if (std::abs(next) < std::abs(this->m_delta.back()))
	    {
	      this->m_delta.push_back(next);
	      this->m_sum += value_type{0.5L} * this->m_delta.back();
	    }
	  else
	    this->m_sum += next;
	}
/*
      lasteps = std::abs(sum - lastval);
      if (lasteps <= eps)
	++m_num_convergences;
      if (m_num_convergences >= 2)
	this->m_converged = true;
//return (lastval = sum);
*/
      return *this;
    }

  /**
   *  Perform a series compression on a monotone series - converting
   *  it to an alternating series - for the regular van Wijngaarden sum.
   *  ADL for ctors anyone?  I'd like to put a lambda in here*
   */
  template<typename TermFn>
    typename VanWijngaardenCompressor<TermFn>::return_t
    VanWijngaardenCompressor<TermFn>::operator[](std::size_t j) const
    {
      using value_type = decltype(this->m_term_fn(j));
      const auto s_min = emsr::lim_min<value_type>();
      const auto s_eps = emsr::epsilon<value_type>();
      // Maximum number of iterations before 2^k overflow.
      const auto k_max = std::numeric_limits<std::size_t>::digits;

      auto sum = value_type{};
      auto two2k = std::size_t{1};
      for (auto k = std::size_t{0}; k < k_max; k += std::size_t{2})
	{
	  // Index for the term in the original series.
	  auto i = std::size_t{0};
	  if (__builtin_mul_overflow(two2k, j + 1, &i))
	    throw std::runtime_error("VanWijngaardenCompressor: index overflow");	  
	  --i;

	  // Increment the sum.
	  auto term = two2k * this->m_term_fn(i);
	  sum += term;

	  // Stop summation if either the sum is zero
	  // or if |term / sum| is below requested accuracy.
	  if (std::abs(sum) <= s_min
	   || std::abs(term / sum) < value_type{1.0e-2} * s_eps)
	    break;

	  if (__builtin_mul_overflow(two2k, std::size_t{2}, &two2k))
	    throw std::runtime_error("VanWijngaardenCompressor: index overflow");
	}

      auto sign = (j % 2 == 1 ? -1 : +1);
      return sign * sum;
    }

  /**
   * Perform one step of the Aitken delta-squared process.
   */
  template<typename Sum>
    void
    AitkenDeltaSquaredSum<Sum>::m_update()
    {
      using Tp = value_type;
      using Val = emsr::num_traits_t<Tp>;
      const auto s_huge = emsr::root_max(Val{5}); // 1.0e+60
      const auto s_tiny = emsr::root_min(Val{5}); // 1.0e-60;

      const auto n = this->m_part_sum.num_terms() - 1;
      const auto s_n = this->m_part_sum();

      this->m_a.push_back(s_n);
      if (n < 2)
	this->m_sum = s_n;
      else
	{
	  auto lowmax = n / 2;
	  for (auto j = 1u; j <= lowmax; ++j)
	    {
	      auto m = n - 2 * j;
	      auto denom = (this->m_a[m + 2] - this->m_a[m + 1])
			   - (this->m_a[m + 1] - this->m_a[m]);
	      if (std::abs(denom) < s_tiny)
		this->m_a[m] = s_huge;
	      else
		{
		  auto del = this->m_a[m] - this->m_a[m + 1];
		  this->m_a[m] -= del * del / denom;
		}
	    }
	  this->m_sum = this->m_a[n % 2];
	}
    }

  /**
   * Perform one step of the Winn epsilon transformation.
   */
  template<typename Sum>
    void
    WinnEpsilonSum<Sum>::m_update()
    {
      using Tp = value_type;
      using Val = emsr::num_traits_t<Tp>;
      const auto s_huge = emsr::root_max(Val{5}); // 1.0e+60
      const auto s_tiny = emsr::root_min(Val{5}); // 1.0e-60;

      const auto n = this->m_part_sum.num_terms() - 1;
      const auto s_n = this->m_part_sum();

      this->m_e.push_back(s_n);
      if (n == 0)
        this->m_sum = s_n;
      else
        {
          auto aux2 = Tp{0};
          for (auto j = n; j >= 1; --j)
            {
	      auto aux1 = aux2;
	      aux2 = this->m_e[j - 1];
	      auto diff = this->m_e[j] - aux2;
	      if (std::abs(diff) < s_tiny)
		this->m_e[j - 1] = s_huge;
	      else
		this->m_e[j - 1] = aux1 + Tp{1} / diff;
	    }
	  this->m_sum = this->m_e[n % 2];
        }
      return;
    }

  /**
   * Perform one step of the Brezinski Theta transformation.
   */
  template<typename Sum>
    void
    BrezinskiThetaSum<Sum>::m_update()
    {
      using Tp = value_type;
      using Val = emsr::num_traits_t<Tp>;
      const auto s_huge = emsr::root_max(Val{5}); // 1.0e+60
      const auto s_tiny = emsr::root_min(Val{5}); // 1.0e-60;

      const auto n = this->m_part_sum.num_terms() - 1;
      const auto s_n = this->m_part_sum();

      this->m_arj.push_back(s_n);
      if (n < 3)
	this->m_sum = s_n;
      else
	{
	  auto lmax = n / 3;
	  auto m = n;
	  for (auto l = 1u; l <= lmax; ++l)
	    {
	      m -= 3;
	      auto diff0 = this->m_arj[m + 1] - this->m_arj[m];
	      auto diff1 = this->m_arj[m + 2] - this->m_arj[m + 1];
	      auto diff2 = this->m_arj[m + 3] - this->m_arj[m + 2];
	      auto denom = diff2 * (diff1 - diff0)
			   - diff0 * (diff2 - diff1);
	      if (std::abs(denom) < s_tiny)
		this->m_arj[m] = s_huge;
	      else
		this->m_arj[m] = this->m_arj[m + 1]
				  - diff0 * diff1
				    * (diff2 - diff1) / denom;
	    }
	  this->m_sum = this->m_arj[n % 3];
	}
    }

  /**
   * Perform one step of the Levin summation process.
   */
  template<typename Sum, typename RemainderModel>
    void
    LevinSum<Sum, RemainderModel>::m_update(value_type r_n)
    {
      using Tp = value_type;
      using Val = emsr::num_traits_t<Tp>;
      const auto s_huge = emsr::root_max(Val{5}); // 1.0e+60
      const auto s_tiny = emsr::root_min(Val{5}); // 1.0e-60;

      const auto n = this->m_part_sum.num_terms() - 1;
      const auto s_n = this->m_part_sum();

      this->m_num.push_back(s_n / r_n);
      this->m_den.push_back(value_type{1} / r_n);
      if (n == 0)
	this->m_sum = s_n;
      else
	{
	  this->m_num[n - 1] = this->m_num[n] - this->m_num[n - 1];
	  this->m_den[n - 1] = this->m_den[n] - this->m_den[n - 1];
	  if (n > 1)
	    {
	      auto bn1 = this->m_beta + Tp(n - 1);
	      auto bn2 = this->m_beta + Tp(n);
	      auto coef = bn1 / bn2;
	      auto coefp = Tp{1};
	      for (auto j = 2u; j <= n; ++j)
		{
		  auto fact = (this->m_beta + Tp(n - j)) * coefp / bn2;
		  this->m_num[n - j] = this->m_num[n - j + 1]
					  - fact * this->m_num[n - j];
		  this->m_den[n - j] = this->m_den[n - j + 1]
					  - fact * this->m_den[n - j];
		  coefp *= coef;
		}
	    }
	  if (std::abs(this->m_den[0]) < s_tiny)
	    this->m_sum = s_huge;
	  else
	    this->m_sum = this->m_num[0] / this->m_den[0];
	}
    }

  /**
   * Perform one step of the Weniger summation process.
   */
  template<typename Sum, typename RemainderModel>
    void
    WenigerSum<Sum, RemainderModel>::m_update(value_type r_n)
    {
      using Tp = value_type;
      using Val = emsr::num_traits_t<Tp>;
      const/*expr*/ auto s_huge = emsr::root_max(Val{5}); // 1.0e+60
      const/*expr*/ auto s_tiny = emsr::root_min(Val{5}); // 1.0e-60;

      const auto n = this->m_part_sum.num_terms() - 1;
      const auto s_n = this->m_part_sum();

      this->m_num.push_back(s_n / r_n);
      this->m_den.push_back(value_type{1} / r_n);
      if (n == 0)
	this->m_sum = s_n;
      else
	{
	  this->m_num[n - 1] = this->m_num[n] - this->m_num[n - 1];
	  this->m_den[n - 1] = this->m_den[n] - this->m_den[n - 1];
	  if (n > 1)
	    {
	      auto bn1 = this->m_beta + Tp(n - 1);
	      auto bn2 = this->m_beta + Tp(n - 2);
	      for (auto j = 2ull; j <= n; ++j)
		{
		  auto fact = (bn1 / (bn1 + Tp(j - 1)))
			      * (bn2 / (bn2 + Tp(j - 1)));
		  this->m_num[n - j] = this->m_num[n - j + 1]
					  - fact * this->m_num[n - j];
		  this->m_den[n - j] = this->m_den[n - j + 1]
					  - fact * this->m_den[n - j];
		}
	    }
	  if (std::abs(this->m_den[0]) < s_tiny)
	    this->m_sum = s_huge;
	  else
	    this->m_sum = this->m_num[0] / this->m_den[0];
	}
    }

} // namespace emsr

#endif // SUMMATION_TCC
