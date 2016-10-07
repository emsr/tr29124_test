// math special functions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
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

/** @file bits/summation.tcc
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SUMMATION_TCC
#define _GLIBCXX_BITS_SUMMATION_TCC 1

#pragma GCC system_header

#include <vector>
#include <utility> // For exchange
#include <bits/complex_util.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Add a new term to a vanWijnGaarden sum.
   */
  template<typename _Tp>
    _VanWijngaardenSum<_Tp>&
    _VanWijngaardenSum<_Tp>::operator+=(value_type __term)
    {
      if (std::__detail::__isnan(__term))
	std::__throw_runtime_error(__N("_VanWijngaardenSum: bad term"));
      if (std::__detail::__isinf(__term))
	std::__throw_runtime_error(__N("_VanWijngaardenSum: infinite term"));

      ++this->_M_num_terms;
      this->_M_term = __term;

      if (this->_M_delta.size() == 0)
	{
	  this->_M_delta.push_back(__term);
	  this->_M_sum += value_type{0.5L} * this->_M_delta.back();
	}
      else
	{
	  auto __temp = this->_M_delta[0];
	  this->_M_delta[0] = __term;
	  auto __n = this->_M_delta.size();
	  for (auto __j = 0; __j < __n - 1; ++__j)
	    __temp = std::exchange(this->_M_delta[__j + 1],
			     value_type{0.5L} * (this->_M_delta[__j] + __temp));
	  auto __next = value_type{0.5L} * (this->_M_delta.back() + __temp);
	  if (std::abs(__next) < std::abs(this->_M_delta.back()))
	    {
	      this->_M_delta.push_back(__next);
	      this->_M_sum += value_type{0.5L} * this->_M_delta.back();
	    }
	  else
	    this->_M_sum += __next;
	}
/*
      lasteps = std::abs(sum - lastval);
      if (lasteps <= eps)
	++_M_num_convergences;
      if (_M_num_convergences >= 2)
	this->_M_converged = true;
//return (lastval = sum);
*/
      return *this;
    }

  /**
   *  Perform a series compression on a monotone series - converting
   *  it to an alternating series - for the regular van Wijngaarden sum.
   *  ADL for ctors anyone?  I'd like to put a lambda in here*
   */
  template<typename _TermFn>
    auto
    _VanWijngaardenCompressor<_TermFn>::operator[](std::size_t __j) const
    {
      using value_type = decltype(this->_M_term_fn(__j));
      constexpr auto _S_min = std::numeric_limits<value_type>::min();
      constexpr auto _S_eps = std::numeric_limits<value_type>::epsilon();
      // Maximum number of iterations before 2^k overflow.
      constexpr auto __k_max = std::numeric_limits<std::size_t>::digits;

      auto __sum = value_type{};
      auto __two2k = std::size_t{1};
      for (auto __k = std::size_t{0}; __k < __k_max; __k += std::size_t{2})
	{
	  // Index for the term in the original series.
	  auto __i = std::size_t{0};
	  if (__builtin_mul_overflow(__two2k, __j + 1, &__i))
	    std::__throw_runtime_error(__N("_VanWijngaardenCompressor: "
					   "index overflow"));	  
	  --__i;

	  // Increment the sum.
	  auto __term = __two2k * this->_M_term_fn(__i);
	  __sum += __term;

	  // Stop summation if either the sum is zero
	  // or if |term / sum| is below requested accuracy.
	  if (std::abs(__sum) <= _S_min
	   || std::abs(__term / __sum) < value_type{1.0e-2} * _S_eps)
	    break;

	  if (__builtin_mul_overflow(__two2k, std::size_t{2}, &__two2k))
	    std::__throw_runtime_error(__N("_VanWijngaardenCompressor: "
					   "index overflow"));
	}

      auto __sign = (__j % 2 == 1 ? -1 : +1);
      return __sign * __sum;
    }

  /**
   * Perform one step of the Aitken delta-squared process.
   */
  template<typename _Sum>
    void
    _AitkenDeltaSquaredSum<_Sum>::_M_update()
    {
      using _Tp = value_type;
      using _Val = std::__detail::__num_traits_t<_Tp>;
      const/*expr*/ auto _S_huge = __gnu_cxx::__root_max(_Val{5}); // 1.0e+60
      const/*expr*/ auto _S_tiny = __gnu_cxx::__root_min(_Val{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      auto& __a = this->_M_a;

      __a.push_back(__s_n);
      if (__n < 2)
	this->_M_sum = __s_n;
      else
	{
	  auto __lowmax = __n / 2;
	  for (auto __j = 1; __j <= __lowmax; ++__j)
	    {
	      auto __m = __n - 2 * __j;
	      auto __denom = (__a[__m + 2] - __a[__m + 1])
			   - (__a[__m + 1] - __a[__m]);
	      if (std::abs(__denom) < _S_tiny)
		__a[__m] = _S_huge;
	      else
		{
		  auto __del = __a[__m] - __a[__m + 1];
		  __a[__m] -= __del * __del / __denom;
		}
	    }
	  this->_M_sum = __a[__n % 2];
	}
    }

  /**
   * Perform one step of the Winn epsilon transformation.
   */
  template<typename _Sum>
    void
    _WinnEpsilonSum<_Sum>::_M_update()
    {
      using _Tp = value_type;
      using _Val = std::__detail::__num_traits_t<_Tp>;
      const auto _S_huge = __gnu_cxx::__root_max(_Val{5}); // 1.0e+60
      const auto _S_tiny = __gnu_cxx::__root_min(_Val{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      auto& __e = this->_M_e;

      __e.push_back(__s_n);
      if (__n == 0)
        this->_M_sum = __s_n;
      else
        {
          auto __aux2 = _Tp{0};
          for (auto __j = __n; __j >= 1; --__j)
            {
	      auto __aux1 = __aux2;
	      __aux2 = __e[__j - 1];
	      auto __diff = __e[__j] - __aux2;
	      if (std::abs(__diff) < _S_tiny)
		__e[__j - 1] = _S_huge;
	      else
		__e[__j - 1] = __aux1 + _Tp{1} / __diff;
	    }
	  this->_M_sum = __e[__n % 2];
        }
      return;
    }

  /**
   * Perform one step of the Brezinski Theta transformation.
   */
  template<typename _Sum>
    void
    _BrezinskiThetaSum<_Sum>::_M_update()
    {
      using _Tp = value_type;
      using _Val = std::__detail::__num_traits_t<_Tp>;
      const auto _S_huge = __gnu_cxx::__root_max(_Val{5}); // 1.0e+60
      const auto _S_tiny = __gnu_cxx::__root_min(_Val{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      auto& __arj = this->_M_arj;

      __arj.push_back(__s_n);
      if (__n < 3)
	this->_M_sum = __s_n;
      else
	{
	  auto __lmax = __n / 3;
	  auto __m = __n;
	  for (auto __l = 1; __l <= __lmax; ++__l)
	    {
	      __m -= 3;
	      auto __diff0 = __arj[__m + 1] - __arj[__m];
	      auto __diff1 = __arj[__m + 2] - __arj[__m + 1];
	      auto __diff2 = __arj[__m + 3] - __arj[__m + 2];
	      auto __denom = __diff2 * (__diff1 - __diff0)
			   - __diff0 * (__diff2 - __diff1);
	      if (std::abs(__denom) < _S_tiny)
		__arj[__m] = _S_huge;
	      else
		__arj[__m] = __arj[__m + 1]
			   - __diff0 * __diff1 * (__diff2 - __diff1) / __denom;
	    }
	  this->_M_sum = __arj[__n % 3];
	}
    }

  /**
   * Perform one step of the Levin summation process.
   */
  template<typename _Sum, typename _RemainderModel>
    void
    _LevinSum<_Sum, _RemainderModel>::_M_update(value_type __r_n)
    {
      using _Tp = value_type;
      using _Val = std::__detail::__num_traits_t<_Tp>;
      const auto _S_huge = __gnu_cxx::__root_max(_Val{5}); // 1.0e+60
      const auto _S_tiny = __gnu_cxx::__root_min(_Val{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      const auto __beta = this->_M_beta;
      auto& __anum = this->_M_num;
      auto& __aden = this->_M_den;

      __anum.push_back(__s_n / __r_n);
      __aden.push_back(value_type{1} / __r_n);
      if (__n == 0)
	this->_M_sum = __s_n;
      else
	{
	  __anum[__n - 1] = __anum[__n] - __anum[__n - 1];
	  __aden[__n - 1] = __aden[__n] - __aden[__n - 1];
	  if (__n > 1)
	    {
	      auto __bn1 = __beta + _Tp(__n - 1);
	      auto __bn2 = __beta + _Tp(__n);
	      auto __coef = __bn1 / __bn2;
	      auto __coefp = _Tp{1};
	      for (auto __j = 2; __j <= __n; ++__j)
		{
		  auto __fact = (__beta + _Tp(__n - __j)) * __coefp / __bn2;
		  __anum[__n - __j] = __anum[__n - __j + 1]
				    - __fact * __anum[__n - __j];
		  __aden[__n - __j] = __aden[__n - __j + 1]
				    - __fact * __aden[__n - __j];
		  __coefp *= __coef;
		}
	    }
	  if (std::abs(__aden[0]) < _S_tiny)
	    this->_M_sum = _S_huge;
	  else
	    this->_M_sum = __anum[0] / __aden[0];
	}
    }

  /**
   * Perform one step of the Weniger summation process.
   */
  template<typename _Sum, typename _RemainderModel>
    void
    _WenigerSum<_Sum, _RemainderModel>::_M_update(value_type __r_n)
    {
      using _Tp = value_type;
      using _Val = std::__detail::__num_traits_t<_Tp>;
      const/*expr*/ auto _S_huge = __gnu_cxx::__root_max(_Val{5}); // 1.0e+60
      const/*expr*/ auto _S_tiny = __gnu_cxx::__root_min(_Val{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      const auto __beta = this->_M_beta;
      auto& __anum = this->_M_num;
      auto& __aden = this->_M_den;

      __anum.push_back(__s_n / __r_n);
      __aden.push_back(value_type{1} / __r_n);
      if (__n == 0)
	this->_M_sum = __s_n;
      else
	{
	  __anum[__n - 1] = __anum[__n] - __anum[__n - 1];
	  __aden[__n - 1] = __aden[__n] - __aden[__n - 1];
	  if (__n > 1)
	    {
	      auto __bn1 = __beta + _Tp(__n - 2);
	      auto __bn2 = __beta + _Tp(__n - 1);
	      for (auto __j = 2; __j <= __n; ++__j)
		{
		  auto __fact = __bn1 * __bn2
			   / ((__bn1 + _Tp(__j - 1)) * (__bn2 + _Tp(__j - 1)));
		  __anum[__n - __j] = __anum[__n - __j + 1]
				    - __fact * __anum[__n - __j];
		  __aden[__n - __j] = __aden[__n - __j + 1]
				    - __fact * __aden[__n - __j];
		}
	    }
	  if (std::abs(__aden[0]) < _S_tiny)
	    this->_M_sum = _S_huge;
	  else
	    this->_M_sum = __anum[0] / __aden[0];
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SUMMATION_TCC
