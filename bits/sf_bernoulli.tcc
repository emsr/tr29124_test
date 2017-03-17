// Special functions -*- C++ -*-

// Copyright (C) 2017 Free Software Foundation, Inc.
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

/** @file bits/sf_bernoulli.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// C++ Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_SF_BERNOULLI_TCC
#define _GLIBCXX_BITS_SF_BERNOULLI_TCC 1

#pragma GCC system_header


namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief This returns Bernoulli numbers from a table or by summation
   *         for larger values.
   *
   *  Upward recursion is unstable.
   *
   *  @param __n the order n of the Bernoulli number.
   *  @return  The Bernoulli number of order n.
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __bernoulli_series(unsigned int __n)
    {
      constexpr unsigned long _S_num_bern_tab = 12;
      constexpr _Tp
      _S_bernoulli_2n[_S_num_bern_tab]
      {
	 _Tp{1ULL},
	 _Tp{1ULL}             / _Tp{6ULL},
	-_Tp{1ULL}             / _Tp{30ULL},
	 _Tp{1ULL}             / _Tp{42ULL},
	-_Tp{1ULL}             / _Tp{30ULL},
	 _Tp{5ULL}             / _Tp{66ULL},
	-_Tp{691ULL}           / _Tp{2730ULL},
	 _Tp{7ULL}             / _Tp{6ULL},
	-_Tp{3617ULL}          / _Tp{510ULL},
	 _Tp{43867ULL}         / _Tp{798ULL},
	-_Tp{174611ULL}        / _Tp{330ULL},
	 _Tp{854513ULL}        / _Tp{138ULL}
      };
      constexpr _Tp _S_2pi = __gnu_cxx::__const_2_pi(_Tp{});

      if (__n == 0)
	return _Tp{1};
      else if (__n == 1)
	return -_Tp{1} / _Tp{2};
      // Take care of the rest of the odd ones.
      else if (__n % 2 == 1)
	return _Tp{0};
      // Take care of some small evens that are painful for the series.
      else if (__n / 2 < _S_num_bern_tab)
	return _S_bernoulli_2n[__n / 2];
      else
	{
	  auto __fact = _Tp{1};
	  if ((__n / 2) % 2 == 0)
	    __fact *= -_Tp{1};
	  for (unsigned int __k = 1; __k <= __n; ++__k)
	    __fact *= __k / _S_2pi;
	  __fact *= _Tp{2};

	 // Riemann zeta function minus-1 for even integer argument.
	  auto __sum = _Tp{0};
	  for (unsigned int __i = 2; __i < 1000; ++__i)
	    {
	      auto __term = std::pow(_Tp(__i), -_Tp(__n));
	      __sum += __term;
	      if (__term < __gnu_cxx::__epsilon<_Tp>() * __sum)
		break;
	    }

	  return __fact + __fact * __sum;
	}
    }

  /**
   *   @brief This returns Bernoulli number @f$ B_n @f$.
   *
   *   @param __n the order n of the Bernoulli number.
   *   @return  The Bernoulli number of order n.
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __bernoulli(unsigned int __n)
    { return __bernoulli_series<_Tp>(__n); }

  /**
   * @brief This returns Bernoulli number @f$ B_2n @f$ at even integer
   * arguments @f$ 2n @f$.
   *
   * @param __n the half-order n of the Bernoulli number.
   * @return  The Bernoulli number of order 2n.
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __bernoulli_2n(unsigned int __n)
    { return __bernoulli_series<_Tp>(2 * __n); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_BERNOULLI_TCC

