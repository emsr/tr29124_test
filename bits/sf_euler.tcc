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

/** @file bits/sf_euler.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// C++ Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_SF_EULER_TCC
#define _GLIBCXX_BITS_SF_EULER_TCC 1

#pragma GCC system_header

#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the Euler number from lookup or by series expansion.
   *
   * The Euler numbers are given by the recursive sum:
   * @f[
   *   E_n = B_n(1) = B_n
   * @f]
   * where @f$ E_0 = 1 @f$, @f$ E_1 = 0 @f$, @f$ E_2 = -1 @f$
   *
   * @todo Find a way to predict the maximum Euler number for a type.
   */
  template<typename _Tp>
    _Tp
    __euler_series(unsigned int __n)
    {
      static constexpr std::size_t _S_len = 22;
      static constexpr _Tp
      _S_num[_S_len]
      {
	 1ll, 0,
	-1ll, 0ll,
	 5ll, 0ll,
	-61ll, 0ll,
	 1385ll, 0ll,
	-50521ll, 0ll,
	 2702765ll, 0ll,
	-199360981ll, 0ll,
	 19391512145ll, 0ll,
	-2404879675441ll, 0ll,
	 370371188237525ll, 0ll,
	//-69348874393137901ll, 0ll,
      };

      if (__n == 0)
	return _Tp{1};
      else if (__n & 1)
	return _Tp{0};
      else if (__n == 2)
        return _Tp{-1};
      else if (__n < _S_len)
	return _S_num[__n];
      else
	{
	  std::vector<_Tp> _En(__n + 1);
	  _En[0] = _Tp{1};
	  _En[1] = _Tp{0};
	  _En[2] = _Tp{-1};

	  for (auto __i = 3u; __i <= __n; ++__i)
	    {
	      _En[__i] = 0;

	      if (__i % 2 == 0)
		{
		  for (auto __j = 2u; __j <= __i; __j += 2u)
		    _En[__i] -= __binomial<_Tp>(__i, __j)
			      * _En[__i - __j];
		}
	    }
	  return _En[__n];
	}
    }

  /**
   * @brief This returns Euler number @f$ E_n @f$.
   *
   * @param __n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename _Tp>
    inline _Tp
    __euler(unsigned int __n)
    { return __euler_series<_Tp>(__n); }

  /**
   * Return the Euler polynomial @f$ E_n(x) @f$ of order n at argument x.
   *
   * The derivative is proportional to the previous polynomial:
   * @f[
   *   E_n'(x) = n E_{n-1}(x)
   * @f]
   *
   * @f[
   *   E_n(1/2) = \frac{E_n}{2^n}, \mbox{ where } E_n
   *             \mbox{ is the n-th Euler number.}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __euler(unsigned int __n, _Tp __x)
    {
      if (__isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __bx1 = __bernoulli(__n + 1, __x );
	  auto __bx2 = __bernoulli(__n + 1, _Tp{0.5L} * __x );

	  auto _E_n = _Tp{2} * (__bx1 - __bx2 * std::pow(_Tp{2}, _Tp(__n + 1)))
		    / _Tp(__n + 1);

	  return _E_n;
	}
    }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers of the first kind are defined by recursion:
   * @f[
   *   \newcommand{\eulerian}[2]{\genfrac{\langle}{\rangle}{0pt}{0}{#1}{#2}}
   *
   *   \eulerian{n}{m} = (n-m)\eulerian{n-1}{m-1} + (m+1)\eulerian{n-1}{m}
   *   \mbox{ for } n > 0
   * @f]
   * Note that @f$ A(n,m) @f$ is a common older notation.
   */
  template<typename _Tp>
    _Tp
    __eulerian_1_recur(unsigned int __n, unsigned int __m)
    {
      if (__m == 0)
	return _Tp{1};
      else if (__m >= __n)
	return _Tp{0};
      else if (__m == __n - 1)
	return _Tp{1};
      else if (__n - __m - 1 < __m) // Symmetry.
	return __eulerian_1_recur<_Tp>(__n, __n - __m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<_Tp> _Aold(__m + 1), _Anew(__m + 1);
	  _Aold[0] = _Tp{1};
	  _Anew[0] = _Tp{1};
	  _Anew[1] = _Tp{1};
	  for (auto __in = 3u; __in <= __n; ++__in)
	    {
	      std::swap(_Aold, _Anew);
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Anew[__im] = (__in - __im) * _Aold[__im - 1]
			    + (__im + 1) * _Aold[__im];
	    }
	  return _Anew[__m];
	}
    }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers of the first kind are defined by recursion:
   * @f[
   *    \genfrac\langle\rangle{0pt}{0}{n}{m}
   *       = (n-m)\genfrac\langle\rangle{0pt}{0}{n-1}{m-1}
   *       + (m+1)\genfrac\langle\rangle{0pt}{0}{n-1}{m}
   *    \mbox{ for } n > 0
   * @f]
   * Note that @f$ A(n,m) @f$ is a common older notation.
   */
  template<typename _Tp>
    inline _Tp
    __eulerian_1(unsigned int __n, unsigned int __m)
    { return __eulerian_1_recur<_Tp>(__n, __m); }

  /**
   * Return the Eulerian number of the second kind by recursion.
   * The recursion is:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __eulerian_2_recur(unsigned int __n, unsigned int __m)
    {
      if (__m == 0)
	return _Tp{1};
      else if (__m >= __n)
	return _Tp{0};
      else if (__n == 0)
	return _Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<_Tp> _Aold(__m + 1), _Anew(__m + 1);
	  _Aold[0] = _Tp{1};
	  _Anew[0] = _Tp{1};
	  _Anew[1] = _Tp{2};
	  for (auto __in = 3u; __in <= __n; ++__in)
	    {
	      std::swap(_Aold, _Anew);
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Anew[__im] = (2 * __in - __im - 1) * _Aold[__im - 1]
			    + (__im + 1) * _Aold[__im];
	    }
	  return _Anew[__m];
	}
    }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    inline _Tp
    __eulerian_2(unsigned int __n, unsigned int __m)
    { return __eulerian_2_recur<_Tp>(__n, __m); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_EULER_TCC
