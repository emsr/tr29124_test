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

/** @file bits/sf_stirling.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// C++ Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_SF_STIRLING_TCC
#define _GLIBCXX_BITS_SF_STIRLING_TCC 1

#pragma GCC system_header

#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the Stirling number of the second kind from lookup
   * or by series expansion.
   *
   * The series is:
   * @f[
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *
   *   \sigma_n^{(m)} = \stirling{n}{m}
   *      = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
   *
   * @todo Find a way to predict the maximum Stirling number for a type.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_series(unsigned int __n, unsigned int __m)
    {
      if (__m > _S_num_factorials<_Tp>)
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      auto __lf1 = __log_factorial<_Tp>(__k);
	      auto __lf2 = __log_factorial<_Tp>(__m - __k);
	      _S2 += ((__m - __k) & 1 ? _Tp{-1} : _Tp{1})
		   * std::exp(__n * std::log(__k) - __lf1 - __lf2);
	    }
	  return _S2;
	}
      else
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      _S2 += ((__m - __k) & 1 ? _Tp{-1} : _Tp{1})
		   * std::pow(__k, __n)
		   / __factorial<_Tp>(__k)
		   / __factorial<_Tp>(__m - __k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(_S2);
	}
    }

  /**
   * Return the Stirling number of the second kind by recursion.
   * The recursion is
   * @f[
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *   \stirling{n}{m} = m \stirling{n-1}{m} + \stirling{n-1}{m-1}
   * @f]
   * with starting values
   * @f[
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *   \stirling{0}{0\rightarrow m} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *   \stirling{0\rightarrow n}{0} = {1, 0, 0, ..., 0}
   * @f]
   *
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> __sigold(__m + 1), __signew(__m + 1);
	  __sigold[1] = _Tp{1};
	  if (__n == 1)
	    return __sigold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      __signew[1] = __sigold[1];
	      for (auto __im = 2u; __im <= __m; ++__im)
		__signew[__im] = __im * __sigold[__im] + __sigold[__im - 1];
	      std::swap(__sigold, __signew);
	    }
	  return __signew[__m];
	}
    }

  /**
   * Return the Stirling number of the second kind from lookup
   * or by series expansion.
   *
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * @todo Find asymptotic solutions for Stirling numbers of the second kind.
   * @todo Develop an iterator model for Stirling numbers of the second kind.
   */
  template<typename _Tp>
    _Tp
    __stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 10)
	return _Tp{0};
      else
	return __stirling_2_recur<_Tp>(__n, __m);
    }

  /**
   * Return the Stirling number of the second kind.
   *
   * @todo Look into asymptotic solutions.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return -std::numeric_limits<_Tp>::infinity();
      else if (__m == __n)
	return _Tp{0};
      else if (__m == 0 && __n >= 1)
	return -std::numeric_limits<_Tp>::infinity();
      else
	return std::log(__stirling_2<_Tp>(__n, __m));
    }

  /**
   * Return the Stirling number of the first kind by series expansion.
   * N.B. This seems to be a total disaster.
   */
  template<typename _Tp>
    _Tp
    __stirling_1_series(unsigned int __n, unsigned int __m)
    {
      using __gnu_cxx::__parity;
      if (2 * __n - __m > _S_num_factorials<_Tp> / 2)
	{
	  auto _S1 = _Tp{0};
	  for (auto __k = 0u; __k <= __n - __m; ++__k)
	    {
	      const auto __nmpk = __n - __m + __k;
	      const auto __nmmk = __n - __m - __k;
	      const auto __lbc1 = __log_binomial<_Tp>(__n - 1 + __k, __nmpk);
	      const auto __slbc1 = __log_binomial_sign<_Tp>(__n - 1 + __k, __nmpk);
	      const auto __lbc2 = __log_binomial<_Tp>(2 * __n - __m, __nmmk);
	      const auto __slbc2 = __log_binomial_sign<_Tp>(2 * __n - __m, __nmmk);
	      _S1 += __parity<_Tp>(__k) * __slbc1 * __slbc2
		   * std::exp(__lbc1 + __lbc2 + __log_stirling_2<_Tp>(__nmpk, __k));
	    }
	  return _S1;
	}
      else
	{
	  auto _S1 = _Tp{0};
	  for (auto __k = 0u; __k <= __n - __m; ++__k)
	    {
	      const auto __nmpk = __n - __m + __k;
	      const auto __nmmk = __n - __m - __k;
	      _S1 += __parity<_Tp>(__k)
		   * __binomial<_Tp>(__n - 1 + __k, __nmpk)
		   * __binomial<_Tp>(2 * __n - __m, __nmmk)
		   * __stirling_2<_Tp>(__nmpk, __k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(_S1);
	}
    }

  /**
   * Return the Stirling number of the first kind by recursion.
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __stirling_1_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> _Sold(__m + 1), _Snew(__m + 1);
	  _Sold[1] = _Tp{1};
	  if (__n == 1)
	    return _Sold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Snew[__im] = _Sold[__im - 1] - __in * _Sold[__im];
	      std::swap(_Sold, _Snew);
	    }
	  return _Snew[__m];
	}
    }

  /**
   * Return the Stirling number of the first kind.
   *
   * The Stirling numbers of the first kind are the coefficients of
   * the Pocchammer polynomials:
   * @f[
   *   (x)_n = \sum_{k=0}^{n} S_n^{(k)} x^k
   * @f]
   *
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   *
   * @todo Find asymptotic solutions for the Stirling numbers of the first kind.
   * @todo Develop an iterator model for Stirling numbers of the first kind.
   */
  template<typename _Tp>
    _Tp
    __stirling_1(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 1)
	return _Tp{0};
      else
        return __stirling_1_recur<_Tp>(__n, __m);
    }

  /**
   * Return the logarithm of the absolute value of Stirling number
   * of the first kind.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_1(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return -std::numeric_limits<_Tp>::infinity();
      else if (__m == __n)
	return _Tp{0};
      else if (__m == 0 && __n >= 1)
	return -std::numeric_limits<_Tp>::infinity();
      else
	return std::log(std::abs(__stirling_1<_Tp>(__n, __m)));
    }

  /**
   * Return the sign of the exponent of the logarithm of the Stirling number
   * of the first kind.
   */
  template<typename _Tp>
    inline _Tp
    __log_stirling_1_sign(unsigned int __n, unsigned int __m)
    { return (__n + __m) & 1 ? _Tp{-1} : _Tp{+1}; }


_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_STIRLING_TCC

